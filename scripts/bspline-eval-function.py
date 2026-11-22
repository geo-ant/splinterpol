# little bit of code generation for evaluating the piecewise functional
# form of a B-Spline of arbitrary (but given) degree.
# The relevant section form the Briand & Monasse paper is 2.3

import sympy as sp
import sys
from fractions import Fraction as frac
from math import comb as n_over_k
from dataclasses import dataclass
from sympy.printing.rust import RustCodePrinter

# compute C^(n)_{k,j} from eq. (32) in the paper. I don't care about
# computational efficiency here, which is why I'm fine not using recursion
# for caching precomputed results.
def compute_cnkj(n:int,k:int,j:int) -> int:
    if j == n:
        sum = 0
        # i = 0,...,k
        for i in range(0,k+1):
            sum += (-1)**i * n_over_k(n+1,i)
        return sum
    else:
        sum = 0
        # i = 0,...,k-1
        for i in range(0,k):
            sum += (-1)**i * n_over_k(n+1,i) * (k-i)**(n-j)

        sum *= n_over_k(n,j)
        return sum

# compute D^(n)_{ntilde,j} from eq 36. Again, I don't care about
# computational efficiency here
def compute_dnkj(n:int, j:int) -> frac:
    sum = 0
    ntilde = n//2
    # i = 0,...,ntilde
    for i in range(0, ntilde+1):
        sum += (-1)**(i+j) * n_over_k(n+1,i) * (n+1-2*i)**(n-j)

    sum *= n_over_k(n,j)
    return frac(sum, 2**(n-j))

@dataclass
class PolyPiece:
    upper_limit: frac
    poly: sp.Poly

# returns the polynomial pieces such that they can be evaluated
# in order
# if |x| < upper_limit[0]
#   //special for this segment! no substitution necessary
#   y = x
# 
#   poly_piece[0](y)
# ...
# if |x| < upper_limit[k]
#   // for all the other pieces we must substitute!
#   y = (n+1)/2 -k -x
#   polypiece[k](y)
def compute_piecewise_betan(n: int) -> list[PolyPiece]:
    if n < 1:
        raise ValueError("This construction requires degree n >= 1")

    ntilde = n // 2
    y = sp.Symbol("y")
    pieces: list[PolyPiece] = []

    # this was a bit annoying for me to understand. The important
    # section is 2.3 in the B&M paper. The beta^(n) (x) can be expressed
    # with eq 33 as a sum of y^j, where y = (n+1)-x-k with 0 < y <= 1.
    # The impotant thing is that the index k is a function of x, meaning
    # an interval x. We can see e.g. that x < (n+1)/2-k. And k = 0,...,ntilde
    # where ntilde = floor(n/2). That means the largest k constrains
    # the first interval, which is why we're traversing in order.
    # So the first interval (k=ntilde) x < (n+1)/2 - ntilde is special
    # because we use the expansion given in eq. 37, which uses coefficient
    # D^(n),{ntilde,j} as defined in 36. In this equation, no substitution
    # of x occurs.
    # For all the other intervals we use eq 33 which is a function of y,
    # not of x, which is why we need to substitude y = (n+1)/2-x-k
    # and use the coefficients C^(n)_{k,j} as defined in (32).
    # Very annoying to understand but well, that's that I guess. The
    # important thing is that the piecewise stuff comes from the fact
    # that the x intervals define the k that we need to use for the
    # coefficients.
    for k in range(ntilde, -1, -1):
        upper_limit = frac(n + 1, 2) - k

        if k == ntilde:
            coefficients = [
                compute_dnkj(n, j) for j in range(n + 1)
            ]
        else:
            coefficients = [
                compute_cnkj(n, k, j) for j in range(n + 1)
           ]

        expression = sum(
            sp.Rational(c.numerator, c.denominator) * y**j
            for j, c in enumerate(coefficients)
        )

        poly = sp.Poly(expression, y, domain=sp.QQ)
        pieces.append(PolyPiece(upper_limit, poly))

    return pieces        

# I realized I was doing something suboptimal, let me explain: I'm calculating
# the weights by first finding the "left hand" coordinate as t = x - (x-R).ceil()
# (this is true for x and y coordinates, I'll just use x here). Where R is the
# radius of the B-Spline support R = (n+1)/2. Then I'm filling
# an array of weights by repeatedly calling the eval function with:
# [eval(t-0), eval(t-1), ..., eval(t-(NW-1))] where NW is the number of weights
# for that B-Spline (which is n+1, where n is the degree of the B-Spline).
# The eval function uses branching logic depending on where it's argument falls,
# to find the right piece of the polynomial. Filling the array like I did by
# repeatedly calling eval is a bit stupid because I have to repeatedly pass
# the branching logic, unless the compiler can understand that there is only
# ever one branch for any particular t-m branch. Why is that?
#
# The reason is that by defintion t = x - (x-R).ceil(). The term (x-R).ceil()
# is in the range [x-R,x-R+1]. That means t \in [R-1,R], i.e. t \in [(n-1)/2,(n+1)/2].
# That means we know the possible range for each t - m expression, m = 0,...,NW-1
# which should give us a way to pick the correct branch at compile time for
# each eval(t-m) and then generate a function that returns the array of weights
# directly PROVIDED THAT IT WAS INDEED CALLED WITH the correct t argument.
#
# This function takes the piecewise beta(n) and re
def compute_weight_array(betan: list[PolyPiece], n: int) -> list[sp.Expr]:
    """Return expressions for [eval(t), eval(t-1), ..., eval(t-n)]."""
    if n < 1 or len(betan) != n // 2 + 1:
        raise ValueError("wrong number of polynomial pieces")

    t = sp.Symbol("t")
    t_lower = frac(n - 1, 2)
    t_upper = frac(n + 1, 2)
    weights = []

    for m in range(n + 1):
        arg = t - m
        cur_lower = t_lower - m
        cur_upper = t_upper - m

        # cur_lower < arg <= cur_upper, but eval branches on |arg| (the
        # absolute value) and here |cur_lower| might be larger than |cur_upper|
        # so we have to pick the correct piece that covers the range
        abs_upper = max(abs(cur_lower), abs(cur_upper))
        piece = next(p for p in betan if abs_upper <= p.upper_limit)

        if cur_lower < 0 < cur_upper:
            # arg can have either sign. This only happens in the central
            # slot of an even-degree spline. Its polynomial contains only
            # even powers, so P(arg) == P(abs(arg)).
            local_y = arg
        else:
            # arg has a known sign throughout this slot, possibly reaching
            # zero at an endpoint. We can resolve abs(arg) during code generation.
            abs_arg = -arg if cur_upper <= 0 else arg

            # The central polynomial takes abs(arg) directly; the outer
            # polynomials take the distance from their upper limit.
            if piece is betan[0]:
                local_y = abs_arg
            else:
                local_y = sp.Rational(piece.upper_limit) - abs_arg

        # y is the polynomial running variable, see above. That's why
        # we use local_y above
        y = piece.poly.gens[0]
        polynomial = sp.horner(piece.poly.as_expr(), y)
        weights.append(polynomial.subs(y, local_y))

    return weights

# heavily LLM assisted with the implementation here, my python
# is too crappy for pretty printing rust code...
class PolynomialRustPrinter(RustCodePrinter):
    def _cast_to_float(self, expr):
        # Avoid the stock printer's implicit casts and lost parentheses.
        return expr

    def _print_Integer(self, expr, _type=False):
        # Polynomial coefficients must be floating-point literals.
        return f"{int(expr)}.0"

    def _print_Pow(self, expr):
        # powi's exponent must remain an integer.
        if not expr.exp.is_Integer:
            raise ValueError("Expected an integer power")
        return f"({self._print(expr.base)}).powi({int(expr.exp)})"

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)

        if q & (q - 1):
            raise ValueError("Expected a power-of-two denominator")

        # q = 2**places, so p/q = (p * 5**places) / 10**places.
        places = q.bit_length() - 1

        if places == 0:
            return f"{p}.0"

        digits = str(abs(p) * 5**places).zfill(places + 1)
        sign = "-" if p < 0 else ""

        return f"{sign}{digits[:-places]}.{digits[-places:]}"

# llm assistance due to crappy python skills
def gen_rust_weight_impl(weights: list[sp.Expr]) -> str:
    printer = PolynomialRustPrinter()
    nw = len(weights)

    # One CSE pass over the entire array, so different weights can share terms.
    temporaries, outputs = sp.cse(weights, symbols=sp.numbered_symbols("tmp"))

    radius = sp.Rational(nw, 2)
    lower = printer.doprint(radius - 1)
    upper = printer.doprint(radius)

    lines = [
        "#[inline(always)]",
        f"pub fn weights(t: f64) -> [f64; {nw}] {{",
        f"    debug_assert!(t >= {lower} && t <= {upper});",
    ]

    for name, expr in temporaries:
        lines.append(f"    let {name}: f64 = {printer.doprint(expr)};")

    lines.append("    [")
    for j, expr in enumerate(outputs):
        lines.append(f"        {printer.doprint(expr)}, // w[{j}]")
    lines.extend(["    ]", "}"])

    return "\n".join(lines)

# also lots of llm assistance because my python is too crappy here
def gen_rust_eval_impl(poly_pieces: list[PolyPiece]) -> str:
    printer = PolynomialRustPrinter()
    y = sp.Symbol("y")

    lines = [
        "pub fn beta(x: f64) -> f64 {",
        "    let x = x.abs();",
    ]

    for index, piece in enumerate(poly_pieces):
        bound = piece.upper_limit
        limit = printer.doprint(
            sp.Rational(bound.numerator, bound.denominator)
        )

        lines.append(f"    if x < {limit} {{")

        if index == 0:
            lines.append("        let y = x;")
        else:
            lines.append(f"        let y = ({limit}) - x;")

        # Ensure the polynomial variable is named y.
        expr = piece.poly.as_expr().subs(piece.poly.gens[0], y)
        horner = sp.horner(expr, y)

        temporaries, expressions = sp.cse(
            horner, symbols=sp.numbered_symbols("tmp")
        )

        for name, value in temporaries:
            lines.append(
                f"        let {name}: f64 = {printer.doprint(value)};"
            )

        lines.append(
            f"        return {printer.doprint(expressions[0])};"
        )
        lines.append("    }")

    lines.extend([
        "    0.0",
        "}",
    ])

    return "\n".join(lines)

if len(sys.argv) != 2:
    sys.exit(f"Usage: python {sys.argv[0]} <spline degree>")

n = int(sys.argv[1])
poly_pieces = compute_piecewise_betan(n);
# for reference how the eval function would look. We're now using the
# weight function directly which just returns [eval(t),eval(t-1),...,eval(t-(NWEIGHTS-1))]
# directly.
# print(gen_rust_eval_impl(poly_pieces))
print(gen_rust_weight_impl(compute_weight_array(poly_pieces,n)))

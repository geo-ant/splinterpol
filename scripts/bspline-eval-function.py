# little bit of code generation for evaluating the piecewise functional
# form of a b-spline of arbitrary (but given) degree.
# The relevant section form the Briand & Monasse paper is 2.3

import sympy as sp
from fractions import Fraction as frac
from math import comb as n_over_k
from dataclasses import dataclass

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

def gen_rust_impl(poly_pieces: list[PolyPiece]):

    for index, piece in enumerate(poly_pieces):
        print(f"if x < {piece.upper_limit} {{")
        if index == 0:
            # no substitution for this bit
            print("let y = x;")
        else:
            print(f"let y = {piece.upper_limit} - x")

    
        print("TODO return poly eval")

    print("return 0.;")


poly_pieces = compute_piecewise_betan(3)

print(poly_pieces)

print("foo")


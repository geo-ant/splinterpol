# a little utility script that can calculate the roots of the \tilde{B}^(n)(z)
# polynomial in Briand and Monasse, appendix A.1 with arbitraty precision.
# I could already observe (through failing tests) that the poles of the
# BSpline of degree 4 are not correct within float64 precision in their reference
# implementation.

from fractions import Fraction as frac
import sympy as sp

# algorithm 9 in the B&M paper
def compute_z_transform_coefficients(n: int) -> list[frac]:

    if n < 2:
        raise Exception("spline degree must be >= 2")
    
    ntilde = n // 2
    # these are b^(n) (k) in the paper and
    # d^(n) (k)
    bn = [frac(0)] * (ntilde+2)
    dn = [frac(0)] * (ntilde+1)

    bn[0] = frac(1)
    dn[0] = frac(1,2)
        
    for m in range(1,n):
        mtilde = m // 2
        # these are the b^(m) (k) and d^(m) (k)
        # temporary arrays in the B&M paper
        bm = [frac(0)] * (ntilde+1)
        dm = [frac(0)] * (ntilde+1)
        for k in range(0,mtilde+1):
            # the last line is a hack to get d^(m)_{-1} = d^(m)_0, see
            # appendix A.1 in the paper
            bm[k] = ( frac(m+1,2) + k ) * dn[k]+ ( frac(m+1,2) - k)* dn[k-1 if k > 0 else 0]
            bm[k] /= m
            dm[k] = ( frac(m+2,2) + k) * bn[k+1] + (frac(m,2) - k) * bn[k]
            dm[k] /= m

        # reassign the temporaries to the original arrays
        for k in range(0,mtilde+1):
            bn[k] = bm[k]
            dn[k] = dm[k]

    for k in range(0,ntilde+1):
        bn[k] = ( frac(n+1,2) + k) * dn[k] + ( frac(n+1,2) - k)* dn[k-1 if k > 0 else 0]
        bn[k] /= n

    # now mirror and expand the coefficients such that
    # coeff[k] is the coefficient for power z^k in the
    # polynomial \tilde{B}^(n)(z) = coeff[k]*z^k
    # in the B&M article appendix A.1

    coeff = [0] * (n+1 if n % 2 == 0 else n)

    for k in range(0,ntilde+1):
        coeff[k] = bn[ntilde-k]
    for k in range(ntilde+1,2*ntilde+1):
        coeff[k] = coeff[2*ntilde-k]
        
    return coeff

# build a sympy polynomial with rational coefficients from coeff[k]
# such that y(x) = \sum_k coeff[k] x^k
def make_polynomial(coeff: list[Fraction]) -> sp.Poly:
    z = sp.Symbol("z")

    expression = sum(
        sp.Rational(c.numerator, c.denominator) * z**k
        for k, c in enumerate(coeff)
    )

    return sp.Poly(expression, z, domain=sp.QQ)


# absolute overkill
SIGNIFICANT_DIGITS = 45
MAX_DEGREE = 13

print("#![allow(clippy::excessive_precision)]")
print("/// *** AUTO GENERATED. DO NOT EDIT ***")
print(f"/// Poles of B-Splines with degrees 2,...,{MAX_DEGREE}")
print(f"/// auto-generated with scripts/calculate_poles.py\n\n")
for degree in range(2,MAX_DEGREE+1):
    coeff = compute_z_transform_coefficients(degree)
    poly = make_polynomial(coeff)
    croots = sp.real_roots(poly)
    # use many more significant digits than are needed or useful for float64
    # just because we can... It should even be enough for float128 if that
    # becomes available in the far future.
    roots = [root.evalf(SIGNIFICANT_DIGITS) for root in croots]
    poles = [root for root in roots if -1 < root < 0]
    print(f"/// poles for B-Spline of degree {degree}")
    print(f"pub const BSPLINE{degree}_POLES: [f64; {len(poles)}] = [")
    for pole in poles:
        print(f"    {pole:e},")
    print("];\n")

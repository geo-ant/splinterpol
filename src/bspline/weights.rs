use crate::{
    bspline::{BSpline3, BSpline4},
    numerics::FusedMuladd,
};
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float};

// NOTE: why the stupid constant for the number of weights,
// even though the number of weights is defined by each degree only?
// See the comments on the other trait for a rant.
pub trait BSplineWeights<const NWEIGHTS: usize>: Copy + Default {
    type Float: Float;
    fn eval(self, coord: Self::Float) -> Self::Float;
    /// this must be exactly `(NWEIGHTS as Float)/2.0`
    fn radius() -> Self::Float;
}

impl<F> BSplineWeights<4> for BSpline3<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    type Float = F;

    #[inline(always)]
    fn eval(self, x: F) -> F {
        // generated this function with `python scripts/bspline-eval-function.py 3`
        // and then did optimizations and cleanups by hand
        let x = x.abs();
        if x < F::ONE {
            let y = x;
            return (y)
                .powi(2)
                .fmuladd(y.fmuladd(3.0.as_(), -6.0.as_()), 4.0.as_());
        }
        if x < 2.0.as_() {
            let y = (2.0.as_()) - x;
            return (y).powi(3);
        }
        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        Self::Float::ONE + Self::Float::ONE
    }
}

impl<F> BSplineWeights<5> for BSpline4<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    type Float = F;

    #[inline]
    fn eval(self, x: F) -> F {
        // NOTE(geo): also created with the script and then had an llm apply the muladd.
        // The muladd refactoring is very tedious but schematic because we
        // can always go from inside to outside. The one thing we should consider
        // is whether it makes sense to split the calculation into multiple
        // dependency chains for larger calculations, which I don't think
        // is neccessarily possible if we use the Horner rewrite for calculating
        // the polynomials, which is what Briand & Monasse also do... gotta think
        // about that.
        let x = x.abs();

        if x < 0.5.as_() {
            let x2 = x * x;
            return x2.fmuladd(x2.fmuladd(6.0.as_(), (-15.0).as_()), 14.375.as_());
        }

        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            let p = y.fmuladd((-4.0).as_(), 4.0.as_());
            let p = y.fmuladd(p, 6.0.as_());
            let p = y.fmuladd(p, 4.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 2.5.as_() {
            let y = 2.5.as_() - x;
            let y2 = y * y;
            return y2 * y2;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        2.5.as_()
    }
}

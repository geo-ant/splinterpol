use crate::{
    bspline::{BSpline, BSpline3, BSpline4, BSpline5},
    numerics::FusedMuladd,
};
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float};

// NOTE: why the stupid constant for the number of weights,
// even though the number of weights is defined by each degree only?
// See the comments on the other trait for a rant.
/// **Note** The number `NWEIGHTS` is the number of elements in the weight
/// vector (`interpn` calls it the _footprint_) used for evaluating the values
/// of the spline around the grid points for interpolation. The number of
/// weights is exactly `n+1` for a BSpline with polynomial degreen `n`.
pub trait BSplineWeights<const NWEIGHTS: usize>: BSpline {
    fn eval(self, coord: Self::Float) -> Self::Float;
    /// this must be exactly `(NWEIGHTS as Float)/2.0`
    fn radius() -> Self::Float;
}

impl<F> BSplineWeights<4> for BSpline3<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
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

impl<F> BSplineWeights<6> for BSpline5<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    fn eval(self, x: Self::Float) -> Self::Float {
        // from the script, then manually refined, then let an llm handle
        // the tedious work of writing this as muladds
        let x = x.abs();

        if x < F::ONE {
            let x2 = x * x;

            // 66 + x² * (-60 + x² * (30 - 10*x))
            let p = x.fmuladd((-10.0).as_(), 30.0.as_());
            let p = x2.fmuladd(p, (-60.0).as_());
            return x2.fmuladd(p, 66.0.as_());
        }

        if x < 2.0.as_() {
            let y = 2.0.as_() - x;

            // 1 + y * (5 + y * (10 + y * (10 + y * (5 - 5*y))))
            let p = y.fmuladd((-5.0).as_(), 5.0.as_());
            let p = y.fmuladd(p, 10.0.as_());
            let p = y.fmuladd(p, 10.0.as_());
            let p = y.fmuladd(p, 5.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 3.0.as_() {
            let y = 3.0.as_() - x;
            let y2 = y * y;
            let y4 = y2 * y2;
            return y4 * y;
        }

        F::ZERO
    }

    fn radius() -> Self::Float {
        3.0.as_()
    }
}

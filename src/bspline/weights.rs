use crate::bspline::{BSpline3, BSpline4};
use num_traits::{ConstOne, ConstZero, Float};

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
    F: Float + ConstOne + ConstZero,
{
    type Float = F;

    #[inline(always)]
    fn eval(self, x: F) -> F {
        let x = x.abs();
        let n2 = F::ONE + F::ONE;
        let n3 = n2 + F::ONE;
        let n4 = n2 + n2;
        let n6 = n4 + n2;
        if x < F::ONE {
            n4 + (-n6 + n3 * x) * x * x
        } else if x < n2 {
            (n2 - x).powi(3)
        } else {
            F::ZERO
        }
    }

    #[inline]
    fn radius() -> Self::Float {
        Self::Float::ONE + Self::Float::ONE
    }
}

impl BSplineWeights<5> for BSpline4<f64> {
    type Float = f64;

    #[inline]
    fn eval(self, coord: f64) -> f64 {
        let mut x = coord.abs();
        if x <= 0.5 {
            x *= x;
            14.375 + (-15. + 6. * x) * x
        } else if x < 1.5 {
            x = 1.5 - x;
            1. + (4. + (6. + (4. - 4. * x) * x) * x) * x
        } else if x < 2.5 {
            x = 2.5 - x;
            x *= x;
            x * x
        } else {
            return 0.;
        }
    }

    #[inline]
    fn radius() -> Self::Float {
        2.5
    }
}

impl BSplineWeights<5> for BSpline4<f32> {
    type Float = f32;

    #[inline]
    fn eval(self, coord: f32) -> f32 {
        let mut x = coord.abs();
        if x <= 0.5 {
            x *= x;
            14.375 + (-15. + 6. * x) * x
        } else if x < 1.5 {
            x = 1.5 - x;
            1. + (4. + (6. + (4. - 4. * x) * x) * x) * x
        } else if x < 2.5 {
            x = 2.5 - x;
            x *= x;
            x * x
        } else {
            return 0.;
        }
    }

    #[inline]
    fn radius() -> Self::Float {
        2.5
    }
}

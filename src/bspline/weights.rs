use crate::bspline::BSpline3;
use num_traits::{ConstOne, ConstZero, Float};

// NOTE: why the stupid constant for the number of weights,
// even though the number of weights is defined by each degree only?
// See the comments on the other trait for a rant.
pub trait BSplineWeights<const NWEIGHTS: usize>: Copy + Default {
    type Float: Float;
    fn eval(self, coord: Self::Float) -> Self::Float;
    /// this must be exactly float(NWEIGTS)/2.0
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

use num_traits::Float;

use crate::bspline::BSpline3;

// NOTE: why the stupid constant for the number of weights,
// even though the number of weights is defined by each degree only?
// See the comments on the other trait for a rant.
pub trait BSplineWeights<const NWEIGHTS: usize>: Copy + Default {
    type Float: Float;
    fn weights(self, coord: Self::Float) -> [Self::Float; NWEIGHTS];
}

impl<F> BSplineWeights<4> for BSpline3<F>
where
    F: Float,
{
    type Float = F;

    #[inline(always)]
    fn weights(self, coord: F) -> [Self::Float; 4] {
        todo!()
    }
}

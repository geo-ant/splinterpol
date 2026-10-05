#![allow(clippy::excessive_precision)]
use core::marker::PhantomData;
use num_traits::Float;

// NOTE: All B-Splines in this crate implement both the trait
// for prefiltering and the trait for weight calculation. I've
// split this into two parts because the traits are completely
// independently used and the current limitations to const
// generic programming in Rust mean that I'd have to stick
// two constants (number of poles and number of weights) on the
// trait that could both be calculated from the degree of the
// B-Spline, but I can't do that in current rust. It would
// just make the interpolation and coefficient calculation function
// signatures more annoying, so I've decided to simply split this.

/// logic and traits that are used for coefficient calculation
/// (aka prefiltering)
mod prefilter;

/// logic and traits that are used for interpolation calculations
/// (aka calculating the weights).
mod weights;

pub use prefilter::BSplinePrefilter;
pub use weights::BSplineWeights;

/// alias for the BSpline of order 3
pub type CubicBSpline<T> = BSpline3<T>;

/// BSpline of order 3
#[derive(Debug, Clone, Copy)]
pub struct BSpline3<T: Float> {
    phantom: PhantomData<T>,
}

impl<T: Float> Default for BSpline3<T> {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Float> BSpline3<T> {
    #[inline]
    pub fn new() -> Self {
        Self {
            phantom: Default::default(),
        }
    }
}

/// BSpline of order 10
#[derive(Default, Debug, Clone, Copy)]
pub struct BSpline10<T> {
    phantom: PhantomData<T>,
}

/// B-spline of order 11
#[derive(Default, Debug, Clone, Copy)]
pub struct BSpline11<T> {
    phantom: PhantomData<T>,
}

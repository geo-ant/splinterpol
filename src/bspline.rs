#![allow(clippy::excessive_precision)]
use core::marker::PhantomData;
use num_traits::Float;

mod prefilter;

pub use prefilter::BSplinePrefilter;

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

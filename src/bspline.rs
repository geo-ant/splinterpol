#![allow(clippy::excessive_precision)]
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
// In a perfect world, we would only have one trait with an
// _associated constant_ DEGREE and then NPOLES = floor(DEGREE/2) and
// NWEIGHTS = DEGREE+1 could be the array dimensions of the returned
// poles and weights, but we can't use DEGREE/2 or DEGREE+1 as array
// boundaries. This is doubly true if DEGREE was an associated constant,
// cause todays Rust wouldn't even allow using DEGREE as the number of
// elements in an array. But even if we made the BSpline trait generic over
// DEGREE, then we still can't use DEGREE/2 or DEGREE+1 as array counts.

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

/// A supertrait for a B-Spline describing the floating point type it takes
/// and it's degree. Due to Rust's limitations around const generic programming,
/// it isn't useful on its own. The useful numerical implementations are in
/// [`BSplinePrefilter`] and [`BSplineWeights`].
pub trait BSpline: Copy + Default {
    /// The floating point type that a particular B-Spline operates on.
    type Float: num_traits::Float;
    /// The polynomial degree of the B-Spline. Some papers call this the _order_,
    /// but other papers use that term differently. Degree is unambiguous.
    const DEGREE: usize;
}

/// BSpline of degree 2
#[derive(Debug, Clone, Copy)]
pub struct BSpline2<T: Float> {
    // NOTE(geo) due to limited rust generics, we have to put the number
    // of poles in here, which must match the number of poles in the prefilter
    // trait. This match will be enforced at compile time, so it's fine.
    poles: [T; 1],
}

/// BSpline of degree 3
#[derive(Debug, Clone, Copy)]
pub struct BSpline3<T: Float> {
    poles: [T; 1],
}

/// BSpline of degree 4
#[derive(Debug, Clone, Copy)]
pub struct BSpline4<T> {
    poles: [T; 2],
}

/// BSpline of degree 5
#[derive(Debug, Clone, Copy)]
pub struct BSpline5<T> {
    poles: [T; 2],
}

/// BSpline of degree 6
#[derive(Debug, Clone, Copy)]
pub struct BSpline6<T> {
    poles: [T; 3],
}

/// BSpline of degree 7
#[derive(Debug, Clone, Copy)]
pub struct BSpline7<T> {
    poles: [T; 3],
}

/// BSpline of degree 8
#[derive(Debug, Clone, Copy)]
pub struct BSpline8<T> {
    poles: [T; 4],
}

/// BSpline of degree 9
#[derive(Debug, Clone, Copy)]
pub struct BSpline9<T> {
    poles: [T; 4],
}

/// BSpline of degree 10
#[derive(Debug, Clone, Copy)]
pub struct BSpline10<T> {
    poles: [T; 5],
}

/// B-spline of degree 11
#[derive(Debug, Clone, Copy)]
pub struct BSpline11<T> {
    poles: [T; 5],
}

macro_rules! impl_bspline {
    ($spline:ident, $degree:literal) => {
        impl<T> BSpline for $spline<T>
        where
            T: num_traits::Float + 'static,
            f64: num_traits::AsPrimitive<T>,
        {
            type Float = T;
            const DEGREE: usize = $degree;
        }
    };
}

impl_bspline!(BSpline2, 2);
impl_bspline!(BSpline3, 3);
impl_bspline!(BSpline4, 4);
impl_bspline!(BSpline5, 5);
impl_bspline!(BSpline6, 6);
impl_bspline!(BSpline7, 7);
impl_bspline!(BSpline8, 8);
impl_bspline!(BSpline9, 9);
impl_bspline!(BSpline10, 10);
impl_bspline!(BSpline11, 11);

#![allow(clippy::excessive_precision)]

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
// elements in an array. But even if we made the B-Spline trait generic over
// DEGREE, then we still can't use DEGREE/2 or DEGREE+1 as array counts.

/// logic and traits that are used for coefficient calculation
/// (aka prefiltering)
mod prefilter;

/// logic and traits that are used for interpolation calculations
/// (aka calculating the weights).
mod weights;

pub use prefilter::BSplinePrefilter;
pub use weights::BSplineWeights;

use crate::sealed::Sealed;

/// alias for the B-Spline of order 3
pub type CubicBSpline<T> = BSpline3<T>;

/// A supertrait for a B-Spline describing the floating point type it takes
/// and it's degree. Due to Rust's limitations around const generic programming,
/// it isn't useful on its own. The useful numerical implementations are in
/// [`BSplinePrefilter`] and [`BSplineWeights`].
pub trait BSpline: Copy + Default + Sealed {
    /// The floating point type that a particular B-Spline operates on.
    type Float: num_traits::Float;
    /// The polynomial degree of the B-Spline. Some papers call this the _order_,
    /// but other papers use that term differently. Degree is unambiguous.
    const DEGREE: usize;
}

#[derive(Debug, Clone, Copy)]
/// Unifying underlying type for B-Splines of all degrees with floating point type
/// `T`. Only B-Splines of degree 2 up to 12 are actually usable in this library.
///
/// Due to limitations in Rust's const generics, we have to make this type
/// generic over the polynomial degree `DEGREE` and the number of poles `NPOLES`,
/// when in fact for any B-Spline `NPOLES = floor(DEGREE/2)`. That's why all
/// usable B-Splines in this library have aliases indicating their degree,
/// such as `BSpline2`, `BSpline3`, ... , `BSpline12`. Still, this type _might_
/// come in handy in interfaces accepting B-Splines generically.
pub struct BasisSpline<T, const DEGREE: usize, const NPOLES: usize> {
    poles: [T; NPOLES],
}

/// B-Spline of polynomial degree 2
pub type BSpline2<T> = BasisSpline<T, 2, 1>;
/// B-Spline of polynomial degree 3
pub type BSpline3<T> = BasisSpline<T, 3, 1>;
/// B-Spline of polynomial degree 4
pub type BSpline4<T> = BasisSpline<T, 4, 2>;
/// B-Spline of polynomial degree 5
pub type BSpline5<T> = BasisSpline<T, 5, 2>;
/// B-Spline of polynomial degree 6
pub type BSpline6<T> = BasisSpline<T, 6, 3>;
/// B-Spline of polynomial degree 7
pub type BSpline7<T> = BasisSpline<T, 7, 3>;
/// B-Spline of polynomial degree 8
pub type BSpline8<T> = BasisSpline<T, 8, 4>;
/// B-Spline of polynomial degree 9
pub type BSpline9<T> = BasisSpline<T, 9, 4>;
/// B-Spline of polynomial degree 10
pub type BSpline10<T> = BasisSpline<T, 10, 5>;
/// B-Spline of polynomial degree 11
pub type BSpline11<T> = BasisSpline<T, 11, 5>;
/// B-Spline of polynomial degree 12
pub type BSpline12<T> = BasisSpline<T, 12, 6>;

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

        impl<F> crate::sealed::Sealed for $spline<F>
        where
            F: num_traits::Float + Copy + 'static,
            f64: num_traits::AsPrimitive<F>,
        {
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
impl_bspline!(BSpline12, 12);

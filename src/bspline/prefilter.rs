use crate::{
    bspline::{
        BSpline2, BSpline3, BSpline4, BSpline5, BSpline6, BSpline7, BSpline8, BSpline9, BSpline10,
        BSpline11,
        prefilter::poles::{
            BSPLINE2_POLES, BSPLINE3_POLES, BSPLINE4_POLES, BSPLINE5_POLES, BSPLINE6_POLES,
            BSPLINE7_POLES, BSPLINE8_POLES, BSPLINE9_POLES, BSPLINE10_POLES, BSPLINE11_POLES,
        },
    },
    util::ArrayConvert,
};
use num_traits::{AsPrimitive, ConstOne, Float};

mod poles;

// fuck I hate everything about this NPOLES trait bound but generic parameters
// may not be used in constant expressions yet and this is so fucking stupid
// I am losing my shit.
pub trait BSplinePrefilter<const NPOLES: usize>: Copy + Default {
    /// Floating point type for computations
    type Float: num_traits::Float;
    // we can't fucking do this because
    // "generic parameters may not be used in const expressions"
    // so we would have to stick the fucking number of poles into the
    // fucking trait for this to work, which is fucking stupid because
    // one implementor will only ever have one number of poles and
    // it wouldn't make any sense to implement this trait twice.
    // const NPOLES: usize;
    // const POLES: [T; Self::NPOLES];
    // Also I hate that I can't make this an associated constant, but
    // the Float trait cannot work with this as trait methods in general
    // cannot be const.
    /// The poles of the BSpline. The number of for a BSpline of degree
    /// `n` is `floor(n/2)`, see [B&M 2018](https://www.ipol.im/pub/art/2018/221/).
    fn poles(&self) -> &[Self::Float; NPOLES];
    /// normalization constant from eq. (75) in the [B&M 2018](https://www.ipol.im/pub/art/2018/221/)
    /// paper, which is `1` for odd degree `n` and `2^n` otherwise.
    fn normalization() -> Self::Float;
}

/// utility macro to implement the constructors and the prefilter trait for
/// a given B-Spline
macro_rules! impl_bspline {
    (
        $spline:ident,
        npoles = $npoles:literal,
        poles = $poles:ident,
        norm = $normalization:literal
    ) => {
        impl<F: num_traits::Float> $spline<F>
        where
            F: Copy + 'static,
            f64: AsPrimitive<F>,
        {
            #[inline]
            /// create a new instance
            pub fn new() -> Self {
                Self {
                    poles: ArrayConvert::new($poles).as_(),
                }
            }
        }

        /// same as `Self::new` constructor
        impl<T> Default for $spline<T>
        where
            f64: AsPrimitive<T>,
            T: Float + Copy + 'static,
        {
            #[inline]
            fn default() -> Self {
                Self::new()
            }
        }

        impl<F> BSplinePrefilter<$npoles> for $spline<F>
        where
            F: num_traits::Float + Copy + 'static + ConstOne,
            f64: AsPrimitive<F>,
        {
            type Float = F;
            #[inline]
            fn normalization() -> Self::Float {
                $normalization.as_()
            }

            #[inline]
            fn poles(&self) -> &[Self::Float; $npoles] {
                &self.poles
            }
        }
    };
}

impl_bspline!(BSpline2, npoles = 1, poles = BSPLINE2_POLES, norm = 4.0);
impl_bspline!(BSpline3, npoles = 1, poles = BSPLINE3_POLES, norm = 1.0);
impl_bspline!(BSpline4, npoles = 2, poles = BSPLINE4_POLES, norm = 16.0);
impl_bspline!(BSpline5, npoles = 2, poles = BSPLINE5_POLES, norm = 1.0);
impl_bspline!(BSpline6, npoles = 3, poles = BSPLINE6_POLES, norm = 64.0);
impl_bspline!(BSpline7, npoles = 3, poles = BSPLINE7_POLES, norm = 1.0);
impl_bspline!(BSpline8, npoles = 4, poles = BSPLINE8_POLES, norm = 256.);
impl_bspline!(BSpline9, npoles = 4, poles = BSPLINE9_POLES, norm = 1.0);
impl_bspline!(
    BSpline10,
    npoles = 5,
    poles = BSPLINE10_POLES,
    norm = 1024.0
);
impl_bspline!(BSpline11, npoles = 5, poles = BSPLINE11_POLES, norm = 1.0);

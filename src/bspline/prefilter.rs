use crate::bspline::{
    BSpline3, BSpline4, BSpline10, BSpline11,
    prefilter::poles::{BSPLINE3_POLES, BSPLINE4_POLES, BSPLINE10_POLES, BSPLINE11_POLES},
};
use num_traits::{AsPrimitive, ConstOne};

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
    /// the poles of the BSpline
    fn poles() -> [Self::Float; NPOLES];
    /// normalization constant TODO gamma something, but not `n!` (fakultaet)
    fn normalization() -> Self::Float;
}

impl<F> BSplinePrefilter<1> for BSpline3<F>
where
    F: num_traits::Float + Copy + 'static + ConstOne,
    f64: AsPrimitive<F>,
{
    type Float = F;
    #[inline]
    fn normalization() -> Self::Float {
        F::ONE
    }

    #[inline]
    fn poles() -> [Self::Float; 1] {
        core::array::from_fn(|idx| BSPLINE3_POLES[idx].as_())
    }
}

impl<F> BSplinePrefilter<2> for BSpline4<F>
where
    F: num_traits::Float + Copy + 'static + ConstOne,
    f64: AsPrimitive<F>,
{
    type Float = F;

    #[inline]
    fn poles() -> [Self::Float; 2] {
        const { assert!(BSPLINE4_POLES.len() == 2) }
        core::array::from_fn(|idx| BSPLINE4_POLES[idx].as_())
    }

    #[inline]
    fn normalization() -> Self::Float {
        16.0.as_()
    }
}

impl<F> BSplinePrefilter<5> for BSpline10<F>
where
    F: num_traits::Float + Copy + 'static + ConstOne,
    f64: AsPrimitive<F>,
{
    type Float = F;

    #[inline]
    fn poles() -> [Self::Float; 5] {
        const { assert!(BSPLINE10_POLES.len() == 5) }
        core::array::from_fn(|idx| BSPLINE10_POLES[idx].as_())
    }

    #[inline]
    fn normalization() -> Self::Float {
        1024.0.as_()
    }
}

impl<F> BSplinePrefilter<5> for BSpline11<F>
where
    F: num_traits::Float + Copy + 'static + ConstOne,
    f64: AsPrimitive<F>,
{
    type Float = F;

    #[inline]
    fn poles() -> [Self::Float; 5] {
        const { assert!(BSPLINE11_POLES.len() == 5) }
        core::array::from_fn(|idx| BSPLINE11_POLES[idx].as_())
    }

    #[inline]
    fn normalization() -> Self::Float {
        F::ONE
    }
}

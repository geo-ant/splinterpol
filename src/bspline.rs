use core::marker::PhantomData;
use num_traits::Float;

// fuck I hate everything about this NPOLES trait bound but generic parameters
// may not be used in constant expressions yet and this is so fucking stupid
// I am losing my shit.
pub trait BSpline<T: Float, const NPOLES: usize> {
    // number of poles of the b-spline
    // const NPOLES: usize;
    /// normalization constant required for coefficient calculations
    const NORMALIZATION: T;
    // we can't fucking do this because
    // "generic parameters may not be used in const expressions"
    // so we would have to stick the fucking number of poles into the
    // fucking trait for this to work, which is fucking stupid because
    // one implementor will only ever have one number of poles and
    // it wouldn't make any sense to implement this trait twice.
    // const NPOLES: usize;
    // const POLES: [T; Self::NPOLES];
    /// returns the slice of poles, which must be of length npoles
    fn poles(&self) -> &[T; NPOLES];
}

/// alias for the spline of order 3
pub type CubicBSpline<T> = BSpline3<T>;

/// B-spline of order 3
#[derive(Default, Debug, Clone, Copy)]
pub struct BSpline3<T: Float> {
    phantom: PhantomData<T>,
}

// TODO: make macro for this for the other splines
impl BSpline3<f32> {
    const POLES: [f32; 1] = [-0.2679491924311227 as _];
}

impl BSpline<f32, 1> for BSpline3<f32> {
    const NORMALIZATION: f32 = 1.;

    #[inline(always)]
    fn poles(&self) -> &[f32; 1] {
        &Self::POLES
    }
}

impl BSpline3<f64> {
    const POLES: [f64; 1] = [-0.2679491924311227 as _];
}

impl BSpline<f64, 1> for BSpline3<f64> {
    const NORMALIZATION: f64 = 1.;

    #[inline(always)]
    fn poles(&self) -> &[f64; 1] {
        &Self::POLES
    }
}

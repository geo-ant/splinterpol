use crate::bspline::{BSpline3, BSpline4, BSpline10, BSpline11};

// fuck I hate everything about this NPOLES trait bound but generic parameters
// may not be used in constant expressions yet and this is so fucking stupid
// I am losing my shit.
pub trait BSplinePrefilter<const NPOLES: usize>: Copy + Default {
    /// Floating point type for computations
    type Float: num_traits::Float;
    /// normalization constant for coefficient computation
    const NORMALIZATION: Self::Float;
    // we can't fucking do this because
    // "generic parameters may not be used in const expressions"
    // so we would have to stick the fucking number of poles into the
    // fucking trait for this to work, which is fucking stupid because
    // one implementor will only ever have one number of poles and
    // it wouldn't make any sense to implement this trait twice.
    // const NPOLES: usize;
    // const POLES: [T; Self::NPOLES];
    /// the poles of the BSpline
    const POLES: [Self::Float; NPOLES];
}

impl BSplinePrefilter<1> for BSpline3<f32> {
    type Float = f32;
    const NORMALIZATION: f32 = 1.;
    const POLES: [f32; 1] = [-0.2679491924311227 as _];
}

impl BSplinePrefilter<1> for BSpline3<f64> {
    type Float = f64;
    const NORMALIZATION: f64 = 1.;
    const POLES: [f64; 1] = [-0.2679491924311227 as _];
}

impl BSplinePrefilter<2> for BSpline4<f64> {
    type Float = f64;
    const NORMALIZATION: Self::Float = 16.;
    // NOTE(geo): the poles in Briand & Monasse are _not correct_ for f64
    // precision.
    const POLES: [f64; 2] = [-0.36134122590022016, -0.013725429297339121];
}

impl BSplinePrefilter<2> for BSpline4<f32> {
    type Float = f32;
    const NORMALIZATION: Self::Float = 16.;
    const POLES: [f32; 2] = [-0.36134122590022016, -0.013725429297339121];
}

impl BSplinePrefilter<5> for BSpline10<f64> {
    type Float = f64;
    const NORMALIZATION: f64 = 1024.;

    const POLES: [f64; 5] = [
        -6.365506639694650e-1,
        -2.381827983775487e-1,
        -6.572703322831758e-2,
        -7.528194675547741e-3,
        -1.698276282327549e-5,
    ];
}

impl BSplinePrefilter<5> for BSpline10<f32> {
    type Float = f32;
    const NORMALIZATION: f32 = 1024.;

    const POLES: [f32; 5] = [
        -6.365506639694650e-1,
        -2.381827983775487e-1,
        -6.572703322831758e-2,
        -7.528194675547741e-3,
        -1.698276282327549e-5,
    ];
}

impl BSplinePrefilter<5> for BSpline11<f64> {
    type Float = f64;
    const NORMALIZATION: f64 = 1.;
    const POLES: [f64; 5] = [
        -6.612660689007345e-1,
        -2.721803492947859e-1,
        -8.975959979371331e-2,
        -1.666962736623466e-2,
        -5.105575344465021e-4,
    ];
}

impl BSplinePrefilter<5> for BSpline11<f32> {
    type Float = f32;
    const NORMALIZATION: f32 = 1.;
    const POLES: [f32; 5] = [
        -6.612660689007345e-1,
        -2.721803492947859e-1,
        -8.975959979371331e-2,
        -1.666962736623466e-2,
        -5.105575344465021e-4,
    ];
}

/// utitilty module for the macro-generated (I mean in future they will
/// be macro-generated) tests for the poles.
// TODO: actually, remove this once I use macros, because this won't matter
// anymore...

// utility function to make sure that f32 and f64 arrays have the same contents
// at compile time
// TODO: make this part of the macro that I'll eventually use to define the
// bspline structs! That's why I'm writing it like this already
#[cfg(test)]
mod comptime_test_bspline3 {
    use super::BSpline3;
    use super::BSplinePrefilter;
    use crate::test_utils::array_f32_f64_eq;
    const ASSRT: () = assert!(array_f32_f64_eq(
        BSpline3::<f32>::POLES,
        BSpline3::<f64>::POLES
    ));

    #[test]
    fn bspline3_poles() {
        let _assertion = ASSRT;
    }
}

#[cfg(test)]
mod comptime_test_bspline11 {
    use super::BSpline11;
    use super::BSplinePrefilter;
    use crate::test_utils::array_f32_f64_eq;
    const ASSRT: () = assert!(array_f32_f64_eq(
        BSpline11::<f32>::POLES,
        BSpline11::<f64>::POLES
    ));

    #[test]
    fn bspline11_poles() {
        let _assertion = ASSRT;
    }
}

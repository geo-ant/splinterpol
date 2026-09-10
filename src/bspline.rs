#![allow(clippy::excessive_precision)]
use core::marker::PhantomData;
use num_traits::Float;

// fuck I hate everything about this NPOLES trait bound but generic parameters
// may not be used in constant expressions yet and this is so fucking stupid
// I am losing my shit.
pub trait BSpline<const NPOLES: usize> {
    /// Floating point type with which the calculations for this spline type
    /// can be performed
    type Float: num_traits::Float;
    // it would be nice if we could also do this for the associated constants, but...
    // see below.
    // number of poles of the b-spline
    // const NPOLES: usize;
    // normalization constant required for coefficient calculations
    const NORMALIZATION: Self::Float;
    // we can't fucking do this because
    // "generic parameters may not be used in const expressions"
    // so we would have to stick the fucking number of poles into the
    // fucking trait for this to work, which is fucking stupid because
    // one implementor will only ever have one number of poles and
    // it wouldn't make any sense to implement this trait twice.
    // const NPOLES: usize;
    // const POLES: [T; Self::NPOLES];
    /// returns the slice of poles, which must be of length npoles
    const POLES: [Self::Float; NPOLES];
}

/// alias for the spline of order 3
pub type CubicBSpline<T> = BSpline3<T>;

/// B-spline of order 3
#[derive(Default, Debug, Clone, Copy)]
pub struct BSpline3<T: Float> {
    phantom: PhantomData<T>,
}

impl BSpline<1> for BSpline3<f32> {
    type Float = f32;
    const NORMALIZATION: f32 = 1.;
    const POLES: [f32; 1] = [-0.2679491924311227 as _];
}

impl BSpline3<f64> {}

impl BSpline<1> for BSpline3<f64> {
    type Float = f64;
    const NORMALIZATION: f64 = 1.;
    const POLES: [f64; 1] = [-0.2679491924311227 as _];
}

/// BSpline of order 10
#[derive(Default, Debug, Clone, Copy)]
pub struct BSpline10<T> {
    phantom: PhantomData<T>,
}

impl BSpline<5> for BSpline10<f64> {
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

impl BSpline<5> for BSpline10<f32> {
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

/// B-spline of order 11
#[derive(Default, Debug, Clone, Copy)]
pub struct BSpline11<T> {
    phantom: PhantomData<T>,
}

impl BSpline<5> for BSpline11<f64> {
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

impl BSpline<5> for BSpline11<f32> {
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
#[cfg(test)]
mod test_util {
    /// compare a f32 and a f64 array for equality
    pub const fn array_eq<const N: usize>(lhs: [f32; N], rhs: [f64; N]) -> bool {
        // can't use for loops in const functions... cause of reasons I guess.
        let mut j = 0;
        loop {
            if j >= N {
                break;
            }
            if lhs[j] != rhs[j] as f32 {
                return false;
            }
            j += 1;
        }
        true
    }
}

// utility function to make sure that f32 and f64 arrays have the same contents
// at compile time
// TODO: make this part of the macro that I'll eventually use to define the
// bspline structs! That's why I'm writing it like this already
#[cfg(test)]
mod comptime_test_bspline3 {
    use super::BSpline;
    use super::BSpline3;
    use super::test_util::array_eq;
    const ASSRT: () = assert!(array_eq(BSpline3::<f32>::POLES, BSpline3::<f64>::POLES));

    #[test]
    fn bspline3_poles() {
        let _assertion = ASSRT;
    }
}

#[cfg(test)]
mod comptime_test_bspline11 {
    use super::BSpline;
    use super::BSpline11;
    use super::test_util::array_eq;
    const ASSRT: () = assert!(array_eq(BSpline11::<f32>::POLES, BSpline11::<f64>::POLES));

    #[test]
    fn bspline11_poles() {
        let _assertion = ASSRT;
    }
}

use std::{fs::File, io::Write};

use crate::{
    bspline::*,
    inplace::{
        TransmittableBoundaryExtension,
        builder::{RowMajor, SplineInterpolator2d},
        interpolation::test::synthetic_image::SynthImage,
    },
    test_utils::{read_f32le_asf64, read_f64le},
};
use float_cmp::{ApproxEq, FloatMargin, approx_eq};
use num_traits::{AsPrimitive, Float, FromPrimitive, float::TotalOrder};
use test_case::test_case;

mod synthetic_image;

// little hack macro to make the spline test function "generic" over the Bspline
// without actually pulling in another generic test dependency for now, because
// I don't even know if that would play nicely with the #[test_case] stuff
// that I have going on. Also I have two generic parameters...
//
// thinking of this, I could have just made this a much simpler impl by
// not parsing the generic function, but just write that once in this
// test module. And then the other functions just call that generic impl...
// Then I wouldn't even need recursion. But now I'm quite proud of this
// macro so that stays in...
macro_rules! gen_bspline_test {
    // top level parsing, tailored very specifically to the bspline test
    (
        // #[splines($Spline:ident<$Scalar:ty> $(,$rest:ty)* $(,)?)]
        #[splines($($splines:tt)*)]
        $(
            #[$meta:meta]
        )+
        fn $fn:ident<B,const NP :usize, const NW: usize> ($($arg:ident : $type:ty),+ $(,)?)
        // this is my lazy as fuck way of not parsing the where clause and still
        // only ever having to write this once... I know, right?
        where ... $body:block
    ) => {
        // the generic function which will be called by all recursively
        // expanded functions internally
        fn $fn<B,const NP :usize, const NW: usize> ($($arg : $type),+)
        where B: BSpline<Float = f64> + BSplineWeights<NW> + BSplinePrefilter<NP>, $body

        gen_bspline_test! {
            @internal,
            splines = {$($splines)*},
            meta = {$(#[$meta])+},
            fn = $fn,
            args = {$($arg: $type),+}

        }
    };

    // make a concrete implementation that calls the generic implementation
    // for the first b-spline in the list.
    (
        @internal,
        splines = {$(,)? $spline:ident<$Scalar:ty> $($rest:tt)*},
        meta = {$(#[$meta:meta])+},
        fn = $fn:ident,
        args = {$($arg:ident : $type:ty),+} $(,)?
    ) => {

        paste::paste! {
            $(
                #[$meta]
            )+
            #[track_caller]
            fn [<$fn _ $spline:lower>]($($arg:$type),+) {
                $fn::<$spline<$Scalar>,_,_>($($arg),+)
            }
        }

        gen_bspline_test! {
            @internal,
            splines = {$($rest)*},
            meta = {$(#[$meta])+},
            fn = $fn,
            args = {$($arg: $type),+}

        }
    };

    // end of the recursion, just do nothing.
    (
        @internal,
        splines = {$(,)?},
        meta = {$(#[$meta:meta])+},
        fn = $fn:ident,
        args = {$($arg:ident : $type:ty),+} $(,)?
    ) => {};
}

// this test makes sure that the interpolation for all implemented spline
// orders _on the original grid coordinates_ stays within the given relative
// error bounds. We're only testing interpolation on the original grid coords,
// since we don't actually know the true interpolated values. I've seletected
// a range of image from the imagej sample images covering cells, brain mri,
// photos and such.
gen_bspline_test! {
#[splines(BSpline3<f64>,BSpline4<f64>,BSpline5<f64>,BSpline6<f64>,
    BSpline7<f64>,BSpline9<f64>,BSpline10<f64>,BSpline11<f64>)]
// NOTE(geo): see also below. The B&M paper specificies `epsilon` as
// a _relative bound_, see top of page 108. The error for each interpolated
// value is err(f_k) < epsilon max(|f|)
// However(!!), in the remarks on p. 137 at the end of appendix C.1,
// the paper also states:
// > For really small values of ε (e.g. ε < 10−12), the machine precision
// > should be considered. Thus  the experimental error may be higher in
// > this case.
//
// This is in line with the observation that for 1e-15, we get small, but
// noticeable floating point errors that violate the relative bound.
// ** this is a synthetic image generated in matlab that has proven pretty useful, so it stays in the suite.
#[test_case("test-assets/coefficients2d/image_512x1024.f64",512,1024,false,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/coefficients2d/image_512x1024.f64",512,1024,false,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/coefficients2d/image_512x1024.f64",512,1024,false,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/coefficients2d/image_512x1024.f64",512,1024,false,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/coefficients2d/image_512x1024.f64",512,1024,false,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/coefficients2d/image_512x1024.f64",512,1024,false,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
// *** blobs
#[test_case("test-assets/imagej/blobs_256x254.f32.raw",256,254,true,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/imagej/blobs_256x254.f32.raw",256,254,true,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/imagej/blobs_256x254.f32.raw",256,254,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/imagej/blobs_256x254.f32.raw",256,254,true,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/imagej/blobs_256x254.f32.raw",256,254,true,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/imagej/blobs_256x254.f32.raw",256,254,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
// *** boats
#[test_case("test-assets/imagej/boats_720x576.f32.raw",720,576,true,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/imagej/boats_720x576.f32.raw",720,576,true,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/imagej/boats_720x576.f32.raw",720,576,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/imagej/boats_720x576.f32.raw",720,576,true,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/imagej/boats_720x576.f32.raw",720,576,true,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/imagej/boats_720x576.f32.raw",720,576,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
// *** clown
#[test_case("test-assets/imagej/clown_320x200.f32.raw",320,200,true,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/imagej/clown_320x200.f32.raw",320,200,true,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/imagej/clown_320x200.f32.raw",320,200,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/imagej/clown_320x200.f32.raw",320,200,true,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/imagej/clown_320x200.f32.raw",320,200,true,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/imagej/clown_320x200.f32.raw",320,200,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
// *** fly-brain
#[test_case("test-assets/imagej/fly-brain_256x256.f32.raw",256,256,true,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/imagej/fly-brain_256x256.f32.raw",256,256,true,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/imagej/fly-brain_256x256.f32.raw",256,256,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/imagej/fly-brain_256x256.f32.raw",256,256,true,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/imagej/fly-brain_256x256.f32.raw",256,256,true,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/imagej/fly-brain_256x256.f32.raw",256,256,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
// *** gel
#[test_case("test-assets/imagej/gel_276x467.f32.raw",276,467,true,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/imagej/gel_276x467.f32.raw",276,467,true,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/imagej/gel_276x467.f32.raw",276,467,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/imagej/gel_276x467.f32.raw",276,467,true,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/imagej/gel_276x467.f32.raw",276,467,true,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/imagej/gel_276x467.f32.raw",276,467,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
// *** mri
#[test_case("test-assets/imagej/mri_186x226.f32.raw",186,226,true,TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case("test-assets/imagej/mri_186x226.f32.raw",186,226,true,TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case("test-assets/imagej/mri_186x226.f32.raw",186,226,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case("test-assets/imagej/mri_186x226.f32.raw",186,226,true,TransmittableBoundaryExtension::Periodic, 1e-13)]
#[test_case("test-assets/imagej/mri_186x226.f32.raw",186,226,true,TransmittableBoundaryExtension::Mirrored, 1e-13)]
#[test_case("test-assets/imagej/mri_186x226.f32.raw",186,226,true,TransmittableBoundaryExtension::WholeSymmetric, 1e-13)]
fn input_is_reproduced_within_relative_error_on_original_grid<B, const NP: usize, const NW: usize>(
    image_path: &str,
    width: usize,
    height: usize,
    is_f32: bool,
    ext: TransmittableBoundaryExtension,
    epsilon: f64,
) where ...
{
    // i know this is hacky...
    let original = if is_f32 {
        read_f32le_asf64(image_path)
    } else {
        read_f64le(image_path)}
    ;
    assert_eq!(original.len(), width*height, "image size probably wrong");
    let fmax = original
        .iter()
        .max_by(|a, b| f64::total_cmp(&a.abs(), &b.abs()))
        .copied()
        .unwrap()
        .abs();
    assert!(fmax.is_finite());
    let data = RowMajor {
        data: original.clone(),
        width,
        height,
    };
    let original = RowMajor {
        data: original,
        width,
        height,
    };
    let interpolator = SplineInterpolator2d::from_image(B::default(), data, ext, epsilon).unwrap();

    let mut abs_sum = 0.;

    let mut discrepancies = Vec::new();
    for r in 0..original.height {
        for c in 0..original.width {
            let actual = interpolator.interp_single(r as f64, c as f64).unwrap();
            let expected = *original.at(r, c);
            abs_sum += expected.abs();
            if !approx_eq!(
                f64,
                actual,
                expected,
                // NOTE: the B&M paper specifies that epsilon is a relative
                // bound so this should be the right bound to enforce here.
                epsilon = epsilon * fmax
            ) {
                discrepancies.push(((r,c),(expected,actual)));

            }
        }
    }

    assert!(discrepancies.is_empty(),"following errors found: {:?}",discrepancies);
    // just a sanity check
    assert!(
        abs_sum > f64::EPSILON,
        "interpolated images is likely all zeros... very suspicious"
    );
}
}

fn interpolation_for_synthetic_image<T, B, const NP: usize, const NW: usize>(
    bspline: B,
    width: usize,
    height: usize,
    ext: TransmittableBoundaryExtension,
    epsilon: T,
) where
    B: BSpline<Float = T> + BSplineWeights<NW> + BSplinePrefilter<NP>,
    T: Float
        + 'static
        + num_traits::ConstOne
        + num_traits::ConstZero
        + std::iter::Sum
        + std::fmt::Debug
        + ApproxEq
        + FromPrimitive
        + TotalOrder,
    f64: AsPrimitive<T>,
    T: AsPrimitive<<<T as ApproxEq>::Margin as FloatMargin>::F>,
    <<T as ApproxEq>::Margin as FloatMargin>::F: Copy,
{
    let synth_image = SynthImage::new(width, height);
    let image = synth_image.as_row_major();
    let fmax = image
        .data
        .iter()
        .map(|f| f.abs())
        .max_by(|a, b| T::total_cmp(&a, &b))
        .unwrap()
        .abs();

    let mut w = File::create(format!(
        "debug_image_{}x{}.double.raw",
        image.width, image.height
    ))
    .unwrap();
    unsafe {
        let as_bytes = core::slice::from_raw_parts(
            image.data.as_ptr() as *const u8,
            core::mem::size_of_val(image.data.as_slice()),
        );
        w.write_all(as_bytes).unwrap();
    }

    assert!(fmax.is_finite());
    let data = image.clone();
    let original = image;
    let interpolator = SplineInterpolator2d::from_image(bspline, data, ext, epsilon).unwrap();

    let mut abs_sum = T::ZERO;

    let mut discrepancies = Vec::new();
    for r in 0..original.height {
        for c in 0..original.width {
            let actual = interpolator
                .interp_single(T::from_usize(r).unwrap(), T::from_usize(c).unwrap())
                .unwrap();
            let expected = *original.at(r, c);
            abs_sum = abs_sum + expected.abs();
            if !approx_eq!(
                T,
                actual,
                expected,
                // NOTE: the B&M paper specifies that epsilon is a relative
                // bound so this should be the right bound to enforce here.
                epsilon = (epsilon * fmax).as_()
            ) {
                discrepancies.push(((r, c), (expected, actual)));
            }
        }
    }

    assert!(
        discrepancies.is_empty(),
        "following errors found: {:?}",
        discrepancies
    );
    // just a sanity check
    assert!(
        abs_sum > T::epsilon(),
        "interpolated images is likely all zeros... very suspicious"
    );
}

#[test]
#[track_caller]
fn test_foo() {
    interpolation_for_synthetic_image(
        BSpline3::<f64>::new(),
        1230,
        1123,
        TransmittableBoundaryExtension::Mirrored,
        1e-15,
    );
}

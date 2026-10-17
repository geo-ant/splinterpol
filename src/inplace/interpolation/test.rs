use crate::{
    bspline::{BSpline3, BSpline4, BSplinePrefilter, BSplineWeights},
    inplace::{
        TransmittableBoundaryExtension,
        builder::{RowMajor, SplineInterpolator2d},
    },
    test_utils::read_f64le,
};
use float_cmp::assert_approx_eq;
use test_case::test_case;

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
        where B: BSplineWeights<NW, Float = f64> + BSplinePrefilter<NP, Float = f64> $(,)? $body:block
    ) => {
        // the generic function which will be called by all recursively
        // expanded functions internally
        fn $fn<B,const NP :usize, const NW: usize> ($($arg : $type),+)
        where B: BSplineWeights<NW, Float = f64> + BSplinePrefilter<NP, Float = f64>, $body

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

gen_bspline_test! {
//TODO(geo): append the other splines here
#[splines(BSpline3<f64>,BSpline4<f64>)]
#[test_case(TransmittableBoundaryExtension::Periodic, 1e-8)]
#[test_case(TransmittableBoundaryExtension::Mirrored, 1e-8)]
#[test_case(TransmittableBoundaryExtension::WholeSymmetric, 1e-8)]
#[test_case(TransmittableBoundaryExtension::Periodic, 1e-15)]
#[test_case(TransmittableBoundaryExtension::Mirrored, 1e-15)]
#[test_case(TransmittableBoundaryExtension::WholeSymmetric, 1e-15)]
fn input_is_reproduced_within_relative_error_on_original_grid<B, const NP: usize, const NW: usize>(
    ext: TransmittableBoundaryExtension,
    epsilon: f64,
) where
    B: BSplineWeights<NW, Float = f64> + BSplinePrefilter<NP, Float = f64>,
{
    let original = read_f64le("test-assets/coefficients2d/image_512x1024.f64");
    let fmax = original
        .iter()
        .max_by(|a, b| f64::total_cmp(a, b))
        .unwrap()
        .abs();
    assert!(fmax.is_finite());
    let data = RowMajor {
        data: original.clone(),
        width: 512,
        height: 1024,
    };
    let original = RowMajor {
        data: original,
        width: 512,
        height: 1024,
    };
    let interpolator = SplineInterpolator2d::from_image(B::default(), data, ext, epsilon).unwrap();

    let mut abs_sum = 0.;
    for r in 0..original.height {
        for c in 0..original.width {
            let actual = interpolator.interp_single(r as f64, c as f64).unwrap();
            let expected = *original.at(r, c);
            abs_sum += expected.abs();
            assert_approx_eq!(
                f64,
                actual,
                expected,
                // NOTE: the B&M paper specifies that epsilon is a relative
                // bound so this should be the right bound to enforce here.
                epsilon = epsilon * fmax
            );
        }
    }
    // just a sanity check
    assert!(
        abs_sum > f64::EPSILON,
        "interpolated images is likely all zeros... very suspicious"
    );
}
}

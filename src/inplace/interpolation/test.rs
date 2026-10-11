use crate::{
    bspline::BSpline3,
    inplace::{
        TransmittableBoundaryExtension,
        builder::{RowMajor, SplineInterpolator2d},
    },
    test_utils::read_f64le,
};
use float_cmp::assert_approx_eq;
use test_case::test_case;

#[test_case(TransmittableBoundaryExtension::Periodic)]
#[test_case(TransmittableBoundaryExtension::Mirrored)]
#[test_case(TransmittableBoundaryExtension::WholeSymmetric)]
fn input_is_reproduced_within_relative_error_on_original_grid(ext: TransmittableBoundaryExtension) {
    let original = read_f64le("test-assets/coefficients2d/image_512x1024.f64");
    let fmax = *original.iter().max_by(|a, b| f64::total_cmp(a, b)).unwrap();
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
    let epsilon = 1e-8;
    let interpolator =
        SplineInterpolator2d::from_image(BSpline3::<f64>::new(), data, ext, epsilon).unwrap();

    for r in 0..original.height {
        for c in 0..original.width {
            let actual = interpolator.interp_single(r as f64, c as f64).unwrap();
            let expected = *original.at(r, c);
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
}

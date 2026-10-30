use crate::{
    RowMajor,
    bspline::BSpline3,
    inplace::{SplineInterpolator2d, TransmittableBoundaryExtension},
    test_utils::read_f32le_asf64,
};

use interpn::MultiBsplineRegular;
use test_case::test_case;
#[test_case(
    "test-assets/imagej/blobs_256x254.f32.raw",
    256,
    254,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case(
    "test-assets/imagej/blobs_256x254.f32.raw",
    256,
    254,
    TransmittableBoundaryExtension::Reflected
)]
#[test_case(
    "test-assets/imagej/blobs_256x254.f32.raw",
    256,
    254,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case(
    "test-assets/imagej/boats_720x576.f32.raw",
    720,
    576,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case(
    "test-assets/imagej/boats_720x576.f32.raw",
    720,
    576,
    TransmittableBoundaryExtension::Reflected
)]
#[test_case(
    "test-assets/imagej/boats_720x576.f32.raw",
    720,
    576,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case(
    "test-assets/imagej/clown_320x200.f32.raw",
    320,
    200,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case(
    "test-assets/imagej/clown_320x200.f32.raw",
    320,
    200,
    TransmittableBoundaryExtension::Reflected
)]
#[test_case(
    "test-assets/imagej/clown_320x200.f32.raw",
    320,
    200,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case(
    "test-assets/imagej/fly-brain_256x256.f32.raw",
    256,
    256,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case(
    "test-assets/imagej/fly-brain_256x256.f32.raw",
    256,
    256,
    TransmittableBoundaryExtension::Reflected
)]
#[test_case(
    "test-assets/imagej/fly-brain_256x256.f32.raw",
    256,
    256,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case(
    "test-assets/imagej/gel_276x467.f32.raw",
    276,
    467,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case(
    "test-assets/imagej/gel_276x467.f32.raw",
    276,
    467,
    TransmittableBoundaryExtension::Reflected
)]
#[test_case(
    "test-assets/imagej/gel_276x467.f32.raw",
    276,
    467,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case(
    "test-assets/imagej/mri_186x226.f32.raw",
    186,
    226,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case(
    "test-assets/imagej/mri_186x226.f32.raw",
    186,
    226,
    TransmittableBoundaryExtension::Reflected
)]
#[test_case(
    "test-assets/imagej/mri_186x226.f32.raw",
    186,
    226,
    TransmittableBoundaryExtension::WholeSymmetric
)]
// NOTE(geo): see the interpolation tests and documentation. Our epsilon
// is actually a relative bound, but here we want to compare the ABSOLUTE error
// against the error that interpn achieves on a range of real world test images.
// Strictly, our error guarantees will always be worse than the ones in interpn
// since we're using truncation in the boundaries to achieve a relative error
// ob epsilon. But in practice this might still produce competitive results
// at (potentially) higher speed. The reason that interpn has (theoretically)
// zero error is that it uses a tridiagonal solve to enforce the interpolation
// condition (spline matches the function at the grid points) exactly.
// There can only be floating point errors that make this inexact.
fn absolute_accuracy_on_gridpoints_vs_interpn(
    image_path: &str,
    width: usize,
    height: usize,
    ext: TransmittableBoundaryExtension,
) {
    let original = read_f32le_asf64(image_path);
    let epsilon = 1e-13;
    assert_eq!(original.len(), width * height, "image size probably wrong");
    let data = RowMajor {
        data: original.clone(),
        width: width.try_into().unwrap(),
        height: height.try_into().unwrap(),
    };
    let original = data.clone();
    let interpolator =
        SplineInterpolator2d::from_image(BSpline3::default(), data, ext, epsilon).unwrap();

    let dims = [width, height];
    let mut scratch_buffer =
        vec![0f64; MultiBsplineRegular::<f64, 2>::construction_scratch_len(dims)];
    let mut coeffs = vec![0f64; MultiBsplineRegular::<f64, 2>::coeff_storage_len(dims)];
    let interpn_interpolator = MultiBsplineRegular::<f64, 2>::from_values_with_workspace(
        dims,
        [0f64, 0f64],
        [1f64, 1f64],
        original.data.as_slice(),
        coeffs.as_mut_slice(),
        scratch_buffer.as_mut_slice(),
        false,
    )
    .unwrap();

    for r in 0..height {
        for c in 0..width {
            let interpn_error = (interpn_interpolator.interp_one([r as _, c as _]).unwrap()
                - original.at(r, c))
            .abs();
            let this_error =
                (interpolator.interp_single(r as _, c as _).unwrap() - original.at(r, c)).abs();
            assert!(
                // see above: interpn has no error on the grid points except
                // for floating point error.
                this_error <= interpn_error + 1e-10,
                "{this_error} <= {}",
                interpn_error
            );
        }
    }
}

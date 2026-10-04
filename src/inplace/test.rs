use float_cmp::assert_approx_eq;
use std::num::NonZero;

use crate::{
    TransmittableBoundaryExtension,
    bspline::{BSpline3, BSpline10},
    inplace::{
        apply_exponential_filter_in_place_1d, apply_exponential_filter_in_place_1d_iter,
        apply_exponential_filter_in_place_strided, coefficients2d, coefficients2d_blocked,
    },
    iter_utils::StridedSliceMut,
    test_utils::read_f64le,
};

/// reference results generated with my fork of the Briand&Monasse repository,
/// see: https://github.com/geo-ant/reference_bspline_interpolation
/// using the `efip` executable for the 1D data and the `prefilter2d`
/// executable for the 2d data

#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    1
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    2
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    7
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic,
    8
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric,
    8
)]
#[allow(unsafe_code)]
fn exponential_filter_reproduces_reference_results(
    data: &str,
    expected: &str,
    truncation_index: usize,
    alpha: f64,
    extension: TransmittableBoundaryExtension,
    blocksize: usize,
) {
    let orig = read_f64le(data);
    let mut data = orig.clone();
    let expected = read_f64le(expected);

    let apply_exponential_filter_in_place_1d_with_blocksize =
        |alpha: f64, data, truncation_index, extension| {
            // little bit of a hack to allow different block sizes as rutime args
            match blocksize {
                1 => apply_exponential_filter_in_place_1d::<_, 1>(
                    alpha,
                    data,
                    truncation_index,
                    extension,
                ),
                2 => apply_exponential_filter_in_place_1d::<_, 2>(
                    alpha,
                    data,
                    truncation_index,
                    extension,
                ),
                7 => apply_exponential_filter_in_place_1d::<_, 7>(
                    alpha,
                    data,
                    truncation_index,
                    extension,
                ),
                8 => apply_exponential_filter_in_place_1d::<_, 8>(
                    alpha,
                    data,
                    truncation_index,
                    extension,
                ),
                _ => {
                    panic!("unsupported dynamic blocksize")
                }
            }
        };

    apply_exponential_filter_in_place_1d_with_blocksize(
        alpha,
        data.as_mut_slice(),
        truncation_index,
        extension,
    )
    .unwrap();

    assert_approx_eq!(
        &[f64],
        expected.as_slice(),
        data.as_slice(),
        // NOTE: we're bit identical to the reference implementation
        // ulps = 0,
        // TODO reset to 0 again...
        epsilon = 5. * f64::EPSILON
    );

    // only run those other blocks in one test since they don't depend
    // on the blocksize itself
    if blocksize == 1 {
        let count = NonZero::new(orig.len()).unwrap();
        let mut data = orig.clone();
        unsafe {
            apply_exponential_filter_in_place_strided(
                alpha,
                data.as_mut_ptr(),
                const { NonZero::new(1).unwrap() },
                count,
                truncation_index,
                extension,
            )
            .unwrap();
        }

        assert_approx_eq!(
            &[f64],
            expected.as_slice(),
            data.as_slice(),
            // NOTE: we're bit identical to the reference implementation
            ulps = 0,
            epsilon = 0.
        );

        let mut data = orig;
        apply_exponential_filter_in_place_1d_iter(
            alpha,
            data.as_mut_slice(),
            truncation_index,
            extension,
        )
        .unwrap();

        assert_approx_eq!(
            &[f64],
            expected.as_slice(),
            data.as_slice(),
            // NOTE: we're bit identical to the reference implementation
            ulps = 0,
            epsilon = 0.
        );
    }
}

#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.hsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Mirrored
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.wsym.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.peri.alpha-0.283.f64",
    100,
    -0.283,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.hsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Mirrored
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.peri.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::Periodic
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n5.wsym.alpha-0.456.f64",
    5,
    -0.456,
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[allow(unsafe_code)]
fn exponential_filter_with_non_unit_stride(
    data: &str,
    expected: &str,
    truncation_index: usize,
    alpha: f64,
    extension: TransmittableBoundaryExtension,
) {
    const STRIDE: usize = 4;
    let orig = read_f64le(data);
    let count = NonZero::new(orig.len()).unwrap();
    let expected = read_f64le(expected);
    // we just pad the data with zeros between each element with a spacing of stride,
    // so abcd becomes a000b000c000d000
    let orig_strided: Vec<_> = orig
        .into_iter()
        .flat_map(|val| {
            let mut arr = [-123f64; STRIDE];
            arr[0] = val;
            arr
        })
        .collect();

    let mut data = orig_strided.clone();
    unsafe {
        apply_exponential_filter_in_place_strided(
            alpha,
            data.as_mut_ptr(),
            const { NonZero::new(STRIDE).unwrap() },
            count,
            truncation_index,
            extension,
        )
        .unwrap();
    }

    // now collapse the strided data again so we can compare
    let data: Vec<_> = data.chunks_exact(STRIDE).map(|chunk| chunk[0]).collect();

    assert_approx_eq!(
        &[f64],
        expected.as_slice(),
        data.as_slice(),
        // NOTE: we're bit identical to the reference implementation
        ulps = 0,
        epsilon = 0.
    );

    let mut data = orig_strided;
    apply_exponential_filter_in_place_1d_iter(
        alpha,
        StridedSliceMut::new(
            &mut data,
            crate::iter_utils::StridedAccessArgs {
                offset: 0,
                stride: NonZero::new(STRIDE).unwrap(),
                count,
            },
        )
        .unwrap(),
        truncation_index,
        extension,
    )
    .unwrap();

    // now collapse the strided data again so we can compare
    let data: Vec<_> = data.chunks_exact(STRIDE).map(|chunk| chunk[0]).collect();

    assert_approx_eq!(
        &[f64],
        expected.as_slice(),
        data.as_slice(),
        // NOTE: we're bit identical to the reference implementation
        ulps = 0,
        epsilon = 0.
    );
}

#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    1,
    1
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    1,
    1
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    1,
    1
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    8,
    8
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    8,
    8
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    8,
    8
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    6,
    7
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    6,
    7
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    6,
    7
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    7,
    6
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    7,
    6
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord3.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    7,
    6
)]
fn test_coefficients2d_bspline3(
    data: &str,
    expected: &str,
    width: usize,
    height: usize,
    ext: TransmittableBoundaryExtension,
    row_blocksize: usize,
    col_blocksize: usize,
) {
    let orig = read_f64le(data);
    let expected = read_f64le(expected);

    let epsilon = 1e-5;

    let mut image_data = orig.clone();

    // see the other tests: hacky, but useful way to allow runtime
    // blocksizes and transform this into compile time block sizes

    let coefficients2d_with_blocksize =
        |bspline, image_data: &mut [f64], width, height, epsilon, ext| match (
            row_blocksize,
            col_blocksize,
        ) {
            (1, 1) => coefficients2d(bspline, image_data, width, height, epsilon, ext),
            (6, 7) => coefficients2d_blocked::<_, _, 6, 7>(
                bspline, image_data, width, height, epsilon, ext,
            ),
            (7, 6) => coefficients2d_blocked::<_, _, 7, 6>(
                bspline, image_data, width, height, epsilon, ext,
            ),
            (8, 8) => coefficients2d_blocked::<_, _, 8, 8>(
                bspline, image_data, width, height, epsilon, ext,
            ),
            (r, c) => {
                panic!("unsupported blocksize: R={r}, C={c}")
            }
        };

    coefficients2d_with_blocksize(
        BSpline3::<f64>::default(),
        &mut image_data,
        width,
        height,
        epsilon,
        ext,
    )
    .unwrap();
    assert_approx_eq!(
        &[f64],
        &expected.as_slice(),
        &image_data.as_slice(),
        // tiny inaccuracies prevent us from being bit-identical, but this
        // is good enough for all intents and purposes...
        epsilon = 20. * f64::EPSILON
    );
}

#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    1,
    1
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    1,
    1
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    1,
    1
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    8,
    8
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    8,
    8
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    8,
    8
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    6,
    7
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    6,
    7
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    6,
    7
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.hsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Mirrored,
    7,
    6
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.peri.f64",
    512,
    1024,
    TransmittableBoundaryExtension::Periodic,
    7,
    6
)]
#[test_case::test_case(
    "test-assets/coefficients2d/image_512x1024.f64",
    "test-assets/coefficients2d/image_512x1024.ord10.eps1e-05.wsym.f64",
    512,
    1024,
    TransmittableBoundaryExtension::WholeSymmetric,
    7,
    6
)]
fn test_coefficients2d_bspline10(
    data: &str,
    expected: &str,
    width: usize,
    height: usize,
    ext: TransmittableBoundaryExtension,
    row_blocksize: usize,
    col_blocksize: usize,
) {
    let orig = read_f64le(data);
    let expected = read_f64le(expected);

    let epsilon = 1e-5;

    let coefficients2d_with_blocksize =
        |bspline, image_data: &mut [f64], width, height, epsilon, ext| match (
            row_blocksize,
            col_blocksize,
        ) {
            (1, 1) => coefficients2d(bspline, image_data, width, height, epsilon, ext),
            (6, 7) => coefficients2d_blocked::<_, _, 6, 7>(
                bspline, image_data, width, height, epsilon, ext,
            ),
            (7, 6) => coefficients2d_blocked::<_, _, 7, 6>(
                bspline, image_data, width, height, epsilon, ext,
            ),
            (8, 8) => coefficients2d_blocked::<_, _, 8, 8>(
                bspline, image_data, width, height, epsilon, ext,
            ),
            (r, c) => {
                panic!("unsupported blocksize: R={r}, C={c}")
            }
        };

    let mut image_data = orig.clone();
    coefficients2d_with_blocksize(
        BSpline10::<f64>::default(),
        &mut image_data,
        width,
        height,
        epsilon,
        ext,
    )
    .unwrap();

    assert_approx_eq!(
        &[f64],
        &expected.as_slice(),
        &image_data.as_slice(),
        epsilon = 2. * f64::EPSILON
    );
}

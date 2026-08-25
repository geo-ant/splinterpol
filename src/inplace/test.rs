use std::num::NonZero;

use float_cmp::assert_approx_eq;

use crate::{
    TransmittableBoundaryExtension, inplace::apply_exponential_filter_in_place,
    test_utils::read_f64le,
};

/// reference results generated with my fork of the Briand&Monasse repository,
/// see: https://github.com/geo-ant/reference_bspline_interpolation

#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64",
    TransmittableBoundaryExtension::Mirrored
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.wsym.alpha-0.283.f64",
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case::test_case(
    "test-assets/expfilter/data1.f64",
    "test-assets/expfilter/data1.filtered.n100.peri.alpha-0.283.f64",
    TransmittableBoundaryExtension::Periodic
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.hsym.alpha-0.283.f64",
    TransmittableBoundaryExtension::Mirrored
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.wsym.alpha-0.283.f64",
    TransmittableBoundaryExtension::WholeSymmetric
)]
#[test_case::test_case(
    "test-assets/expfilter/data2.f64",
    "test-assets/expfilter/data2.filtered.n100.peri.alpha-0.283.f64",
    TransmittableBoundaryExtension::Periodic
)]
fn exponential_filter_reproduces_reference_results(
    data: &str,
    expected: &str,
    extension: TransmittableBoundaryExtension,
) {
    let mut data = read_f64le(data);
    let expected = read_f64le(expected);

    let data_len = data.len().try_into().unwrap();
    apply_exponential_filter_in_place(
        -0.283,
        data.as_mut_slice(),
        0,
        NonZero::new(1).unwrap(),
        data_len,
        100.try_into().unwrap(),
        extension,
    )
    .unwrap();

    assert_approx_eq!(&[f64], expected.as_slice(), data.as_slice(), ulps = 4);
}

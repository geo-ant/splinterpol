use std::num::NonZero;

use float_cmp::{Ulps, approx_eq, assert_approx_eq};

use crate::{inplace::apply_exponential_filter_in_place, test_utils::read_f64le};

/// reference results generated with my fork of the Briand&Monasse repository,
/// see: https://github.com/geo-ant/reference_bspline_interpolation

#[test]
fn exponential_filter_reproduces_reference_results() {
    let mut data = read_f64le("test-assets/expfilter/data1.f64");
    let expected = read_f64le("test-assets/expfilter/data1.filtered.n100.hsym.alpha-0.283.f64");

    let data_len = data.len().try_into().unwrap();
    apply_exponential_filter_in_place(
        -0.283,
        data.as_mut_slice(),
        0,
        NonZero::new(1).unwrap(),
        data_len,
        100.try_into().unwrap(),
        crate::TransmittableBoundaryExtension::Mirrored,
    )
    .unwrap();

    assert_approx_eq!(&[f64], expected.as_slice(), data.as_slice(), ulps = 4);
}

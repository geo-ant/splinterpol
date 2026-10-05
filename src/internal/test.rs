use crate::bspline::{BSplinePrefilter, BSpline10, BSpline11};
use crate::{bspline::BSpline3, internal::compute_truncation_indices};

#[test]
fn bspline3_trunctation_indices() {
    // I calculated the trunctation indices ground truth values
    // with the https://github.com/geo-ant/reference_bspline_interpolation
    // codebase
    let eps = 0.001;
    let expected = [8];
    let calc = compute_truncation_indices(&BSpline3::<f64>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.001;
    let calc = compute_truncation_indices(&BSpline3::<f32>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.00001;
    let expected = [11];
    let calc = compute_truncation_indices(&BSpline3::<f64>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.00001;
    let calc = compute_truncation_indices(&BSpline3::<f32>::POLES, eps).unwrap();
    assert_eq!(calc, expected);
}

#[test]
fn bspline11_trunctation_indices() {
    // I calculated the trunctation indices ground truth values
    // with the https://github.com/geo-ant/reference_bspline_interpolation
    // codebase
    let eps = 0.001;
    let expected = [42, 15, 8, 5, 3];
    let calc = compute_truncation_indices(&BSpline11::<f64>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.001;
    let calc = compute_truncation_indices(&BSpline11::<f32>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.00001;
    let expected = [53, 18, 10, 6, 4];
    let calc = compute_truncation_indices(&BSpline11::<f64>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.00001;
    let calc = compute_truncation_indices(&BSpline11::<f32>::POLES, eps).unwrap();
    assert_eq!(calc, expected);
}

#[test]
fn bspline10_trunctation_indices() {
    // I calculated the trunctation indices ground truth values
    // with the https://github.com/geo-ant/reference_bspline_interpolation
    // codebase
    let eps = 0.001;
    let expected = [36, 13, 7, 4, 2];
    let calc = compute_truncation_indices(&BSpline10::<f64>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.001;
    let calc = compute_truncation_indices(&BSpline10::<f32>::POLES, eps).unwrap();
    assert_eq!(calc, expected);

    let eps = 0.00001;
    let expected = [46, 16, 9, 5, 3];
    let calc = compute_truncation_indices(&BSpline10::<f64>::POLES, eps).unwrap();
    assert_eq!(calc, expected);
}

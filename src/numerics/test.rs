use crate::numerics::{FloatVector, FusedMuladd};

#[test]
//NOTE: for this test to meaningfully test something, we have to enable
// the `fma` or `fma,std` feature flags
fn f32_fma() {
    let a = 2f32;
    let b = 7f32;
    let c = -13f32;

    assert_eq!(a.fmuladd(b, c), a * b + c);
}

#[test]
fn f64_fma() {
    let a = 2f32;
    let b = 7f32;
    let c = -13f32;

    assert_eq!(a.fmuladd(b, c), a * b + c);
}

#[test]
fn float_vector_for_array() {
    let first = [12f64, 20., 30.];
    let second = [-10., -4., 31.];

    assert_eq!(first.axpy(2., second), [14., 36., 91.]);
    assert_eq!(second.axpy(-3., first), [42., 32., -63.]);

    assert_eq!(first.scal(4.), [48., 80., 120.]);
    assert_eq!(second.scal(-2.), [20., 8., -62.]);

    assert_eq!(first.sub(second), [22., 24., -1.]);
    assert_eq!(second.sub(first), [-22., -24., 1.]);

    assert_eq!(first.dot(second), second.dot(first));
    assert_eq!(first.dot(second), 730.);
}

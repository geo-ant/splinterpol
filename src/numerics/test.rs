use crate::numerics::FusedMuladd;

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

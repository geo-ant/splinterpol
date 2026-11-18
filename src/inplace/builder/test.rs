use crate::{
    bspline::BSpline3,
    inplace::{InplaceBoundaryExtension, SplineInterpolator1d, SplineInterpolator2d},
};

#[test]
fn spline_interpolator_uses_whole_symmetric_extension_by_default() {
    assert_eq!(
        SplineInterpolator2d::builder(BSpline3::<f64>::new()).boundary_extension,
        InplaceBoundaryExtension::WholeSymmetric
    );
    assert_eq!(
        SplineInterpolator1d::builder(BSpline3::<f64>::new()).boundary_extension,
        InplaceBoundaryExtension::WholeSymmetric
    );
}

#[test]
fn spline_interpolator_setters() {
    // NOTE: this is the only setter so far that affects a runtime value.
    assert_eq!(
        SplineInterpolator2d::builder(BSpline3::<f64>::new())
            .extension(InplaceBoundaryExtension::Periodic)
            .boundary_extension,
        InplaceBoundaryExtension::Periodic
    );
    assert_eq!(
        SplineInterpolator1d::builder(BSpline3::<f64>::new())
            .extension(InplaceBoundaryExtension::Periodic)
            .boundary_extension,
        InplaceBoundaryExtension::Periodic
    );
}

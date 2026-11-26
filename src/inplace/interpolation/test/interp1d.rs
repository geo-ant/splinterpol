// this is a long ass test...
use crate::inplace::SplineInterpolator1d;
use crate::test_utils::Const;
use crate::{bspline::*, inplace::InplaceBoundaryExtension};
use float_cmp::assert_approx_eq;
use test_case::test_case;
// BSpline3
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline3::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline3::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline4
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline4::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline4::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline5
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline5::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline5::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline6
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline6::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline6::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline7
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline7::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline7::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline8
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline8::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline8::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline9
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline9::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline9::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline10
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline10::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline10::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline11
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline11::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline11::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
// BSpline12
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-6, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Periodic, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::WholeSymmetric, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<1>)]
#[test_case(BSpline12::<f64>::new(),1021, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),1024, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2069, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
#[test_case(BSpline12::<f64>::new(),2048, InplaceBoundaryExtension::Reflected, 1e-12, Const::<8>)]
fn input_is_reproduced_on_original_grid_1d<
    B,
    const NP: usize,
    const NW: usize,
    const ROW_BLOCKSIZE: usize,
>(
    bspline: B,
    len: usize,
    ext: InplaceBoundaryExtension,
    epsilon: f64,
    // just so it can be given as an argument
    _rb: Const<ROW_BLOCKSIZE>,
) where
    B: BSpline<Float = f64> + BSplinePrefilter<NP> + BSplineWeights<NW>,
{
    // we use a synthethic image here for 1d signals
    // which is a sine with chirp that's multiplied with an exponential
    let fstart = 0.001;
    let fend = 0.45;
    let flen = len as f64;
    let f = |x| {
        let x = x as f64;
        (100. * f64::exp(-x / (flen / 3.)))
            * f64::sin(
                2. * core::f64::consts::PI
                    * (fstart * x + (fend - fstart) / (2. * flen - 2.) * x.powi(2)),
            )
    };
    let data: Vec<_> = (0..len).map(f).collect();

    let fmax = data
        .iter()
        .copied()
        .map(f64::abs)
        .max_by(|a, b| f64::total_cmp(a, b))
        .unwrap();

    let interp = SplineInterpolator1d::builder(bspline)
        .extension(ext)
        .row_blocksize::<ROW_BLOCKSIZE>()
        .build(data.clone(), epsilon)
        .unwrap();

    for pos in 0..data.len() {
        assert_approx_eq!(
            f64,
            data[pos],
            interp.interp_single(pos as _).unwrap(),
            // NOTE: the B&M paper specifies that epsilon is a relative
            // bound so this should be the right bound to enforce here.
            epsilon = epsilon * fmax
        );
    }
}

#[cfg(feature = "__compare_mkl")]
#[test_case(128,  1e-6, Const::<1>)]
#[test_case(128,  1e-12, Const::<1>)]
#[test_case(128,  1e-6, Const::<8>)]
#[test_case(128,  1e-12, Const::<8>)]
#[test_case(129,  1e-6, Const::<1>)]
#[test_case(129,  1e-12, Const::<1>)]
#[test_case(129,  1e-6, Const::<8>)]
#[test_case(129,  1e-12, Const::<8>)]
#[test_case(1024, 1e-6, Const::<1>)]
#[test_case(1024, 1e-12, Const::<1>)]
#[test_case(1024, 1e-6, Const::<8>)]
#[test_case(1024, 1e-12, Const::<8>)]
#[test_case(1025, 1e-6, Const::<1>)]
#[test_case(1025, 1e-12, Const::<1>)]
#[test_case(1025, 1e-6, Const::<8>)]
#[test_case(1025, 1e-12, Const::<8>)]
fn bspline3_interpolation_wholesym_vs_mkl<const ROW_BLOCKSIZE: usize>(
    len: usize,
    epsilon: f64,
    _rb: Const<ROW_BLOCKSIZE>,
) {
    use mkl_spline::*;

    let bspline = BSpline3::new();
    let ext = InplaceBoundaryExtension::WholeSymmetric;
    // we use a synthethic image here for 1d signals
    // which is a sine with chirp that's multiplied with an exponential
    let fstart = 0.001;
    let fend = 0.45;
    let flen = len as f64;
    let f = |x| {
        let x = x as f64;
        (100. * f64::exp(-x / (flen / 3.)))
            * f64::sin(
                2. * core::f64::consts::PI
                    * (fstart * x + (fend - fstart) / (2. * flen - 2.) * x.powi(2)),
            )
    };
    let data: Vec<_> = (0..len).map(f).collect();

    let fmax = data
        .iter()
        .copied()
        .map(f64::abs)
        .max_by(|a, b| f64::total_cmp(a, b))
        .unwrap();

    let splinterp = SplineInterpolator1d::builder(bspline)
        .extension(ext)
        .row_blocksize::<ROW_BLOCKSIZE>()
        .build(data.clone(), epsilon)
        .unwrap();

    let mut mkl_coeffs = vec![0.; 4 * (len - 1)];
    let mut mklinterp = MklInterpF64::new(&data, &mut mkl_coeffs).unwrap();
    mklinterp.calculate_coefficients().unwrap();

    // interpolating with step of 0.1
    let positions: Vec<f64> = (0..=(10 * (len - 1))).map(|j| (j as f64) / 10.).collect();
    let mut splinterpol_results = vec![0.; positions.len()];
    let mut mkl_results = vec![-1.; positions.len()];

    splinterp
        .interp_multi(&positions, &mut splinterpol_results)
        .unwrap();
    mklinterp
        .interp_multi_generic(&positions, &mut mkl_results)
        .unwrap();

    // this is still a relative comparison, but the cool part is that this is
    // true even for interpolation
    assert_approx_eq!(
        &[f64],
        splinterpol_results.as_slice(),
        mkl_results.as_slice(),
        epsilon = epsilon * fmax
    );
}

#[test]
#[cfg(target_arch = "x86_64")]
#[cfg(not(feature = "__compare_mkl"))]
fn enable_mkl_on_my_machine() {
    if std::env::var_os("CI").is_none() {
        // NOTE(geo-ant) I know this is a horrible horrible hack, but it will remind me
        // to enable that feature on my PC. If you're not me, feel free to disable
        // the test...
        panic!(
            "Running on my local x64 machine? Enable the `__compare_mkl` feature (sorry if you're not me)"
        );
    }
}

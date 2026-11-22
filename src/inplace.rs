//! Contains low-level and high-level logic for spline-interpolation of
//! 1D and 2D images using an _in-place_ calculation strategy. It's
//! recommended to use the high-level interface unless you can't for some reason.
//! Both ways of interacting with this library are completely **allocation-free**.
//!
//! ## High-level Interface
//!
//! The [`SplineInterpolator1d`] and [`SplineInterpolator2d`] are a
//! zero-overhead, high level interface to the 1D and 2D interpolation functionality
//! of this library. It's the recommended way of interacting with this library
//! and offers the same configuration options as the low-level interface.
//! The low-level interface is exposed for users who want direct exposure to
//! the 2-step process that this library uses for interpolation.
//!
//! ## Low-level Interface
//!
//! Spline interpolation is a two-step process, where the first step is to calculate
//! the spline coefficients on the whole image. Then, the interpolation
//! can be performed on arbitrary row and column indices inside the image.
//! In this module, the calculation of the spline coefficients happens _in-place_,
//! which means no extra allocations are necessary, but the original image
//! data will be overwritten with the coefficients.
//!
//! ### Step 1: Coefficient Calculation
//!
//! To calculate the spline coefficients for a 2D image, the most low-level
//! way is to use the [`coefficients2d`] function directly, whose
//! docs explain the parameters and give some recommendations.
//!
//! For 1D data, the process is analogous, we just use [`coefficients1d`].
//!
//! ```rust
//! # fn foo() {
//! # use splinterpol::*;
//! # use splinterpol::bspline::*;
//! # let mut some_data = [10.;90];
//! # let epsilon = 1e-1;
//! # let width = core::num::NonZero::new(9).unwrap();
//! # let height = core::num::NonZero::new(10).unwrap();
//! # use splinterpol::inplace::InplaceBoundaryExtension;
//! let mut image_inout = RowMajor {
//!     data: some_data.as_mut(),
//!     width,
//!     height,
//! };
//! inplace::coefficients2d(
//!     &BSpline3::<f64>::new(),
//!     image_inout,
//!     epsilon,
//!     // default is whole-symmetric, same as scipy
//!     InplaceBoundaryExtension::default(),
//!     // choose this according to your machine, see function docs for guidance
//!     RowBlockSize::<8>,
//!     // default row and column block-sizes are 1
//!     ColBlockSize::default(),
//! )
//! .unwrap();
//!
//! # }
//! ```
//!
//! At this point, the image data in `some_data` is overwritten with the
//! spline coefficients and can be used for interpolation.
//!
//! ### Step 2: Interpolation
//!
//! Now that we have the coefficients, we can use them for interpolation using
//! the second half of the low-level interface. The general purpose entry point
//! for interpolating an arbitraty floating point coordinate on the grid are
//! [`interp1d_single`] and [`interp2d_single`]. The
//! following code shows how to use the previously calculated coefficients
//! for interpolation.
//!
//! **Attention**: the same B-Spline type `B` and boundary
//! extension that was used to create the coefficients must be used for
//! interpolation or the results will be incorrect. Using the aforementioned
//! high-level interface takes care of this automatically.
//!
//! ```rust
//! # use splinterpol::*;
//! # use splinterpol::bspline::*;
//! # use core::num::NonZeroUsize;
//! # use splinterpol::inplace::interp2d_single;
//! # let epsilon = 1e-1;
//! # let width = core::num::NonZero::new(9).unwrap();
//! # let height = core::num::NonZero::new(10).unwrap();
//! # use splinterpol::inplace::InplaceBoundaryExtension;
//! # fn foo(some_data:& [f64], width: NonZeroUsize, height: NonZeroUsize, boundary_ext: InplaceBoundaryExtension) {
//! let mut coeffs = RowMajor {
//!   data: some_data.as_ref(),
//!   width,
//!   height
//! };
//!
//! //                          must be the same type as was used
//! //                          for calculating the coefficients
//! //                           |
//! //                           |              always inferred, leave as _
//! //                           |               |
//! //                           v               v
//! let val = interp2d_single::<BSpline3::<f64>, _>(
//!   coeffs,
//!   // boundary extension must be the same that was
//!   // used for computing coefficients
//!   boundary_ext,
//!   Coord::row(3.7).col(1.3)
//! ).unwrap();
//!
//! # }
//! ```
//!
//! ## References
//!
//! * [B&M](https://www.ipol.im/pub/art/2018/221/) T. Briand, P. Monasse: "Theory
//!   and Practice of Image B-Spline Interpolation"

#[cfg(test)]
mod test;

/// some helper code that contains builder logic to call the function because
/// generic parameters and const generics can make the function slightly awkward
/// to call.
mod builder;

/// contains the logic for coefficient calculation in-place
pub(crate) mod coefficients;

/// interpolation logic... duh
pub(crate) mod interpolation;

/// boundary extension enum and logic
mod boundary;

pub use boundary::InplaceBoundaryExtension;
pub use builder::SplineInterpolator1d;
pub use builder::SplineInterpolator2d;
pub use coefficients::coefficients1d;
pub use coefficients::coefficients2d;
pub use interpolation::interp1d_single;
pub use interpolation::interp2d_single;

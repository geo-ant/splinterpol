//! Contains low-level and high-level logic for spline-interpolation of
//! 1D and 2D images using an _in-place_ calculation strategy. It's
//! recommended to use the high-level interface unless you can't for some reason.
//! Both ways of interacting with this library are completely **allocation-free**.
//!
//! ## High-level Interface
//!
//! The [`SplineInterpolator1d`](crate::inplace::SplineInterpolator1d) and
//! [`SplineInterpolator2d`](crate::inplace::SplineInterpolator2d) are a
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
//! way is to use the [`coefficients2d_blocked`] function directly, whose
//! docs explain the parameters and give some recommendations. However, that
//! function is quite clunky to use, so prefer to use the builder style
//! [`Coefficients2d`] wrapper for better ergonomics. Again, the parameters
//! are documented in [`coefficients2d_blocked`], so please consult its docs.
//!
//! For 1D data, the process is analogous, we just use [`coefficients1d_blocked`]
//! or [`Coefficients1d`] for calculating the coefficients.
//!
//! ```rust
//! # fn foo() {
//! # use splinterpol::*;
//! # use splinterpol::bspline::*;
//! # let mut my_image_data = [10.;90];
//! # let epsilon = 1e-1;
//! # let width = core::num::NonZero::new(9).unwrap();
//! # let height = core::num::NonZero::new(10).unwrap();
//! # let extension = splinterpol::inplace::InplaceBoundaryExtension::Periodic;
//!
//! inplace::Coefficients2d::new(BSpline3::<f64>::new())
//!     // optional, uses block-size 1 if not given
//!     .with_row_blocksize::<4>()
//!     // optional, uses block-size 1 if not given
//!     .with_column_blocksize::<8>()
//!     .compute(
//!         RowMajor {
//!           data: &mut my_image_data,
//!           width,
//!           height,
//!         },
//!         epsilon,
//!         extension,
//!     ).unwrap();
//! # }
//! ```
//!
//! At this points, `my_image_data` is overwritten with the spline coefficients
//! and can be used for interpolation.
//!
//! ### Step 2: Interpolation
//!
//! Now that we have the coefficients, we can use them for interpolation in this
//! low-level interface. The general purpose entry point for interpolating
//! a point on the grid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!
//! TODO TODO
//! TODO:documentation
//!
//!
//!
//! ## References
//!
//! * [B&M](https://www.ipol.im/pub/art/2018/221/) T. Briand, P. Monasse: "Theory
//! and Practice of Image B-Spline Interpolation"

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
pub use builder::Coefficients2d;
pub use builder::SplineInterpolator1d;
pub use builder::SplineInterpolator2d;
pub use coefficients::coefficients1d_blocked;
pub use coefficients::coefficients2d_blocked;
pub use interpolation::interp2d_single;

//! Contains low-level and high-level logic for spline-interpolation of
//! 1D and 2D images using an _in-place_ calculation strategy. It's
//! recommended to use the high-level interface unless you can't for some reason.
//! Both ways of interacting with this library are completely **allocation-free**.
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
//! ```rust
//! # fn foo() {
//! # use spleen::*;
//! # use spleen::bspline::*;
//! # let mut my_image_data = [10.;90];
//! # let epsilon = 1e-1;
//! # let width = core::num::NonZero::new(9).unwrap();
//! # let height = core::num::NonZero::new(10).unwrap();
//! # let extension = spleen::inplace::InplaceBoundaryExtension::Periodic;
//!
//! inplace::Coefficients2d::new(BSpline3::<f64>::new())
//!     .with_row_blocksize::<2>()
//!     .with_column_blocksize::<3>()
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
//! TODO!!
//!
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
pub use builder::SplineInterpolator2d;
pub use coefficients::coefficients2d_blocked;
pub use interpolation::interp2d_single;

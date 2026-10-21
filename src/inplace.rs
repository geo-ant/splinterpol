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

pub mod boundary;

pub use boundary::TransmittableBoundaryExtension;
pub use builder::Coefficients2d;
pub use builder::RowMajor;
pub use builder::SplineInterpolator2d;
pub use coefficients::coefficients2d;
pub use coefficients::coefficients2d_blocked;
pub use interpolation::interp2d_single;

#![doc=include_str!("../README.md")]
#![warn(missing_docs)]
// lesson learned from interpn: needless range loops can sometimes be pretty
// neat for performance
#![allow(clippy::needless_range_loop)]
#![cfg_attr(not(any(test, feature = "std")), no_std)]
#![deny(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::panic_in_result_fn,
    clippy::panicking_unwrap,
    unsafe_code,
    clippy::alloc_instead_of_core
)]

// the `no_std` feature isn't actually used in the code, but it has significance
// for the dependency selection in the toml file. On `no_std` we use the `libm`
// feature of the `num_traits` crate whereas on `std` we use the `std` feature
// of that crate. Since `no_std` is a default feature of _this_ crate, we
// can't simply pass `features = ["std"]` when adding this crate as a dependency
// in a downstream `Cargo.toml` file. That's why we have to pass
// `default-features=false` as well. Hope that clears things up.
#[cfg(all(feature = "std", feature = "no_std"))]
compile_error!(
    "The `std` and `no_std` feature must not be used together. Use `default-features=false` and enable the desired feature explicitly."
);

use core::num::NonZero;

// test the absolute accuracy of this crate's interpolation on a couple
// of test images (from ImageJ/FIJI) against interpn.
#[cfg(test)]
#[cfg(feature = "__interpn")]
mod interpn_accuracy;

// #[cfg(feature = "__bench_internal")]
#[doc(hidden)]
pub mod bench_exports;

pub mod inplace;
pub(crate) mod shared;
pub use shared::ColBlockSize;
pub use shared::Coord;
pub use shared::RowBlockSize;
pub use shared::RowMajor;
pub use shared::UnfinishedCoord;

mod internal;
mod iter_utils;
mod util;

mod sealed;

/// B-splines for interpolation
pub mod bspline;

#[cfg(test)]
#[macro_use]
mod test_utils;

mod numerics;
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
/// This crate's error type.
pub enum Error {
    /// not enough data to perform the interpolation with the required precision
    NotEnoughData,
    /// no data at all
    NoData,
    /// some given indices or lengths imply out of bounds access to the signal
    OutOfBounds,
    /// data size limits exceeded
    Limits,
    /// internal logic error
    Internal,
    /// data doesn't match the given dimensions
    Dimensions,
    /// overflow occurred
    Overflow,
    /// precision set to a unreasonably small (or forbidden) value
    Precision,
    /// given number is unrepresentable as a float or integer in the
    /// current number type
    Unrepresentable,
}

/// re-export some dependencies with significance for the public interface
/// to this crate
pub mod reexports {
    /// the `num-traits` crate
    pub use num_traits;
}

#[inline(never)]
#[allow(clippy::unwrap_used, missing_docs)]
pub fn blocked_expf64_strided_iter(
    // the pole, -1 < alpha < 0
    alpha: f64,
    signal: &mut [f64],
    stride: NonZero<usize>,
    // NOTE: in theory, this N can be larger than the amount of elements
    // but we're restricting it here
    // TODO: if we want to allow this to be larger, we would have to have
    // temporary storage for the original elements and implement some wraparound,
    // which makes things annoying and complicated I think and defeats the purpose
    // of keeping all this nice and in place. As-is, we can just use one variable
    // for temporary storage.
    // NOTE: in the Briand and Monasse paper N is the truncation INDEX,
    // not the number of elements to be truncated. The element at index N is
    // counted for the initialization.
    truncation_idx: usize,
    boundary_extension: inplace::InplaceBoundaryExtension,
) -> Result<(), Error> {
    let col = iter_utils::BlockedStridedSliceMut::<8, f64>::new(
        signal,
        crate::iter_utils::StridedAccessArgs {
            offset: 0,
            stride,
            count: NonZero::new(signal.len() / stride.get()).unwrap(),
        },
    )
    .ok_or(Error::Dimensions)?;
    // inplace::apply_exponential_filter_in_place_1d(alpha, signal, truncation_idx, boundary_extension)
    inplace::coefficients::apply_exponential_filter_in_place_1d_iter_vector(
        alpha,
        col,
        truncation_idx,
        boundary_extension,
    )
}

#[inline(never)]
#[allow(clippy::unwrap_used, missing_docs)]
pub fn expf64_strided_iter(
    // the pole, -1 < alpha < 0
    alpha: f64,
    signal: &mut [f64],
    stride: NonZero<usize>,
    //NOTE: in theory, this N can be larger than the amount of elements
    // but we're restricting it here
    // TODO: if we want to allow this to be larger, we would have to have
    // temporary storage for the original elements and implement some wraparound,
    // which makes things annoying and complicated I think and defeats the purpose
    // of keeping all this nice and in place. As-is, we can just use one variable
    // for temporary storage.
    // NOTE: in the Briand and Monasse paper N is the truncation INDEX,
    // not the number of elements to be truncated. The element at index N is
    // counted for the initialization.
    truncation_idx: usize,
    boundary_extension: inplace::InplaceBoundaryExtension,
) -> Result<(), Error> {
    let col = iter_utils::StridedSliceMut::new(
        signal,
        crate::iter_utils::StridedAccessArgs {
            offset: 0,
            stride,
            count: NonZero::new(signal.len() / stride.get()).unwrap(),
        },
    )
    .ok_or(Error::Dimensions)?;
    // inplace::apply_exponential_filter_in_place_1d(alpha, signal, truncation_idx, boundary_extension)
    inplace::coefficients::apply_exponential_filter_in_place_1d_iter(
        alpha,
        col,
        truncation_idx,
        boundary_extension,
    )
}

// #[inline(never)]
// pub fn expf64_strided(
//     // the pole, -1 < alpha < 0
//     alpha: f64,
//     signal: &mut [f64],
//     stride: core::num::NonZero<usize>,
//     count: core::num::NonZero<usize>,
//     //NOTE: in theory, this N can be larger than the amount of elements
//     // but we're restricting it here
//     // TODO: if we want to allow this to be larger, we would have to have
//     // temporary storage for the original elements and implement some wraparound,
//     // which makes things annoying and complicated I think and defeats the purpose
//     // of keeping all this nice and in place. As-is, we can just use one variable
//     // for temporary storage.
//     // NOTE: in the Briand and Monasse paper N is the truncation INDEX,
//     // not the number of elements to be truncated. The element at index N is
//     // counted for the initialization.
//     truncation_idx: usize,
//     boundary_extension: TransmittableBoundaryExtension,
// ) -> Result<(), Error> {
//     unsafe {
//         inplace::apply_exponential_filter_in_place_strided(
//             alpha,
//             signal.as_mut_ptr(),
//             stride,
//             count,
//             truncation_idx,
//             boundary_extension,
//         )
//     }
// }

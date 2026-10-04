// TODO make this into a `deny` lint conditional on not having some feature
// `with_unsafe` or `allow_unsafe` in there...
#![warn(unsafe_code)]
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
    // unsafe_code,
    clippy::alloc_instead_of_core
)]

#[cfg(all(feature = "std", feature = "no_std"))]
compile_error!(
    "The `std` and `no_std` feature may not be enabled together. Use `default-features=false` and enable the desired feature explicitly."
);

use core::num::NonZero;

#[cfg(feature = "__bench_internal")]
#[doc(hidden)]
pub mod bench_exports;

pub mod inplace;

mod internal;
mod iter_utils;

/// B-splines for interpolation
pub mod bspline;

#[cfg(test)]
#[macro_use]
mod test_utils;

mod numerics;
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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
}

/// Describes a way to extend the signal beyond the measured domain in a way
/// that can be _transmitted_ like described in section 3.2.2 of the
/// Briand & Monasse paper (cf also Algorithm 5 therein). Note that constant
/// extension can _not_ be written as a transmittable boundary condition.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum TransmittableBoundaryExtension {
    /// Signal `abcde` is extended to `...cba|abcde|edc` and is called
    /// _half-symmetric_ in the Briand & Monasse paper.
    Mirrored = 1,
    /// Signal `abcde` is extended to `...dcb|abcde|dcb...`
    WholeSymmetric = 2,
    /// Signal `abcde` is extended to `...cde|abcde|abc...`
    Periodic = 3,
}

/// re-export some dependencies with significance for the public interface
/// to this crate
pub mod reexports {
    /// the `num-traits` crate
    pub use num_traits;
}

#[inline(never)]
#[allow(clippy::unwrap_used)]
pub fn blocked_expf64_strided_iter(
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
    boundary_extension: TransmittableBoundaryExtension,
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
    inplace::apply_exponential_filter_in_place_1d_iter_vector(
        alpha,
        col,
        truncation_idx,
        boundary_extension,
    )
}

#[inline(never)]
#[allow(clippy::unwrap_used)]
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
    boundary_extension: TransmittableBoundaryExtension,
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
    inplace::apply_exponential_filter_in_place_1d_iter(
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

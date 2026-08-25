#![cfg_attr(not(any(test, feature = "std")), no_std)]

use num_traits::Float;

mod inplace;
mod strided;

#[cfg(test)]
mod test_utils;

/// trait exposing a fused mul-add or a fallback if the instruction doesn't exist
/// Idea taken from the interpN crate.
pub(crate) trait FusedMuladd {
    fn fmuladd(self, factor: Self, addend: Self) -> Self;
}

impl<T> FusedMuladd for T
where
    T: Float,
{
    #[inline(always)]
    fn fmuladd(self, factor: Self, addend: Self) -> Self {
        #[cfg(feature = "fma")]
        {
            self.mul_add(factor, addend)
        }
        #[cfg(not(feature = "fma"))]
        {
            self * factor + addend
        }
    }
}

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
}

/// Describes a way to extend the signal beyond the measured domain in a way
/// that can be _transmitted_ like described in section 3.2.2 of the
/// Briand & Monasse paper (cf also Algorithm 5 therein). Note that constant
/// extension can _not_ be written as a transmittable boundary condition.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum TransmittableBoundaryExtension {
    /// Signal `abcde` is extended to `...cba|abcde|edc` and is called
    /// _half-symmetric_ in the Briand & Monasse paper.
    Mirrored,
    /// Signal `abcde` is extended to `...dcb|abcde|dcb...`
    WholeSymmetric,
    /// Signal `abcde` is extended to `...cde|abcde|abc...`
    Periodic,
}

/// re-export some dependencies with significance for the public interface
/// to this crate
pub mod reexports {
    /// the `num-traits` crate
    pub use num_traits;
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}

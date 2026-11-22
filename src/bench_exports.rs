use crate::{Error, inplace::InplaceBoundaryExtension, numerics::FusedMuladd};
use core::iter::Sum;
use num_traits::{ConstOne, ConstZero, Float};

pub use crate::internal::compute_truncation_indices;

pub fn apply_exponential_filter_in_place<F, const BLOCKSIZE: usize>(
    // the pole, -1 < alpha < 0
    alpha: F,
    signal: &mut [F],
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
    boundary_extension: InplaceBoundaryExtension,
) -> Result<(), Error>
where
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
{
    crate::inplace::coefficients::apply_exponential_filter_in_place_1d::<_, BLOCKSIZE>(
        alpha,
        signal,
        truncation_idx,
        boundary_extension,
    )
}

pub fn apply_exponential_filter_in_place_iter<F>(
    // the pole, -1 < alpha < 0
    alpha: F,
    signal: &mut [F],
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
    boundary_extension: InplaceBoundaryExtension,
) -> Result<(), Error>
where
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
{
    crate::inplace::coefficients::apply_exponential_filter_in_place_1d_iter(
        alpha,
        signal,
        truncation_idx,
        boundary_extension,
    )
}

use core::iter::Sum;
use num_traits::{ConstOne, ConstZero, Float};

use crate::{
    Error, FusedMuladd, TransmittableBoundaryExtension,
    inplace::apply_exponential_filter_in_place_1d as apply_exponential_filter_in_place_internal,
};

pub fn apply_exponential_filter_in_place<F>(
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
    boundary_extension: TransmittableBoundaryExtension,
) -> Result<(), Error>
where
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
{
    // crate::inplace::apply_exponential_filter_in_place_strided(
    //     alpha,
    //     signal.as_mut_ptr(),
    //     const { NonZero::new(1).unwrap() },
    //     NonZero::new(signal.len()).unwrap(),
    //     truncation_idx,
    //     boundary_extension,
    // )
    apply_exponential_filter_in_place_internal(alpha, signal, truncation_idx, boundary_extension)
}

#[deny(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::panic_in_result_fn,
    clippy::panicking_unwrap
)]
use core::iter::Sum;
use core::num::NonZero;
use num_traits::{ConstOne, ConstZero, Float};

use super::Error;
use crate::{
    FusedMuladd, TransmittableBoundaryExtension,
    strided::{StridedSliceArgs, StridedSliceMut},
};

#[cfg(test)]
mod test;

// TODO: figure out the data layout here.. should I use slices and steps, pointers
// and steps, or iterators? Iterators are probably the most idiomatic, but
// I'm not sure I can make it work and I'm not sure if this optimizes well...
// TODO #2: I think for the 2D case, there might be a better way of doing that
// by doing the "y-axis" alongside the x-axis, because every load of the
// non-contiguous axis loads a number of elements from the contiguous axis
// into the cache which should allow us to make progress on the contiguous axis
// as well...
// NOTE: apply exponential filter in the original domain
/// This is algorithm 3 from the Briand & Monasse paper for approach 2 with a
/// transmitted boundary condition, meaning the data is processed only on the
/// original domain.
#[expect(dead_code)]
pub(crate) fn apply_exponential_filter_in_place<F>(
    // the pole, -1 < alpha < 0
    alpha: F,
    signal: &mut [F],
    first: usize,
    stride: NonZero<usize>,
    count: NonZero<usize>,
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
    F: Float + ConstOne + ConstZero + Sum,
{
    if signal.is_empty() {
        return Err(Error::NoData);
    }

    let mut signal = StridedSliceMut::new(
        StridedSliceArgs {
            first,
            stride,
            count,
        },
        signal,
    )
    .ok_or(Error::Internal)?;

    if truncation_idx >= count.get() {
        return Err(Error::NotEnoughData);
    }

    // this should not happen in practice but this is important to ensure that
    // the i32 cast of the index based on the truncation number doesn't
    // exceed the i32 limit
    if truncation_idx > i32::MAX as usize {
        return Err(Error::Limits);
    }

    debug_assert!(
        -F::ONE < alpha && alpha < F::ZERO,
        "internal error: pole must be in (-1,0)"
    );

    // eq (48) from the B&M paper
    // NOTE: this also takes care of line 9 in algorithm 5 because this
    // calculates the initialization value for the coefficients as well using the
    // boundary conditions. It doesn't matter whether the given signal is
    // the original signal or one of the filtering steps for the coefficients
    // since those obey the same (transmitted) boundary conditions.
    let s_lini = {
        let signal_first = signal.first().copied().ok_or(Error::NotEnoughData)?;
        match boundary_extension {
            TransmittableBoundaryExtension::Mirrored => {
                // NOTE: this is true for all the boundary code. I'm reusing the
                // same logic of recursively calculating a^j = a * a^(j-1) that
                // the reference C-code uses. If I use a^j = a.powi(j), I might
                // get better ILP (and actually better floating point accuracy
                // I think, but the code doesn't reproduce the float-identical
                // results anymore.
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
            TransmittableBoundaryExtension::WholeSymmetric => {
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.iter().copied().skip(1).take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
            TransmittableBoundaryExtension::Periodic => {
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.iter().rev().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
        }
    };

    signal[0] = s_lini;

    // causal filtering, see algo 3 in the B&M paper lines 2-4
    // NOTE: this is probably a good place to optimize since we have pretty
    // bad instruction level parallelism here because of the recursion.
    // TODO: think about that. I think a good paper that does this for
    // CUBIC splines is Champagnat et. al, but they take it many steps further
    // by writing the filtering as a convolution with a 15-element kernel,
    // see: https://www.tandfonline.com/doi/abs/10.1080/2165347X.2013.824736
    // or author link: https://www.onera.fr/sites/default/files/Departements-scientifiques/DTIM/jgt.pdf
    //
    // But maybe this paper should also be interesting: https://dl.acm.org/doi/10.1145/2070781.2024210
    //
    // Even for CPU there should be a way to make this work on MxM tiles
    // (where M is chosen as something like max(truncation_n, next multiple of number of float elements in cacheline)
    // and then we observe that we can apply the filtering for rows and columns
    // of 2D images as causal filter partial row, causal filter partial col,
    // anti-causal filter row, anti-causal filter col.
    //
    // The first 2 filter steps and the last 2 can be performed on tiles,
    // and there are probably better ways of writing the recurrence relation
    // inside those tiles by unrolling it explicitly to take advantage of
    // more SIMD and ILP.
    //
    // NOTE 2: there might even be a way to take advantage of this inside
    // the rows themselves, but then we lose the advantage of cache-efficiency
    // over rows and columns... yes we can definitely do that where we "block"
    // in 1D, then calculate the recurrence relation inside the 1D block and then
    // transmit the boundardy values across to all other blocks meaning boundary
    // value for block i is added to all blocks j with j > i.
    // NOTE 3: this is probably very simd-unfriendly since we have a stride
    // in there as well... not only does this have terrible ILP, the compiler
    // likely can't even use SIMD instructions here

    for idx in 1..signal.len() {
        // eq. (39) in the Briand and Monasse paper, cf also algo 3, line 3.
        signal[idx] = signal[idx - 1].fmuladd(alpha, signal[idx]);
    }

    // line 5 in B&M algorithm 3, but we're not using eq (49) directly but the
    // special cases in eqs. (55)-(66)
    let s_lend = {
        let signal_last = signal.last().copied().ok_or(Error::NotEnoughData)?;
        match boundary_extension {
            TransmittableBoundaryExtension::Periodic => {
                let mut acc = signal_last;
                let mut pow = F::ONE;
                for s in signal.iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                -alpha * acc
            }
            TransmittableBoundaryExtension::Mirrored => alpha / (alpha - F::ONE) * signal_last,
            TransmittableBoundaryExtension::WholeSymmetric => {
                let signal_2nd_to_last = signal
                    .get(signal.len() - 2)
                    .copied()
                    .ok_or(Error::NotEnoughData)?;
                alpha / (alpha * alpha - F::ONE) * (signal_last + alpha * signal_2nd_to_last)
            }
        }
    };

    let l_end = signal.len() - 1;
    signal[l_end] = s_lend;

    for idx in (0..l_end).rev() {
        signal[idx] = alpha * (signal[idx + 1] - signal[idx]);
    }

    Ok(())
}

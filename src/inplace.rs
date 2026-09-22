use core::iter::Sum;
use core::num::NonZero;
use num_traits::{ConstOne, ConstZero, Float};

use super::Error;
use crate::{
    TransmittableBoundaryExtension,
    bspline::{BSpline, BSpline3},
    fma::FusedMuladd,
    internal::compute_truncation_indices,
    iter_utils::{StridedPointerAccess, StridedSliceMut},
};

#[cfg(test)]
mod test;

// just to check that the function can be called easily
fn _foo() {
    _ = coefficients2d(
        BSpline3::<f64>::default(),
        &mut [1.],
        100,
        100,
        0.1,
        TransmittableBoundaryExtension::Mirrored,
    );
}

pub fn coefficients2d<B, const NPOLES: usize>(
    _bspline: B,
    image_data: &mut [B::Float],
    width: usize,
    height: usize,
    epsilon: B::Float,
    boundary_extension: TransmittableBoundaryExtension,
) -> Result<(), Error>
where
    B::Float: Float + ConstOne + ConstZero + Sum + FusedMuladd,
    B: BSpline<NPOLES>,
{
    let expected_len = match width.overflowing_mul(height) {
        (val, false) => val,
        (_, true) => return Err(Error::Overflow),
    };

    if image_data.len() != expected_len {
        return Err(Error::Dimensions);
    }

    let width = NonZero::new(width).ok_or(Error::NoData)?;
    let height = NonZero::new(height).ok_or(Error::NoData)?;

    // NOTE: this is almost certainly not the most efficient way of doing things,
    // since the columns have massive stride typically,so a blocking strategy
    // is the way to go here in general

    let truncation_indices = compute_truncation_indices(&B::POLES, epsilon)?;

    // column filtering
    // TODO: remove this horrible unsafe code. We should not have to do any of this
    // for 2D sequential row-by row iteration. We can write an abstraction that
    // safely returns column iterators for a row major mutable view of an array...
    // but for parallel iterators we have to use unsafe code. But maybe I
    // don't have to do any iterator stuff and I can use a blocking strategy
    // with manual indexing right away anyways, so let's see if we can get this
    // working and see how well it performs compared to the reference C impl

    for (pole, truncation_index) in B::POLES
        .iter()
        .copied()
        .zip(truncation_indices.iter().copied())
    {
        for col in 0..width.get() {
            let col = StridedSliceMut::new(
                image_data,
                crate::iter_utils::StridedAccessArgs {
                    offset: col,
                    stride: width,
                    count: height,
                },
            )
            .ok_or(Error::Dimensions)?;
            apply_exponential_filter_in_place_1d_iter(
                pole,
                col,
                truncation_index,
                boundary_extension,
            )?;

            // TODO PERF: with the access checks for at and at_mut disabled
            // this is almost at C reference performance
            // unsafe {
            //     let col = image_data.as_mut_ptr().add(col);
            //     apply_exponential_filter_in_place_strided(
            //         pole,
            //         col,
            //         width,
            //         height,
            //         truncation_index,
            //         boundary_extension,
            //     )?;
            // }
        }
    }

    // row filtering
    for (pole, truncation_index) in B::POLES
        .iter()
        .copied()
        .zip(truncation_indices.iter().copied())
    {
        image_data
            .chunks_exact_mut(width.get())
            .try_for_each(|row| {
                // NOTE PERF: interestingly, this performs better than the
                // slice based method in 2d benchmarks...
                apply_exponential_filter_in_place_1d(
                    pole,
                    row,
                    truncation_index,
                    boundary_extension,
                )
            })?;
    }

    if B::NORMALIZATION != B::Float::ONE {
        let normalization = B::NORMALIZATION.powi(2);
        image_data
            .iter_mut()
            .for_each(|val| *val = *val * normalization);
    }
    Ok(())
    // Normalization
}

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
// #[expect(dead_code)]
#[allow(unsafe_code)]
// #[cfg_attr(not(test), expect(dead_code))]
#[cfg_attr(not(test), expect(dead_code))]
pub(crate) unsafe fn apply_exponential_filter_in_place_strided<F>(
    // the pole, -1 < alpha < 0
    alpha: F,
    signal: *mut F,
    // NOTE(geo)... we either have to go full pointer here or we have to
    // split the 1D and 2D logic
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
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
{
    let mut signal =
        unsafe { StridedPointerAccess::new(signal, stride, count).ok_or(Error::Internal)? };

    if truncation_idx >= count.get() {
        return Err(Error::NotEnoughData);
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
        let signal_first = signal.first();
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

    signal.at_mut(0).write(s_lini);

    // causal filtering, see algo 3 in the B&M paper lines 2-4
    let mut last = signal.at(0);
    let mut idx = 1;
    let end_idx = count.get() - 1;
    if end_idx > signal.len() {
        return Err(Error::OutOfBounds);
    }

    loop {
        if idx > end_idx {
            break;
        }
        // let cur = signal.as_slice()[idx];
        let cur = signal.at(idx);
        let new = last.fmuladd(alpha, cur);
        signal.at_mut(idx).write(new);
        last = new;
        idx += 1;
    }

    // this is the same loop as above rewritten with iterators, which
    // makes things not necessarily clearer
    // signal
    //     .iter_mut()
    //     .skip(1)
    //     .scan(last, |last, cur| {
    //         let new = last.fmuladd(alpha, *cur);
    //         *cur = new;
    //         *last = *cur;
    //         Some(())
    //     })
    //     .for_each(drop);

    // line 5 in B&M algorithm 3, but we're not using eq (49) directly but the
    // special cases in eqs. (55)-(66)
    let s_lend = {
        let signal_last = signal.last();
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
                if signal.len() < 2 {
                    return Err(Error::NotEnoughData);
                }
                let signal_2nd_to_last = signal.at(signal.len() - 2);
                alpha / (alpha * alpha - F::ONE) * (signal_last + alpha * signal_2nd_to_last)
            }
        }
    };

    let l_end = signal.len() - 1;
    signal.at_mut(l_end).write(s_lend);

    // for idx in (0..l_end).rev() {
    //     signal[idx] = alpha * (signal[idx + 1] - signal[idx]);
    // }
    // anti-causal filter: signal[idx] = alpha * (signal[idx+1] - signal[idx])
    // let mut last = signal[l_end];
    // for idx in (0..l_end).rev() {
    //     let cur = signal[idx];
    //     let new = alpha * (last - cur);
    //     signal[idx] = new;
    //     last = new;
    // }

    let mut idx = end_idx - 1;
    let mut last = signal.at(end_idx);
    loop {
        let cur = signal.at(idx);
        let new = alpha * (last - cur);
        signal.at_mut(idx).write(new);
        last = new;
        if idx == 0 {
            break;
        }
        idx -= 1;
    }

    // signal
    //     .iter_mut()
    //     .rev()
    //     .skip(1)
    //     .scan(last, |last, cur| {
    //         let new = alpha * (*last - *cur);
    //         *cur = new;
    //         *last = *cur;
    //         Some(())
    //     })
    //     .for_each(drop);
    Ok(())
}

//TODO doc: apply the filter in place from elements {0, len-1}
// this only works for 1D signals becaue it requires the signal
// to be contiguous.
pub(crate) fn apply_exponential_filter_in_place_1d<F>(
    // the pole, -1 < alpha < 0
    alpha: F,
    signal: &mut [F],
    // NOTE: like in the Briand and Monasse paper this is the truncation INDEX,
    // not the number of elements to be truncated. The element at index N is
    // counted for the initialization.
    truncation_idx: usize,
    boundary_extension: TransmittableBoundaryExtension,
) -> Result<(), Error>
where
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
{
    if signal.is_empty() {
        return Err(Error::NoData);
    }

    if truncation_idx >= signal.len() {
        return Err(Error::NotEnoughData);
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
    //
    // eq. (39) in the Briand and Monasse paper, cf also algo 3, line 3.
    // NOTE PERF: here's some godbolt: https://godbolt.org/z/vrv5TGT18

    // NOTE PERF: this loop might perform very very slighly faster
    // than the iterator based implementation but not enough for me to abandon
    // the iterators

    // // --- original loop ---
    // let mut last = signal[0];
    // let mut idx = 1;
    // let end_idx = signal.len() - 1;
    // if end_idx >= signal.len() {
    //     return Err(Error::OutOfBounds);
    // }

    // loop {
    //     if idx > end_idx {
    //         break;
    //     }
    //     // let cur = signal.as_slice()[idx];
    //     let cur = signal[idx];
    //     let new = last.fmuladd(alpha, cur);
    //     signal[idx] = new;
    //     last = new;
    //     idx += 1;
    // }

    // --- level-2 algebraic unrolling ---
    let mut last = signal[0];
    let mut idx = 2;
    let end_idx = signal.len() - 1;
    if end_idx >= signal.len() {
        return Err(Error::OutOfBounds);
    }
    let alpha_sq = alpha * alpha;

    loop {
        if idx > end_idx {
            break;
        }
        let x0 = signal[idx];
        let x1 = signal[idx - 1];
        // y[idx-1] = x[idx-1] + alpha*y[idx-2];
        //          = x1 + alpha * last
        let y1 = last.fmuladd(alpha, x1);
        // y[idx] = x[idx] + alpha*y[idx-1] = x[idx] + alpha*x[idx-1] + alpha^2 * y[idx-2]
        //        = x0 + alpha*x1 + alpha^2 * last
        let y0 = last.fmuladd(alpha_sq, x1.fmuladd(alpha, x0));
        last = y0;
        signal[idx] = y0;
        signal[idx - 1] = y1;
        idx += 2;
    }
    // tail must be processed as well when unrolling! The last up to step-1 elements
    // must be processed, which for step 2 is at most 1 element
    if idx - 1 == end_idx {
        signal[end_idx] = signal[end_idx] + alpha * signal[end_idx - 1];
    }

    // // Iterator version
    // // this is the same loop as above rewritten with iterators, which
    // // makes things not necessarily clearer
    // let end_idx = signal.len() - 1;
    // let last = signal[0];
    // signal
    //     .iter_mut()
    //     .skip(1)
    //     .scan(last, |last, cur| {
    //         let new = last.fmuladd(alpha, *cur);
    //         *cur = new;
    //         *last = *cur;
    //         Some(())
    //     })
    //     .for_each(drop);

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

    signal[end_idx] = s_lend;

    // NOTE PERF: this is the original loop that I replaced with iterators below.
    // It performs maybe very very slighly better than the iterator based loop,
    // but not enough for me to abandon the iterator loop.

    // // --- original loop ---
    // let mut idx = end_idx - 1;
    // let mut last = signal[end_idx];
    // loop {
    //     let cur = signal[idx];
    //     let new = alpha * (last - cur);
    //     signal[idx] = new;
    //     last = new;
    //     if idx == 0 {
    //         break;
    //     }
    //     idx -= 1;
    // }

    // TODO TAIL HANDLING
    // --- algebraically unrolled ---
    const STEP: usize = 2;
    let mut idx;
    if end_idx > STEP {
        idx = end_idx - STEP;
        let mut last = signal[end_idx];
        loop {
            let x0 = signal[idx];
            let x1 = signal[idx + 1];
            let y1 = alpha * (last - x1);
            // NOTE PERF: fmuladd doesn't much help here... in my benchmarks
            // let y1 = last.fmuladd(alpha, neg_alpha * x1);
            // let y0 = alpha * (y1 - x0);
            let y0 = alpha_sq * last - alpha_sq * x1 - alpha * x0;
            // let y0 = alpha_sq * (last - x1) - alpha * x0;
            // let y0 = alpha * (alpha * (last - x1) - x0);
            // NOTE PERF: fmuladd doesn't much help here... in my benchmarks
            // let y0 = x0.fmuladd(neg_alpha, last.fmuladd(alpha_sq, neg_alpha_sq * x1));
            signal[idx + 1] = y1;
            signal[idx] = y0;
            last = y0;
            if idx < STEP {
                break;
            }
            idx -= STEP;
        }
    } else {
        idx = end_idx - 1;
    }
    if idx != 0 {
        signal[0] = alpha * (signal[1] - signal[0]);
    }

    // // --- iterator based loop ---
    // let last = signal[end_idx];
    // signal
    //     .iter_mut()
    //     .rev()
    //     .skip(1)
    //     .scan(last, |last, cur| {
    //         let new = alpha * (*last - *cur);
    //         *cur = new;
    //         *last = *cur;
    //         Some(())
    //     })
    //     .for_each(drop);

    Ok(())
}

pub(crate) fn apply_exponential_filter_in_place_1d_iter<S, F>(
    // the pole, -1 < alpha < 0
    alpha: F,
    mut signal: S,
    // NOTE: like in the Briand and Monasse paper this is the truncation INDEX,
    // not the number of elements to be truncated. The element at index N is
    // counted for the initialization.
    truncation_idx: usize,
    boundary_extension: TransmittableBoundaryExtension,
) -> Result<(), Error>
where
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
    S: crate::iter_utils::StridedRandomAccess<Item = F>,
{
    // use crate::iter_utils::StridedSlice;
    if signal.strided_is_empty() {
        return Err(Error::NoData);
    }

    if truncation_idx >= signal.strided_len() {
        return Err(Error::NotEnoughData);
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
        let signal_first = signal.strided_first().ok_or(Error::NoData)?;
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
                for s in signal.strided_iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
            TransmittableBoundaryExtension::WholeSymmetric => {
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.strided_iter().copied().skip(1).take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
            TransmittableBoundaryExtension::Periodic => {
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.strided_iter_rev().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
        }
    };

    // signal[0] = s_lini;
    signal.strided_write_first(s_lini);

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
    //
    // eq. (39) in the Briand and Monasse paper, cf also algo 3, line 3.
    // NOTE PERF: here's some godbolt: https://godbolt.org/z/vrv5TGT18

    // NOTE PERF: this loop might perform very very slighly faster
    // than the iterator based implementation but not enough for me to abandon
    // the iterators
    // let mut last = signal[0];
    // let mut idx = 1;
    // let end_idx = signal.len() - 1;
    // if end_idx >= signal.len() {
    //     return Err(Error::OutOfBounds);
    // }

    // loop {
    //     if idx > end_idx {
    //         break;
    //     }
    //     // let cur = signal.as_slice()[idx];
    //     let cur = signal[idx];
    //     let new = last.fmuladd(alpha, cur);
    //     signal[idx] = new;
    //     last = new;
    //     idx += 1;
    // }

    // this is the same loop as above rewritten with iterators, which
    // makes things not necessarily clearer
    let last = signal.strided_first().ok_or(Error::NotEnoughData)?;
    signal
        .strided_iter_mut()
        .skip(1)
        .scan(last, |last, cur| {
            let new = last.fmuladd(alpha, *cur);
            *cur = new;
            *last = *cur;
            Some(())
        })
        .for_each(drop);

    // line 5 in B&M algorithm 3, but we're not using eq (49) directly but the
    // special cases in eqs. (55)-(66)
    let s_lend = {
        let signal_end = signal.strided_last().ok_or(Error::NoData)?;
        match boundary_extension {
            TransmittableBoundaryExtension::Periodic => {
                let mut acc = signal_end;
                let mut pow = F::ONE;
                for s in signal.strided_iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                -alpha * acc
            }
            TransmittableBoundaryExtension::Mirrored => alpha / (alpha - F::ONE) * signal_end,
            TransmittableBoundaryExtension::WholeSymmetric => {
                let signal_2nd_to_last = signal.strided_2nd_last().ok_or(Error::NotEnoughData)?;
                alpha / (alpha * alpha - F::ONE) * (signal_end + alpha * signal_2nd_to_last)
            }
        }
    };

    signal.strided_write_last(s_lend);

    // NOTE PERF: this is the original loop that I replaced with iterators below.
    // It performs maybe very very slighly better than the iterator based loop,
    // but not enough for me to abandon the iterator loop.
    // let mut idx = end_idx - 1;
    // let mut last = signal[end_idx];

    // loop {
    //     let cur = signal[idx];
    //     let new = alpha * (last - cur);
    //     signal[idx] = new;
    //     last = new;
    //     if idx == 0 {
    //         break;
    //     }
    //     idx -= 1;
    // }

    let last = signal.strided_last().ok_or(Error::NoData)?;
    signal
        .strided_iter_mut_rev()
        .skip(1)
        .scan(last, |last, cur| {
            let new = alpha * (*last - *cur);
            *cur = new;
            *last = *cur;
            Some(())
        })
        .for_each(drop);

    Ok(())
}

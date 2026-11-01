use core::iter::Sum;
use num_traits::{ConstOne, ConstZero, Float};

use crate::Error;
use crate::RowMajor;
use crate::{
    bspline::BSplinePrefilter,
    inplace::InplaceBoundaryExtension,
    internal::compute_truncation_indices,
    iter_utils::{BlockedStridedSliceMut, StridedSliceMut},
    numerics::{FloatVector, FusedMuladd},
};

/// Calculate the spline coefficients for the given `bspline` _in place_. This
/// overwrites the input array with the coefficients necessary to perform
/// interpolation with the given B-Spline, see e.g.
/// [`interp2d_single`](crate::inplace::interpolation::interp2d_single). This step is also called _prefiltering_.
///
/// The function is a bit awkward to call because two of the generic arguments
/// (`ROW_BLOCKSIZE` and `COL_BLOCKSIZE`) must be explicitly given, which means
/// that the other arguments have to be written as well (but can be inferred
/// with `_`). This is why the [`Coefficients2D`](crate::inplace::Coefficients2d)
/// structure exists, which provides a builder style interface to this function
/// call just for ergonomics.
///
/// Note also that a more high level interface to the workflow of calculating
/// the coefficients and then using them for interpolation exists with
/// [`crate::inplace::SplineInterpolator2d`].
///
/// ## Arguments
///
/// The B-Spline order for interpolation is defined by `bspline`. The same B-Spline
/// order must be used both for coefficient calculation and interpolation.
///
/// The 2D image data must be given in row-major order in the `image_data`
/// argument, where the length of a row is the _width_ of the image.
/// Consequently, the number of rows (=the number of elements in a column) is
/// the height of the image.
///
/// To get a unique solution to the interpolation problem, this implementation
/// requires us to specify a continuation of the given function beyond the measured
/// grid, such as mirrored or periodic extention. This is what `boundary_extension`
/// is for. Though this library doesn't deal in extrapolation, this implicitly defines
/// an extrapolation beyond the boundaries.
///
/// The parameter `epsilon` controls how faithfully the library reproduces
/// the given function values _on the sampled grid points_. The interpolation
/// error on the grid points is less than `epsilon * max(abs(f))`, where `f`
/// is the image to interpolate (B&M sections 3.2 & 5.1). _Do not_ think of
/// `epsilon` as meaningfully controlling the interpolation accuracy between
/// grid points, because that depends primarily how well the spline model
/// can describe the sampled function at all.
///
/// ### Blocking
///
/// The function has two additional compile time (const generic) arguments
/// `ROW_BLOCKSIZE` and `COL_BLOCKSIZE` which allow us to specify a blocksize
/// for operations on the rows and columns. Those blocksizes can enable SIMD
/// accelaration for calculations, which _can_ make the calculations faster,
/// but the only way to find out on _your system_ is to measure. Block sizes
/// are in _numbers of elements_, not bytes, i.e. a block size of `8`
/// for `f64` data means 64 bytes.
///
/// ## Recommendations
///
/// While B&M don't give a general recommendation for chosing `epsilon`, I'll
/// state a couple helpful things here, based on their results (B&M section 6.2,
/// B&M remarks in appendix C.1).
///
/// * For `f32` precision, consider `epsilon = 1e-6`.
/// * Higher epsilon will slightly (but not dramatically) increase computational
///   complexity of the coefficient calculation.
/// * For typical use cases, the coefficient calculation is much faster than
///   the interpolation anyways, so giving a small epsilon here doesn't matter
///   in the grand scheme of things, unless you have to perform _many_
///   interpolations.
/// * Even for `f64`, due to numerical precision, `epsilon = 1e-12` is a
///   lower limit at which the relative precision guarantee holds exactly.
///   Lower epsilons will often still give lower absolute error, but the error
///   might not stricly be less than `epsilon*max(|f|)`.
///
/// In terms of `boundary_extension` it's hard to make general recommendations
/// and I'll just say that `scipy` uses the whole-symmetric extension in its
/// [`spline_filter`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter)
/// function.
///
/// ## References
///
/// * [b&m](https://www.ipol.im/pub/art/2018/221/) t. briand, p. monasse: "theory
///   and practice of image b-spline interpolation"
/// * [`scipy.ndimage.spline_filter`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter).
pub fn coefficients2d_blocked<
    B,
    const NPOLES: usize,
    const ROW_BLOCKSIZE: usize,
    const COL_BLOCKSIZE: usize,
>(
    bspline: B,
    image_data: RowMajor<&mut [B::Float]>,
    epsilon: B::Float,
    boundary_extension: InplaceBoundaryExtension,
) -> Result<(), Error>
where
    B::Float: Float + ConstOne + ConstZero + Sum,
    B: BSplinePrefilter<NPOLES>,
{
    const { assert!(COL_BLOCKSIZE >= 1) };
    let expected_len = match image_data
        .width
        .get()
        .overflowing_mul(image_data.height.get())
    {
        (val, false) => val,
        (_, true) => return Err(Error::Overflow),
    };

    if image_data.data.len() != expected_len {
        return Err(Error::Dimensions);
    }

    let width = image_data.width;
    let height = image_data.height;

    // NOTE: this is almost certainly not the most efficient way of doing things,
    // since the columns have massive stride typically,so a blocking strategy
    // is the way to go here in general

    let truncation_indices = compute_truncation_indices(bspline.poles(), epsilon)?;

    // column filtering
    // TODO: remove this horrible unsafe code. We should not have to do any of this
    // for 2D sequential row-by row iteration. We can write an abstraction that
    // safely returns column iterators for a row major mutable view of an array...
    // but for parallel iterators we have to use unsafe code. But maybe I
    // don't have to do any iterator stuff and I can use a blocking strategy
    // with manual indexing right away anyways, so let's see if we can get this
    // working and see how well it performs compared to the reference C impl

    for (pole, truncation_index) in bspline
        .poles()
        .iter()
        .copied()
        .zip(truncation_indices.iter().copied())
    {
        // process the columns in chunks of BLOCKSIZE
        // this if is probably not necessary, but since I'm doing all the
        // premature optimization I can think of, I might just as well
        // give the compiler the chance to remove that code at compile
        let remaining_start_idx = if COL_BLOCKSIZE > 1 {
            let mut idx = 0;
            for col in (0..width.get()).step_by(COL_BLOCKSIZE) {
                let Some(col) = BlockedStridedSliceMut::<COL_BLOCKSIZE, B::Float>::new(
                    image_data.data,
                    crate::iter_utils::StridedAccessArgs {
                        offset: col,
                        stride: width,
                        count: height,
                    },
                ) else {
                    break;
                };
                // .ok_or(Error::Dimensions)?;
                apply_exponential_filter_in_place_1d_iter_vector(
                    pole,
                    col,
                    truncation_index,
                    boundary_extension,
                )?;
                idx += COL_BLOCKSIZE;
            }
            idx
        } else {
            0
        };

        // Now process the rest of the columns one by one.
        // This might be optimized by trying to find the largest blocksize
        // by trial and error (at compile time by using unrolled ifs) that
        // can help here, but that's an optimization I wont do because it's
        // crazy... and probably current const generics aren't super helpful
        // here... but maybe I could use crunchy and unroll the loop by
        // halving the index not sure...
        for col in remaining_start_idx..width.get() {
            let col = StridedSliceMut::new(
                image_data.data,
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
        }
    }

    // row filtering
    for (pole, truncation_index) in bspline
        .poles()
        .iter()
        .copied()
        .zip(truncation_indices.iter().copied())
    {
        image_data
            .data
            .chunks_exact_mut(width.get())
            .try_for_each(|row| {
                // NOTE PERF: interestingly, this performs better than the
                // slice based method in 2d benchmarks...
                apply_exponential_filter_in_place_1d::<_, ROW_BLOCKSIZE>(
                    pole,
                    row,
                    truncation_index,
                    boundary_extension,
                )
            })?;
    }

    if B::normalization() != B::Float::ONE {
        let normalization = B::normalization().powi(2);
        image_data
            .data
            .iter_mut()
            .for_each(|val| *val = *val * normalization);
    }
    Ok(())
    // Normalization
}

//TODO doc: apply the filter in place from elements {0, len-1}
// this only works for 1D signals becaue it requires the signal
// to be contiguous.
pub(crate) fn apply_exponential_filter_in_place_1d<F, const BLOCKSIZE: usize>(
    // the pole, -1 < alpha < 0
    alpha: F,
    signal: &mut [F],
    // NOTE: like in the Briand and Monasse paper this is the truncation INDEX,
    // not the number of elements to be truncated. The element at index N is
    // counted for the initialization.
    truncation_idx: usize,
    boundary_extension: InplaceBoundaryExtension,
) -> Result<(), Error>
where
    F: Float + ConstOne + ConstZero + Sum + FusedMuladd,
{
    const { assert!(BLOCKSIZE > 0) };
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
            InplaceBoundaryExtension::Reflected => {
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
            InplaceBoundaryExtension::WholeSymmetric => {
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.iter().copied().skip(1).take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
            InplaceBoundaryExtension::Periodic => {
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

    // --- algebraic unrolling/unrolled ---
    let mut prologue = signal[0];
    // let mut idx = STEP;
    let end_idx = signal.len() - 1;
    if end_idx >= signal.len() {
        return Err(Error::OutOfBounds);
    }
    // let alpha_sq = alpha * alpha;

    // Always 1, independent of STEP!
    // NOTE: if I want to use aligned simd loads/stores in the future (which
    // also requires me to give the users the option of row padding), then
    // it's important to note that this starts at an unaligned index. I need
    // to process a number of head elements, before I'm aligned again. Index
    // 0 is aligned, but index 1 is not, index 7 might be aligned again.
    let mut start = 1;
    if BLOCKSIZE > 1 {
        // -- "generic" unrolled --
        while let Some(chunk) = signal[start..].first_chunk_mut::<BLOCKSIZE>() {
            *chunk = causal_blocked(*chunk, alpha, prologue);
            prologue = chunk[BLOCKSIZE - 1];
            start += BLOCKSIZE;
        }
    }

    // // tail must be processed as well when unrolling! The last up to step-1 elements
    // // must be processed, which for step 2 is at most 1 element
    // if idx - 1 == end_idx {
    //     signal[end_idx] = signal[end_idx] + alpha * signal[end_idx - 1];
    // }

    // tail handling.
    // NOTE PERF(geo) for some extra overengineering we could try and go through
    // the block sizes STEP-1, STEP-2, ..., 1 until the first one hits and
    // then use that block size.
    while start <= end_idx {
        let cur = signal[start];
        let new = prologue.fmuladd(alpha, cur);
        signal[start] = new;
        prologue = new;
        start += 1;
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
            InplaceBoundaryExtension::Periodic => {
                let mut acc = signal_last;
                let mut pow = F::ONE;
                for s in signal.iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                -alpha * acc
            }
            InplaceBoundaryExtension::Reflected => alpha / (alpha - F::ONE) * signal_last,
            InplaceBoundaryExtension::WholeSymmetric => {
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

    // // --- original loop (anticausal) ---
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

    // // TODO TAIL HANDLING
    // // --- manually algebraically unrolled ---
    // const STEP: usize = 2;
    // let mut idx;
    // if end_idx > STEP {
    //     idx = end_idx - STEP;
    //     let mut last = signal[end_idx];
    //     loop {
    //         let x0 = signal[idx];
    //         let x1 = signal[idx + 1];
    //         let y1 = alpha * (last - x1);
    //         // NOTE PERF: fmuladd doesn't much help here... in my benchmarks
    //         // let y1 = last.fmuladd(alpha, neg_alpha * x1);
    //         // let y0 = alpha * (y1 - x0);
    //         let y0 = alpha_sq * last - alpha_sq * x1 - alpha * x0;
    //         // let y0 = alpha_sq * (last - x1) - alpha * x0;
    //         // let y0 = alpha * (alpha * (last - x1) - x0);
    //         // NOTE PERF: fmuladd doesn't much help here... in my benchmarks
    //         // let y0 = x0.fmuladd(neg_alpha, last.fmuladd(alpha_sq, neg_alpha_sq * x1));
    //         signal[idx + 1] = y1;
    //         signal[idx] = y0;
    //         last = y0;
    //         if idx < STEP {
    //             break;
    //         }
    //         idx -= STEP;
    //     }
    // } else {
    //     idx = end_idx - 1;
    // }
    // if idx != 0 {
    //     signal[0] = alpha * (signal[1] - signal[0]);
    // }

    // TODO TAIL HANDLING
    // --- algebraically unrolled ---

    // idx is the first element of the block with `STEP` elements
    let mut idx = end_idx;
    let mut epilogue = signal[end_idx];
    if BLOCKSIZE > 1 {
        if end_idx > BLOCKSIZE {
            idx -= BLOCKSIZE;
            while let Some(chunk) = signal[idx..].first_chunk_mut::<BLOCKSIZE>() {
                *chunk = anticausal_blocked(*chunk, alpha, epilogue);
                epilogue = chunk[0];

                if idx < BLOCKSIZE {
                    break;
                }
                idx -= BLOCKSIZE;
            }
        }
    }

    // tail handling
    // NOTE(geo) for some (probably massive overengineering) we could revisit
    // the same note as in the causal block, where we use the
    if idx != 0 {
        idx -= 1;
        loop {
            let cur = signal[idx];
            let new = alpha * (epilogue - cur);
            signal[idx] = new;
            epilogue = new;
            if idx == 0 {
                break;
            }
            idx -= 1;
        }
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
    boundary_extension: InplaceBoundaryExtension,
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
            InplaceBoundaryExtension::Reflected => {
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
            InplaceBoundaryExtension::WholeSymmetric => {
                let mut acc = signal_first;
                let mut pow = F::ONE;
                for s in signal.strided_iter().copied().skip(1).take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                acc
            }
            InplaceBoundaryExtension::Periodic => {
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
            InplaceBoundaryExtension::Periodic => {
                let mut acc = signal_end;
                let mut pow = F::ONE;
                for s in signal.strided_iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    acc = acc + s * pow;
                }
                -alpha * acc
            }
            InplaceBoundaryExtension::Reflected => alpha / (alpha - F::ONE) * signal_end,
            InplaceBoundaryExtension::WholeSymmetric => {
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

pub(crate) fn apply_exponential_filter_in_place_1d_iter_vector<S, F>(
    // the pole, -1 < alpha < 0
    alpha: F::Scalar,
    mut signal: S,
    // NOTE: like in the Briand and Monasse paper this is the truncation INDEX,
    // not the number of elements to be truncated. The element at index N is
    // counted for the initialization.
    truncation_idx: usize,
    boundary_extension: InplaceBoundaryExtension,
) -> Result<(), Error>
where
    F: FloatVector,
    F::Scalar: num_traits::ConstOne + num_traits::ConstZero,
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
        -F::Scalar::ONE < alpha && alpha < F::Scalar::ZERO,
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
            InplaceBoundaryExtension::Reflected => {
                // NOTE: this is true for all the boundary code. I'm reusing the
                // same logic of recursively calculating a^j = a * a^(j-1) that
                // the reference C-code uses. If I use a^j = a.powi(j), I might
                // get better ILP (and actually better floating point accuracy
                // I think, but the code doesn't reproduce the float-identical
                // results anymore.
                let mut acc = signal_first;
                let mut pow = F::Scalar::ONE;
                for s in signal.strided_iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    // acc = acc + s * pow;
                    acc = s.axpy(pow, acc)
                }
                acc
            }
            InplaceBoundaryExtension::WholeSymmetric => {
                let mut acc = signal_first;
                let mut pow = F::Scalar::ONE;
                for s in signal.strided_iter().copied().skip(1).take(truncation_idx) {
                    pow = pow * alpha;
                    // acc = acc + s * pow;
                    acc = s.axpy(pow, acc);
                }
                acc
            }
            InplaceBoundaryExtension::Periodic => {
                let mut acc = signal_first;
                let mut pow = F::Scalar::ONE;
                for s in signal.strided_iter_rev().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    // acc = acc + s * pow;
                    acc = s.axpy(pow, acc);
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
            let new = last.axpy(alpha, *cur);
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
            InplaceBoundaryExtension::Periodic => {
                let mut acc = signal_end;
                let mut pow = F::Scalar::ONE;
                for s in signal.strided_iter().copied().take(truncation_idx) {
                    pow = pow * alpha;
                    // acc = acc + s * pow;
                    acc = s.axpy(pow, acc);
                }
                acc.scal(-alpha)
            }
            InplaceBoundaryExtension::Reflected => {
                signal_end.scal(alpha / (alpha - F::Scalar::ONE))
            }
            InplaceBoundaryExtension::WholeSymmetric => {
                let signal_2nd_to_last = signal.strided_2nd_last().ok_or(Error::NotEnoughData)?;
                signal_2nd_to_last
                    .axpy(alpha, signal_end)
                    .scal(alpha / (alpha * alpha - F::Scalar::ONE))
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
            let new = (last.sub(*cur)).scal(alpha);
            *cur = new;
            *last = *cur;
            Some(())
        })
        .for_each(drop);

    Ok(())
}

#[inline(always)]
pub(crate) fn causal_blocked<const N: usize, F: Float + Sum + ConstZero>(
    x: [F; N],
    alpha: F,
    prologue: F,
) -> [F; N] {
    const { assert!(N < i32::MAX as _) };
    const { assert!(N > 0) };
    // causal relation is eqn 39 in the B&M paper
    //
    // y[i] = x[i] + alpha*y[i-1]
    //
    // where y[-1] = p = prologue
    //
    // NOTE that the causal recurrence goes left to right and
    // y[0] is the left end of the block and prologue is before that.
    //
    // manually unrolling the recurrence relation also follows the left
    // to right logic, meaning
    //
    // y[0] = x[0] + alpha*p  [all of those terms are known]
    // y[1] = x[1] + alpha*y[0] = x[1] + alpha x[0] + alpha^2 * p
    // y[2] = ...
    //
    // y[k] = ...
    let mut y = [F::ZERO; N];

    // TODO PERF: investigate the pre-calculation based approach: is this really useful!!??
    // TODP PERF: if so, make this into a struct that's initialized once
    // and carries that state.
    // let mut alpha_powi = [F::ZERO; N];
    // for n in 0..N {
    //     alpha_powi[n] = alpha.powi(n as _);
    // }
    // let alpha_powi = alpha_powi;
    // for j in 0..N {
    //     for k in 0..=j {
    //         y[j] = x[j - k].mul_add(alpha_powi[k], y[j]);
    //     }

    //     //                         hacky because we need alpha_powi[j+1] but we don't have j+1 here
    //     y[j] = prologue.mul_add(alpha * alpha_powi[j], y[j]);
    // }
    // y

    for j in 0..N {
        for k in 0..=j {
            y[j] = x[j - k].fmuladd(alpha.powi(k as i32), y[j]);
        }

        y[j] = prologue.fmuladd(alpha.powi((j + 1) as i32), y[j]);
    }
    y
}

#[inline(always)]
pub(crate) fn anticausal_blocked<const N: usize, F: Float + Sum + ConstZero>(
    y: [F; N],
    alpha: F,
    epilogue: F,
) -> [F; N] {
    const { assert!(N < i32::MAX as _) };

    // anticausal relation is eqn 47 in the B&M paper
    // z[i] = alpha (z[i+1] - y[i])
    //
    // where the recurrence goes from RIGHT TO LEFT. And since
    // we have N block elements indexed [0,...,N-1]
    // that means our epilogue comes into play as:
    //
    // z[N-1 + 1] = z[N] = epilogue = e
    //
    // and then it makes sense to first write the recurrence
    // relation starting at the last element z[N-1]:
    //
    // z[N-1] = alpha (e - y[N-1]) = alpha*e - alpha*y[N-1]
    // z[N-1 -1] = alpha * (z[N-1] - y[N-1-1]) = alpha^2 * e - alpha^2 * y[N-1] - alpha y[N-1-1]
    //
    // this leads us to a formulation for z[k]:
    // z[k] = alpha^(N-k)*e - \sum_{j=k}^{N-1} (alpha^{j-k+1}y[j])

    // TODO perf, see if there are better ways of expressing the sums,
    // maybe by hand-coding the terms for certain block sizes rather than
    // using the loops here. Those can (in theory) be unrolled by the compiler,
    // but I'm not sure that it'll take advantage of this...
    let mut z = [F::ZERO; N];
    for k in 0..N {
        let mut sum = F::ZERO;
        for j in k..N {
            sum = y[j].fmuladd(alpha.powi((j - k + 1) as i32), sum);
        }

        z[k] = epilogue.fmuladd(alpha.powi((N - k) as i32), -sum);
    }

    z
}

use crate::Coord;
use crate::RowMajor;
use crate::{
    Error, bspline::BSplineWeights, inplace::InplaceBoundaryExtension, numerics::FloatVector,
};
use num_traits::{ConstOne, ConstZero, Float, FromPrimitive, NumCast, ToPrimitive};

#[cfg(test)]
mod test;

/// Interpolate a single point at floating point coordinate `coord` of the image.
///
/// **Attention (1)**: Using this function on its own requires passing the same
/// boundary extension and spline type that was used to compute the coefficients.
/// If you can use the [`SplineInterpolator2d`](crate::inplace::SplineInterpolator2d)
/// high level interface to this library, please do so as it provides a much more
/// convenient zero-cost interface that manages this mental load for you.
///
/// The function accepts the 2D B-Spline coefficients `coeffs2d` in row-major
/// order, that must have been calculated with
/// [`coefficients2d`](`crate::inplace::coefficients2d`) with a B-Spline of
/// type `B`, precision `epsilon`, and boundary extension `boundary_ext`.
///
/// ## Extrapolation
///
/// When coordinates outside of the image are given, the continuation implied
/// by `boundary_condition` is used to extrapolate the values beyond the grid.
/// While this makes the library internally consistent, it should be obvious
/// that this is _not useful_ for actual prediction. Don't use it for that!
///
/// ## Example
///
/// This assumes that the coefficients were already generate
///
/// ```rust
/// # use splinterpol::inplace::interp2d_single;
/// # use core::num::NonZeroUsize;
/// # use splinterpol::*;
/// # use splinterpol::bspline::*;
/// # let epsilon = 1e-1;
/// # let width = core::num::NonZero::new(9).unwrap();
/// # let height = core::num::NonZero::new(10).unwrap();
/// # use splinterpol::inplace::InplaceBoundaryExtension;
/// # fn foo(some_data:& [f64], width: NonZeroUsize, height: NonZeroUsize, boundary_ext: InplaceBoundaryExtension) {
/// let mut coeffs = RowMajor {
///   data: some_data.as_ref(),
///   width,
///   height
/// };
///
/// //                          must be the same type as was used
/// //                          for calculating the coefficients
/// //                           |
/// //                           |              always inferred, leave as _
/// //                           |               |
/// //                           v               v
/// let val = interp2d_single::<BSpline3::<f64>, _>(
///   coeffs,
///   // boundary extension must be the same that was
///   // used for computing coefficients
///   boundary_ext,
///   Coord::x(1.3).y(3.7)
/// ).unwrap();
///
/// # }
/// ```
#[inline]
pub fn interp2d_single<B, const NW: usize>(
    coeffs2d: RowMajor<&[B::Float]>,
    boundary_extension: InplaceBoundaryExtension,
    // NOTE: just keep in mind that `x` is the column-index and `y` is the row index
    coord: Coord<B::Float>,
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float:
        Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast + core::iter::Sum,
{
    let Coord { row, col } = coord;
    // NOTE: see the B::weights function. The left hand size index must be (x-R).ceil(), where
    // `x` is the floating point row or column index and R is the radius of the spline.
    let fc0 = (col - B::radius()).ceil();
    let fr0 = (row - B::radius()).ceil();
    let c0 = fc0.to_isize().ok_or(Error::Unrepresentable)?;
    let r0 = fr0.to_isize().ok_or(Error::Unrepresentable)?;
    let wx = B::weights(col - fc0);
    let wy = B::weights(row - fr0);

    // let patch = get_row_patch::<NW, _>(coeffs, width, height, ix, iy, boundary_extension);

    let row_indices: [usize; NW] = core::array::from_fn(|j| {
        let idx = r0.wrapping_add(j as isize);
        if (0..coeffs2d.height.get() as isize).contains(&idx) {
            idx as usize
        } else {
            boundary_extension.extended_index(idx, coeffs2d.height.get())
        }
    });

    // NOTE(geo): ON THE CALCULATIONS HERE
    //
    // So the interpolation works like this. To interpolate a value
    // at the given coordinates, we have to find the product
    //
    // f = w_y^T P w_x
    //
    // where `f` is the scalar interpolated value and w_x, w_y in R^NW are the weight
    // vectors above and P is a NW x NW matrix representing the coefficient
    // patch, which is just
    //
    //     | -- p_0^T -- |
    // P = | -- p_1^T -- |
    //     |     ...     |
    //
    // where p_j^T are the NW-vectors of the row patches corresponding to the
    // row indices above (0 -> row_indices[0], 1->row_indices[1],...)
    //
    // This is a vector-matrix-vector product. There are 2 ways I can think
    // of how to actually calculate this product in a step-wise fashion (that's
    // without completely unrolling it I guess).
    //
    // # VARIANT 1:
    //
    // Step 1: q^T = w_y^T P = w_y[0] * p_0^T + w_y[1] * p_1^T + ...
    // Step 2: f = q^T w_x
    //
    // Note that the vector q^T is formed via a scalar*vector sum. Only the
    // last step needs a dot product.
    //
    // # VARIANT 2:
    //
    // Step 1: v = P w_x
    // Step 2: f = w_y^T v
    //
    // where we can calculate v[j] = p_j^T w_x independenly but via dot products
    // and then use a final dot product.
    //
    // NOTE(PERF)
    //
    // The performance results are a bit inconclusive here and warrant some
    // further investigation. On my machine (TM), I can see that variant 2
    // beats variant 1 on larger images (4096x4096) by roughly 5% but shows
    // a regression by about 2% on the small benchmarked images (128x128).

    // // # VARIANT 1 GROUPING
    // //
    // let mut sum = [B::Float::ZERO; NW];
    // for j in 0..NW {
    //     let r_j = row_indices[j];
    //     let wy_j = wy[j];
    //     let start = coeffs2d.width.get() * r_j;
    //     let row = coeffs2d
    //         .data
    //         .get(start..start.wrapping_add(coeffs2d.width.get()))
    //         .ok_or(Error::OutOfBounds)?;
    //     let coeffs = get_row_patch::<NW, _>(row, c0, boundary_extension);
    //     sum = coeffs.axpy(wy_j, sum);
    // }
    // Ok(sum.dot(wx))

    // # VARIANT 2 GROUPING
    let mut row_values = [B::Float::ZERO; NW];
    for j in 0..NW {
        let r_j = row_indices[j];
        let start = coeffs2d.width.get() * r_j;
        let row = coeffs2d
            .data
            .get(start..start.wrapping_add(coeffs2d.width.get()))
            .ok_or(Error::OutOfBounds)?;
        let coeffs = get_row_patch::<NW, _>(row, c0, boundary_extension);
        row_values[j] = coeffs.dot(wx);
    }

    Ok(wy.dot(row_values))
}

#[inline]
/// This works just like [`interp2d_single`], except that it works in 1D. That
/// means the coefficients need to have been generated with [`coefficients1d`](crate::inplace::coefficients1d)
/// using the same spline type `B` and the same `boundary_extension`.
///
/// ```rust
/// # use splinterpol::inplace::interp1d_single;
/// # use core::num::NonZeroUsize;
/// # use splinterpol::*;
/// # use splinterpol::bspline::*;
/// # let epsilon = 1e-1;
/// # let width = core::num::NonZero::new(9).unwrap();
/// # let height = core::num::NonZero::new(10).unwrap();
/// # use splinterpol::inplace::InplaceBoundaryExtension;
/// # fn foo(coeffs:& [f64], width: NonZeroUsize, height: NonZeroUsize, boundary_ext: InplaceBoundaryExtension) {
/// //                          must be the same type as was used
/// //                          for calculating the coefficients
/// //                           |
/// //                           |              always inferred, leave as _
/// //                           |               |
/// //                           v               v
/// let val = interp1d_single::<BSpline3::<f64>, _>(
///   coeffs,
///   // boundary extension must be the same that was
///   // used for computing coefficients
///   boundary_ext,
///   12.34
/// ).unwrap();
///
/// # }
/// ```
pub fn interp1d_single<B, const NW: usize>(
    coeffs: &[B::Float],
    boundary_extension: InplaceBoundaryExtension,
    pos: B::Float,
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float:
        Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast + core::iter::Sum,
{
    let fx0 = (pos - B::radius()).ceil();
    let x0 = fx0.to_isize().ok_or(Error::Unrepresentable)?;
    let w = B::weights(pos - fx0);
    let c = get_row_patch::<NW, _>(coeffs, x0, boundary_extension);

    Ok(c.dot(w))
}

#[inline(always)]
#[allow(clippy::unwrap_used)]
fn get_row_patch<const N: usize, F: Float + ConstZero>(
    row: &[F],
    idx: isize,
    boundary_extension: InplaceBoundaryExtension,
) -> [F; N] {
    const { assert!(N > 0 && N < isize::MAX as usize) };
    debug_assert!(
        !row.is_empty(),
        "getting a row patch from an empty row points to a logic error"
    );
    // TODO PERF check this

    // Fast path: entire patch is inside the row.
    if let Ok(start) = usize::try_from(idx)
        && let Some(patch) = row.get(start..).and_then(|slice| slice.first_chunk::<N>())
    {
        return *patch;
    }

    // NOTE(geo) we use wrapping adds here, which will produce incorrect
    // (but not unsafe) results for very pathological indices, aber as we
    // say in Germany einen Tod muss man sterben

    // Slow path: one or both boundaries are crossed.
    let mut patch = [F::ZERO; N];
    if row.is_empty() {
        return patch;
    }

    for j in 0..N {
        let i = idx.wrapping_add(j as isize);
        if let Some(value) = row.get(i as usize).copied() {
            patch[j] = value;
        } else {
            patch[j] = row[boundary_extension.extended_index(i, row.len())];
        }
    }
    patch
}

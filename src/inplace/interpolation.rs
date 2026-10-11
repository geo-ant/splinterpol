use crate::{
    Error, bspline::BSplineWeights, inplace::TransmittableBoundaryExtension, numerics::FloatVector,
};
use num_traits::{ConstOne, ConstZero, Float, FromPrimitive, NumCast, ToPrimitive};

#[cfg(test)]
mod test;

#[inline]
pub fn interp2d_single<const NW: usize, B>(
    bspline: B,
    coeffs2d: &[B::Float],
    width: usize,
    height: usize,
    boundary_extension: TransmittableBoundaryExtension,
    // NOTE: just keep in mind that `x` is the column-index and `y` is the row index
    (rowpos, colpos): (B::Float, B::Float),
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float: Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast,
{
    let fc0 = (colpos - B::radius()).ceil();
    let fr0 = (rowpos - B::radius()).ceil();
    let c0 = fc0.to_isize().ok_or(Error::UnrepresentableNumerically)?;
    let r0 = fr0.to_isize().ok_or(Error::UnrepresentableNumerically)?;
    let wx = compute_weights(bspline, colpos - fc0).ok_or(Error::UnrepresentableNumerically)?;
    let wy = compute_weights(bspline, rowpos - fr0).ok_or(Error::UnrepresentableNumerically)?;

    // let patch = get_row_patch::<NW, _>(coeffs, width, height, ix, iy, boundary_extension);

    let row_indices: [usize; NW] = core::array::from_fn(|j| {
        let idx = r0.wrapping_add(j as isize);
        if (0..height as isize).contains(&idx) {
            idx as usize
        } else {
            boundary_extension.extended_index(idx, height)
        }
    });

    let mut sum = [B::Float::ZERO; NW];
    for j in 0..NW {
        let r_j = row_indices[j];
        let wy_j = wy[j];
        let start = width * r_j;
        let row = coeffs2d
            .get(start..start.wrapping_add(width))
            .ok_or(Error::OutOfBounds)?;
        let coeffs = get_row_patch::<NW, _>(row, c0, boundary_extension);
        sum = coeffs.axpy(wy_j, sum);
    }

    Ok(sum.dot(wx))
}

#[inline]
pub fn interp1d_single<const NW: usize, B>(
    bspline: B,
    coeffs: &[B::Float],
    pos: B::Float,
    boundary_extension: TransmittableBoundaryExtension,
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float: Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast,
{
    let fx0 = (pos - B::radius()).ceil();
    let x0 = fx0.to_isize().ok_or(Error::UnrepresentableNumerically)?;
    let w = compute_weights(bspline, pos - fx0).ok_or(Error::OutOfBounds)?;
    let c = get_row_patch::<NW, _>(coeffs, x0, boundary_extension);

    Ok(c.dot(w))
}

#[inline(always)]
fn compute_weights<const NW: usize, B: BSplineWeights<NW>>(
    bspline: B,
    t: B::Float,
) -> Option<[B::Float; NW]>
where
    B::Float: Float + ConstOne + ConstZero + NumCast,
{
    // TODO PERF: check if this is optimized, maybe use unsafe code here
    // Maybe we can also use some casting code that makes sure at compile
    // time that the integers can indeed be casted..., but the num_traits
    // crate doesn't offer that out of the box, but maybe this is optimized
    // anyways
    let mut weights = [B::Float::ZERO; NW];
    for j in 0..NW {
        weights[j] = bspline.eval(t - NumCast::from(j)?);
    }
    Some(weights)
}

#[inline]
#[allow(clippy::unwrap_used)]
fn get_row_patch<const N: usize, F: Float + ConstZero>(
    row: &[F],
    idx: isize,
    boundary_extension: TransmittableBoundaryExtension,
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

use crate::{
    Error, bspline::BSplineWeights, inplace::TransmittableBoundaryExtension, numerics::FloatVector,
};
use num_traits::{ConstOne, ConstZero, Float, FromPrimitive, NumCast, ToPrimitive};

#[inline]
pub fn interp2d_single<const NW: usize, B>(
    bspline: B,
    coeffs: &[B::Float],
    width: usize,
    height: usize,
    boundary_extension: TransmittableBoundaryExtension,
    // NOTE: just keep in mind that `x` is the column-index and `y` is the row index
    (x, y): (B::Float, B::Float),
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float: Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast,
{
    let xoffset = (x - B::radius()).ceil();
    let yoffset = (y - B::radius()).ceil();
    let ix = xoffset.to_isize().ok_or(Error::UnrepresentableAsFloat)?;
    let iy = yoffset.to_isize().ok_or(Error::UnrepresentableAsFloat)?;
    let wx = compute_weights(bspline, x - xoffset).ok_or(Error::UnrepresentableAsFloat)?;
    let wy = compute_weights(bspline, y - yoffset).ok_or(Error::UnrepresentableAsFloat)?;

    // let patch = get_row_patch::<NW, _>(coeffs, width, height, ix, iy, boundary_extension);

    let row_indices: [usize; NW] = core::array::from_fn(|j| {
        let idx = iy.wrapping_add(j as isize);
        if (0..height as isize).contains(&idx) {
            idx as usize
        } else {
            boundary_extension.extended_index(idx, height)
        }
    });

    let mut sum = [B::Float::ZERO; NW];
    for (r_j, wy_j) in row_indices.into_iter().zip(wy) {
        let row = coeffs
            .get(width * r_j..width * (r_j + 1))
            .ok_or(Error::OutOfBounds)?;
        let coeffs = get_row_patch::<NW, _>(row, ix, boundary_extension);
        sum = coeffs.axpy(wy_j, sum);
    }

    Ok(sum.dot(wx))
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
fn get_row_patch<const N: usize, F: Float>(
    row: &[F],
    idx: isize,
    boundary_extension: TransmittableBoundaryExtension,
) -> [F; N] {
    const { assert!(N > 0) };
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
    core::array::from_fn(|j| {
        let i = idx.wrapping_add(j as isize);

        if let Some(&value) = row.get(i as usize) {
            value
        } else {
            row[boundary_extension.extended_index(i, row.len())]
        }
    })
}

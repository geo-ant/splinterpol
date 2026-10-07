use crate::{Error, bspline::BSplineWeights, numerics::FloatVector};
use num_traits::{ConstOne, ConstZero, Float, NumCast};

pub fn interpolate_one<const NW: usize, B>(
    bspline: B,
    coeffs: &[B::Float],
    width: usize,
    height: usize,
    coord: B::Float,
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float: Float + ConstOne + ConstZero + NumCast,
{
    let foo = [B::Float::ZERO; NW];
    let offset = B::Float::ONE;
    let weights = compute_weights(bspline, offset).ok_or(Error::UnrepresentableAsFloat)?;
    foo.axpy(B::Float::ONE, weights);
    todo!()
}

#[inline(always)]
fn compute_weights<const NW: usize, B: BSplineWeights<NW>>(
    bspline: B,
    offset: B::Float,
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
        weights[j] = bspline.eval(offset + NumCast::from(j)?);
    }
    Some(weights)
}

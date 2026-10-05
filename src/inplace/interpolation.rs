use crate::{Error, bspline::BSplineWeights, numerics::FloatVector};
use num_traits::{ConstOne, ConstZero, Float};

pub fn interpolate_one<const NW: usize, B>(
    bspline: B,
    coeffs: &[B::Float],
    width: usize,
    height: usize,
    coord: B::Float,
) -> Result<B::Float, Error>
where
    B: BSplineWeights<NW>,
    B::Float: Float + ConstOne + ConstZero,
{
    let foo = [B::Float::ZERO; NW];
    let bar = bspline.weights(coord);
    foo.axpy(B::Float::ONE, bar);
    todo!()
}

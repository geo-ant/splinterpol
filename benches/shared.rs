use num_traits::{Float, FromPrimitive};
use rand::{
    RngExt, SeedableRng,
    distr::uniform::{SampleRange, SampleUniform},
    rngs::{ThreadRng, Xoshiro256PlusPlus},
};
use std::ops::Range;
pub fn make_random_data<F: Float + FromPrimitive + SampleUniform>(
    len: usize,
    range: Range<F>,
    seed: Option<u64>,
) -> Vec<F>
where
    Range<F>: SampleRange<F>,
{
    let mut rng = match seed {
        Some(seed) => Xoshiro256PlusPlus::seed_from_u64(seed),
        None => {
            let mut rng = ThreadRng::default();
            Xoshiro256PlusPlus::from_rng(&mut rng)
        }
    };
    (0..len).map(|_| rng.random_range(range.clone())).collect()
}

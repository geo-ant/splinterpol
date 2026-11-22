use num_traits::{Float, FromPrimitive};
use rand::{
    RngExt, SeedableRng,
    distr::uniform::{SampleRange, SampleUniform},
    rngs::{ThreadRng, Xoshiro256PlusPlus},
};
use std::ops::Range;

// make random data with a configurable seed
pub fn make_random_data<F: Float + FromPrimitive + SampleUniform>(
    len: usize,
    range: Range<F>,
    // if the seed is not specified in the environment, we take this one.
    // if this and the environment seed is not specified we pick one at
    // (pseudo)random
    rng_seed: Option<u64>,
) -> Vec<F>
where
    Range<F>: SampleRange<F>,
{
    // environment overrides the seed even if given
    let env_seed = std::env::var("SP_BENCH_SEED").ok().map(|val| {
        let num: u64 = val.parse().expect("random seed must be a u64 number");
        num
    });

    let seed = {
        if let Some(s) = env_seed {
            Some(s)
        } else {
            rng_seed
        }
    };

    let mut rng = match seed {
        Some(seed) => Xoshiro256PlusPlus::seed_from_u64(seed),
        None => {
            let mut rng = ThreadRng::default();
            Xoshiro256PlusPlus::from_rng(&mut rng)
        }
    };
    (0..len).map(|_| rng.random_range(range.clone())).collect()
}

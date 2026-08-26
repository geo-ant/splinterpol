use core::num::NonZero;
use core::ops::Range;
use criterion::{Criterion, criterion_group, criterion_main};
use num_traits::{Float, FromPrimitive};
use rand::{
    RngExt, SeedableRng,
    distr::uniform::{SampleRange, SampleUniform},
    rngs::{ThreadRng, Xoshiro256PlusPlus},
};
use spleen::bench_exports::apply_exponential_filter_in_place;

fn make_random_data<F: Float + FromPrimitive + SampleUniform>(
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

fn bench_exponential_filter(c: &mut Criterion) {
    // pristine source data; each iteration filters a fresh copy so we always
    // measure filtering the same unfiltered input (not re-filtering smoothed data)
    let data = make_random_data(65535, -100f64..100f64, Some(12345678910));

    let alpha = -0.28;
    c.bench_function("foo", move |b| {
        let data_len = data.len();
        b.iter_batched(
            // setup: produce a fresh, unfiltered copy (not timed)
            || data.clone(),
            // routine: filter the fresh copy (timed)
            |mut fresh| {
                apply_exponential_filter_in_place(
                    alpha,
                    &mut fresh,
                    0,
                    NonZero::new(1).unwrap(),
                    NonZero::new(data_len).unwrap(),
                    1000,
                    spleen::TransmittableBoundaryExtension::Periodic,
                )
                .unwrap()
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, bench_exponential_filter);
criterion_main!(benches);

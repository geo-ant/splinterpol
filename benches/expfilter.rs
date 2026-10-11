use core::num::NonZero;
use core::ops::Range;
use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main};
use num_traits::{Float, FromPrimitive};
use rand::{
    RngExt, SeedableRng,
    distr::uniform::{SampleRange, SampleUniform},
    rngs::{ThreadRng, Xoshiro256PlusPlus},
};
use spleen::{
    bench_exports::apply_exponential_filter_in_place, bench_exports::compute_truncation_indices,
    inplace::TransmittableBoundaryExtension,
};

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
    let sizes = [100, 128, 1000, 1024, 10000, 16384, 60000, 65536];
    // this is a small pole taken from the actual spline implementations
    let pole = -1.698276282327549e-5;
    // pristine source data; each iteration filters a fresh copy so we always
    // measure filtering the same unfiltered input (not re-filtering smoothed data)
    let data = make_random_data(
        sizes.iter().max().copied().unwrap(),
        -100f64..100f64,
        Some(12345678910),
    );

    let mut group = c.benchmark_group("exp-filter");
    let epsilon = 1e-8;

    let trunc_idx = compute_truncation_indices(&[pole], epsilon).unwrap()[0];

    for size in sizes {
        for ext in [
            TransmittableBoundaryExtension::Mirrored,
            TransmittableBoundaryExtension::WholeSymmetric,
            TransmittableBoundaryExtension::Periodic,
        ] {
            group
                .bench_with_input(
                    BenchmarkId::new(format!("unblocked/{:?}", ext), size),
                    &data,
                    |b, data| {
                        b.iter_batched(
                            || data[0..size].to_vec(),
                            |mut data| {
                                let len = data.len();
                                apply_exponential_filter_in_place::<_, 1>(
                                    std::hint::black_box(pole),
                                    data.as_mut_slice(),
                                    trunc_idx,
                                    TransmittableBoundaryExtension::Periodic,
                                )
                                .unwrap();
                                std::hint::black_box(data)
                            },
                            criterion::BatchSize::LargeInput,
                        );
                    },
                )
                .throughput(criterion::Throughput::Elements(size as _));

            group
                .bench_with_input(
                    BenchmarkId::new(format!("blocked-8/{:?}", ext), size),
                    &data,
                    |b, data| {
                        b.iter_batched(
                            || data[0..size].to_vec(),
                            |mut data| {
                                let len = data.len();
                                apply_exponential_filter_in_place::<_, 8>(
                                    std::hint::black_box(pole),
                                    data.as_mut_slice(),
                                    len / 5,
                                    TransmittableBoundaryExtension::Periodic,
                                )
                                .unwrap();
                                std::hint::black_box(data)
                            },
                            criterion::BatchSize::LargeInput,
                        );
                    },
                )
                .throughput(criterion::Throughput::Elements(size as _));
        }
    }
}

criterion_group!(benches, bench_exponential_filter);
criterion_main!(benches);

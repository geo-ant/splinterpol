use core::num::NonZero;
use core::ops::Range;
use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main};
use num_traits::{Float, FromPrimitive};
use splinterpol::{
    bench_exports::apply_exponential_filter_in_place, bench_exports::compute_truncation_indices,
    inplace::InplaceBoundaryExtension,
};

use shared::make_random_data;
mod shared;

fn bench_exponential_filter(c: &mut Criterion) {
    let sizes = [127, 128, 129, 4095, 4096, 4097];
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
        group.throughput(criterion::Throughput::Elements(size as u64));
        for ext in [
            InplaceBoundaryExtension::Reflected,
            InplaceBoundaryExtension::WholeSymmetric,
            InplaceBoundaryExtension::Periodic,
        ] {
            group.bench_with_input(
                BenchmarkId::new(format!("unblocked/{:?}", ext), size),
                &data,
                |b, data| {
                    b.iter_batched(
                        || data[0..size].to_vec(),
                        |mut data| {
                            apply_exponential_filter_in_place::<_, 1>(
                                std::hint::black_box(pole),
                                data.as_mut_slice(),
                                trunc_idx,
                                ext,
                            )
                            .unwrap();
                            std::hint::black_box(data)
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );

            group.bench_with_input(
                BenchmarkId::new(format!("blocked-8/{:?}", ext), size),
                &data,
                |b, data| {
                    b.iter_batched(
                        || data[0..size].to_vec(),
                        |mut data| {
                            apply_exponential_filter_in_place::<_, 8>(
                                std::hint::black_box(pole),
                                data.as_mut_slice(),
                                trunc_idx,
                                ext,
                            )
                            .unwrap();
                            std::hint::black_box(data)
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );
        }
    }
}

criterion_group!(benches, bench_exponential_filter);
criterion_main!(benches);

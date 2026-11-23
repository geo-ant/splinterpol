use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use interpn::{
    MultiBsplineRegular,
    multibspline::regular::{coefficients, coefficients_par},
};

use shared::make_random_data;
use splinterpol::{
    ColBlockSize, RowBlockSize,
    bspline::BSpline3,
    inplace::{InplaceBoundaryExtension::*, coefficients1d, coefficients2d},
};
mod shared;

// NOTE: those sizes are for a small dataset and a large dataset. In the 2D case,
// those sizes are squared (width=height=size). We're using those weird clusters
// because I've observed very strange effects in both libraries where the
// 4096x4096 image is actually slower than the 4097x4097 by a significant
// amount ON MY MACHINE (TM). I'm suspecting some really unlucky cache eviction
// behavior. This affects both splinterpol and interpn and maybe (??) suggests
// that padding the data could have an advantage. Not sure, have to investigate.
const SIZES: [usize; 6] = [127usize, 128, 129, 4095, 4096, 4097];

fn bench_coefficients_1d(c: &mut Criterion) {
    // creating the largest images size once and then "clippin" the images
    // saves us having to regenerate the values.
    let max_len = *SIZES.iter().max().unwrap();
    let data = make_random_data(max_len, -100f64..100f64, Some(313373));

    let mut group = c.benchmark_group("Coefficients-1D");
    // see also the comments in the tests. While smaller values of epsilon will
    // lead to higher precision, the hard guarantee on the relative error being
    // less than eps*max(|f|) doesn't hold anymore due to numerical errors.
    // So this is the strictest reasonable value.
    let epsilon = 1e-12;

    for len in SIZES {
        group.throughput(Throughput::Elements(len as u64));

        // interpn get's to reuse the cache and scratch in the hope that this
        // will produce a benchmark that is most favorable to interpn, since
        // it's not my crate and I don't want to be biased.
        let dims = [len];
        let mut coeffs = vec![0.0; MultiBsplineRegular::<f64, 1>::coeff_storage_len(dims)];
        let mut scratch = vec![0.0; MultiBsplineRegular::<f64, 1>::construction_scratch_len(dims)];
        group.bench_with_input(
            BenchmarkId::new("interpn/regular", len),
            &data[..len],
            |b, samples| {
                b.iter(|| {
                    coefficients(dims, samples, &mut coeffs, &mut scratch).unwrap();
                    std::hint::black_box(coeffs.as_slice());
                });
            },
        );

        #[cfg(feature = "__compare_mkl")]
        {
            let mut mkl_coeffs = vec![0.; 4 * (len - 1)];
            let mut mkl_interp =
                mkl_spline::MklInterpF64::new(&data[..len], &mut mkl_coeffs).unwrap();
            group.bench_function(BenchmarkId::new("mkl/WholeSymmetric", len), move |b| {
                b.iter(|| {
                    mkl_interp.calculate_coefficients().unwrap();
                    std::hint::black_box(mkl_interp.coefficients());
                });
            });
        }

        for ext in [WholeSymmetric, Periodic, Reflected] {
            group.bench_with_input(
                BenchmarkId::new(format!("splinterpol/unblocked/{ext:?}"), len),
                &data[..len],
                |b, samples| {
                    b.iter_batched_ref(
                        || samples.to_vec(),
                        |work| {
                            coefficients1d(&BSpline3::new(), work, epsilon, ext, RowBlockSize::<1>)
                                .unwrap();
                            std::hint::black_box(work.as_slice());
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );

            group.bench_with_input(
                BenchmarkId::new(format!("splinterpol/blocked-8/{ext:?}"), len),
                &data[..len],
                |b, samples| {
                    b.iter_batched_ref(
                        || samples.to_vec(),
                        |work| {
                            coefficients1d(&BSpline3::new(), work, epsilon, ext, RowBlockSize::<8>)
                                .unwrap();
                            std::hint::black_box(work.as_slice());
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );
        }
    }
}

fn bench_coefficients_2d(c: &mut Criterion) {
    // Initialize the pool outside timing and respect RAYON_NUM_THREADS.
    let max_threads = rayon::current_num_threads();
    let max_len = *SIZES.iter().max().unwrap();
    let image = make_random_data(max_len * max_len, -100f64..100f64, Some(313373));

    let mut group = c.benchmark_group("Coefficients-2D");
    let epsilon = 1e-12;

    for len in SIZES {
        group.throughput(Throughput::Elements(len as u64 * len as u64));

        // interpn get's to reuse the cache and scratch in the hope that this
        // will produce a benchmark that is most favorable to interpn, since
        // it's not my crate and I don't want to be biased.
        let dims = [len, len];
        let mut coeffs = vec![0.0; MultiBsplineRegular::<f64, 2>::coeff_storage_len(dims)];
        let mut scratch = vec![0.0; MultiBsplineRegular::<f64, 2>::construction_scratch_len(dims)];
        group.bench_with_input(
            BenchmarkId::new("interpn/regular", format!("{len}x{len}")),
            &image[..(len * len)],
            |b, samples| {
                b.iter(|| {
                    coefficients(dims, samples, &mut coeffs, &mut scratch).unwrap();
                    std::hint::black_box(coeffs.as_slice());
                });
            },
        );

        let mut parallel_scratch = vec![
            0.0;
            MultiBsplineRegular::<f64, 2>::parallel_construction_scratch_len(dims, max_threads)
        ];
        group.bench_with_input(
            BenchmarkId::new("interpn/parallel", format!("{len}x{len}")),
            &image[..(len * len)],
            |b, samples| {
                b.iter(|| {
                    coefficients_par(
                        dims,
                        samples,
                        &mut coeffs,
                        &mut parallel_scratch,
                        max_threads,
                    )
                    .unwrap();
                    std::hint::black_box(coeffs.as_slice());
                });
            },
        );

        for ext in [WholeSymmetric, Periodic, Reflected] {
            group.bench_with_input(
                BenchmarkId::new(
                    format!("splinterpol/unblocked/{ext:?}"),
                    format!("{len}x{len}"),
                ),
                &image[..(len * len)],
                |b, samples| {
                    b.iter_batched_ref(
                        || samples.to_vec(),
                        |work| {
                            coefficients2d(
                                &BSpline3::new(),
                                splinterpol::RowMajor {
                                    data: work,
                                    width: len.try_into().unwrap(),
                                    height: len.try_into().unwrap(),
                                },
                                epsilon,
                                ext,
                                RowBlockSize::<1>,
                                ColBlockSize::<1>,
                            )
                            .unwrap();
                            std::hint::black_box(work.as_slice());
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );

            group.bench_with_input(
                BenchmarkId::new(
                    format!("splinterpol/blocked-8x1/{ext:?}"),
                    format!("{len}x{len}"),
                ),
                &image[..(len * len)],
                |b, samples| {
                    b.iter_batched_ref(
                        || samples.to_vec(),
                        |work| {
                            coefficients2d(
                                &BSpline3::new(),
                                splinterpol::RowMajor {
                                    data: work,
                                    width: len.try_into().unwrap(),
                                    height: len.try_into().unwrap(),
                                },
                                epsilon,
                                ext,
                                RowBlockSize::<8>,
                                ColBlockSize::<1>,
                            )
                            .unwrap();
                            std::hint::black_box(work.as_slice());
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );

            group.bench_with_input(
                BenchmarkId::new(
                    format!("splinterpol/blocked-8x8/{ext:?}"),
                    format!("{len}x{len}"),
                ),
                &image[..(len * len)],
                |b, samples| {
                    b.iter_batched_ref(
                        || samples.to_vec(),
                        |work| {
                            coefficients2d(
                                &BSpline3::new(),
                                splinterpol::RowMajor {
                                    data: work,
                                    width: len.try_into().unwrap(),
                                    height: len.try_into().unwrap(),
                                },
                                epsilon,
                                ext,
                                RowBlockSize::<8>,
                                ColBlockSize::<8>,
                            )
                            .unwrap();
                            std::hint::black_box(work.as_slice());
                        },
                        criterion::BatchSize::LargeInput,
                    );
                },
            );
        }
    }
}

criterion_group!(coeffs, bench_coefficients_1d, bench_coefficients_2d);
criterion_main!(coeffs);

//NOTE(geo) I know there's a lot of code duplication in there, but I'd like
// to have one end-to-end benchmark with random resampling in there

// NOTE(geo) #2
// I've decided to do the randomness in this benchmark a little differently
// than in the specialized interpolation and coefficient tests. Here, we
// use a fresh random access pattern in each iteration so that the average
// result includes many random access patterns. This way I hope to have a more
// truthful picture of random sampling at the expense of possibly some additional
// measurement noise. Also I can't enforce anymore that all the different
// backends get the same random sequences in their tests. However, I think
// this should still be fine (especially in the 1D case since there are 100k+
// up to millions of runs). In the 2D case it might be noisier.
//
// Results show that interpn and splinterpol are similarly affected but
// that MKL stays pretty agnostic about this change. If I understand MKL correctly
// that makes sense, since MKL stores a 4-tuple of coefficients per sampling
// site, so it's likely more insensitive to caching effects. Might be wrong
// about this though.

use criterion::{BatchSize, BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use interpn::MultiBsplineRegular;

const SIZES: [usize; 6] = [127usize, 128, 129, 4095, 4096, 4097];

mod shared;
use rand::{
    SeedableRng,
    rngs::{ThreadRng, Xoshiro256PlusPlus},
};
use shared::{make_random_data, make_random_data_with_rng};
use splinterpol::{
    Coord, RowMajor,
    bspline::BSpline3,
    inplace::{SplineInterpolator1d, SplineInterpolator2d},
};

// this tests a random resampling of ALL the points in the 1D dataset
// buffers are reused as much as possible and allocations happens outside the
// times loops
fn bench_end2end_1d_random_resample(c: &mut Criterion) {
    // NOTE: see the coefficients benchmarks for some rationale on reusing the
    // data.
    let max_len = *SIZES.iter().max().unwrap();
    let data = make_random_data(max_len, -100f64..100f64, Some(313373));

    let epsilon = 1e-12;
    let mut group = c.benchmark_group("End2End-1D/random-resample");

    // TODO: consider making this seedable (?)
    let mut rng = Xoshiro256PlusPlus::from_rng(&mut ThreadRng::default());

    for len in SIZES {
        group.throughput(Throughput::Elements(len as u64));

        let dims = [len];
        let mut coeffs = vec![0.0; MultiBsplineRegular::<f64, 1>::coeff_storage_len(dims)];
        let mut scratch = vec![0.0; MultiBsplineRegular::<f64, 1>::construction_scratch_len(dims)];
        let mut out = vec![0.0; len];
        group.bench_function(BenchmarkId::new("interpn/regular", len), |b| {
            b.iter_batched_ref(
                || make_random_data_with_rng(&mut rng, len, 0.0..((len - 1) as f64)),
                |positions| {
                    let interpn = MultiBsplineRegular::from_values_with_workspace(
                        dims,
                        [0.0],
                        [1.0],
                        &data[..len],
                        &mut coeffs,
                        &mut scratch,
                        false,
                    )
                    .unwrap();

                    interpn.interp(&[positions.as_slice()], &mut out).unwrap();
                    std::hint::black_box(out.as_slice());
                },
                BatchSize::PerIteration,
            );
        });

        #[cfg(feature = "__compare_mkl")]
        {
            let mut mkl_coeffs = vec![0.; 4 * (len - 1)];
            let mut mkl_out = vec![0.0; len];
            group.bench_function(BenchmarkId::new("mkl/WholeSymmetric", len), |b| {
                b.iter_batched_ref(
                    || make_random_data_with_rng(&mut rng, len, 0.0..((len - 1) as f64)),
                    |positions| {
                        let mut mkl_interp =
                            mkl_spline::MklInterpF64::new(&data[..len], &mut mkl_coeffs).unwrap();
                        mkl_interp.calculate_coefficients().unwrap();
                        mkl_interp
                            .interp_multi_generic(positions.as_slice(), &mut mkl_out)
                            .unwrap();
                        std::hint::black_box(mkl_out.as_slice());
                    },
                    BatchSize::PerIteration,
                );
            });
        }

        let mut splinterpol_coeffs = vec![0.0; len];
        group.bench_function(BenchmarkId::new("splinterpol/WholeSymmetric", len), |b| {
            b.iter_batched_ref(
                || make_random_data_with_rng(&mut rng, len, 0.0..((len - 1) as f64)),
                |positions| {
                    splinterpol_coeffs.copy_from_slice(&data[..len]);
                    let interp1d = SplineInterpolator1d::builder(BSpline3::<f64>::new())
                        .extension(splinterpol::inplace::InplaceBoundaryExtension::WholeSymmetric)
                        .row_blocksize::<8>()
                        .build(splinterpol_coeffs.as_mut_slice(), epsilon)
                        .unwrap();

                    interp1d.interp_multi(positions, &mut out).unwrap();
                    std::hint::black_box(out.as_slice());
                },
                BatchSize::PerIteration,
            );
        });
    }
}

fn bench_end2end_2d_random_resample(c: &mut Criterion) {
    let max_len = *SIZES.iter().max().unwrap();
    let data = make_random_data(max_len * max_len, -100f64..100f64, Some(313373));
    let epsilon = 1e-12;
    let mut group = c.benchmark_group("End2End-2D/random-resample");

    let mut rng = Xoshiro256PlusPlus::from_rng(&mut ThreadRng::default());

    for len in SIZES {
        let pixel_count = len * len;
        group.throughput(Throughput::Elements(pixel_count as u64));
        let dims = [len, len];
        let mut coeffs = vec![0.0; MultiBsplineRegular::<f64, 2>::coeff_storage_len(dims)];
        let mut scratch = vec![0.0; MultiBsplineRegular::<f64, 2>::construction_scratch_len(dims)];
        let mut out = vec![0.0; pixel_count];
        group.bench_function(
            BenchmarkId::new("interpn/regular", format!("{len}x{len}")),
            |b| {
                b.iter_batched_ref(
                    || {
                        (
                            make_random_data_with_rng(
                                &mut rng,
                                pixel_count,
                                0.0..((len - 1) as f64),
                            ),
                            make_random_data_with_rng(
                                &mut rng,
                                pixel_count,
                                0.0..((len - 1) as f64),
                            ),
                        )
                    },
                    |positions| {
                        let interpn = MultiBsplineRegular::from_values_with_workspace(
                            dims,
                            [0.0, 0.0],
                            [1.0, 1.0],
                            &data[..pixel_count],
                            &mut coeffs,
                            &mut scratch,
                            false,
                        )
                        .unwrap();

                        interpn
                            .interp(&[positions.0.as_slice(), positions.1.as_slice()], &mut out)
                            .unwrap();
                        std::hint::black_box(out.as_slice());
                    },
                    BatchSize::PerIteration,
                );
            },
        );
        let mut splinterpol_coeffs = vec![0.0; pixel_count];
        group.bench_function(
            BenchmarkId::new("splinterpol/WholeSymmetric", format!("{len}x{len}")),
            |b| {
                b.iter_batched_ref(
                    || {
                        let access_indices = make_random_data_with_rng(
                            &mut rng,
                            2 * pixel_count,
                            0.0..((len - 1) as f64),
                        );
                        let (rows, cols) = access_indices.split_at(pixel_count);
                        rows.iter()
                            .zip(cols)
                            .map(|(&row, &col)| Coord::row(row).col(col))
                            .collect::<Vec<_>>()
                    },
                    |positions| {
                        splinterpol_coeffs.copy_from_slice(&data[..pixel_count]);
                        let interp2d = SplineInterpolator2d::builder(BSpline3::<f64>::new())
                            .extension(
                                splinterpol::inplace::InplaceBoundaryExtension::WholeSymmetric,
                            )
                            .row_blocksize::<8>()
                            .col_blocksize::<8>()
                            .build(
                                RowMajor {
                                    data: splinterpol_coeffs.as_mut_slice(),
                                    width: len.try_into().unwrap(),
                                    height: len.try_into().unwrap(),
                                },
                                epsilon,
                            )
                            .unwrap();

                        for (&position, value) in positions.iter().zip(out.iter_mut()) {
                            *value = interp2d.interp_single(position).unwrap();
                        }
                        std::hint::black_box(out.as_slice());
                    },
                    BatchSize::PerIteration,
                );
            },
        );
    }
}

criterion_group!(
    end2end,
    bench_end2end_1d_random_resample,
    bench_end2end_2d_random_resample
);
criterion_main!(end2end);

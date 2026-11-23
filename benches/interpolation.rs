use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use interpn::MultiBsplineRegular;

const SIZES: [usize; 6] = [127usize, 128, 129, 4095, 4096, 4097];

mod shared;
use shared::make_random_data;
use splinterpol::{
    Coord, RowMajor,
    bspline::BSpline3,
    inplace::{SplineInterpolator1d, SplineInterpolator2d},
};

// this tests a random resampling of ALL the points in the 1D dataset
fn bench_interpolation_1d_random_resample(c: &mut Criterion) {
    // NOTE: see the coefficients benchmarks for some rationale on reusing the
    // data.
    let max_len = *SIZES.iter().max().unwrap();
    let data = make_random_data(max_len, -100f64..100f64, Some(313373));

    let epsilon = 1e-12;
    let mut group = c.benchmark_group("Interpolation-1D/random-resample");

    for len in SIZES {
        let access_indices = make_random_data(len, 0f64..((len - 1) as f64), Some(1234874623497));
        let interp1d = SplineInterpolator1d::builder(BSpline3::<f64>::new())
            .extension(splinterpol::inplace::InplaceBoundaryExtension::WholeSymmetric)
            .row_blocksize::<8>()
            .build(data[0..len].to_vec(), epsilon)
            .unwrap();
        let dims = [len];
        let mut coeffs = vec![0.0; MultiBsplineRegular::<f64, 1>::coeff_storage_len(dims)];
        let mut scratch = vec![0.0; MultiBsplineRegular::<f64, 1>::construction_scratch_len(dims)];
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
        let mut out = vec![0.0; len];

        group.throughput(Throughput::Elements(len as u64));

        group.bench_with_input(
            BenchmarkId::new("interpn/regular", len),
            &access_indices,
            |b, positions| {
                b.iter(|| {
                    interpn.interp(&[positions.as_slice()], &mut out).unwrap();
                    std::hint::black_box(out.as_slice());
                });
            },
        );

        #[cfg(feature = "__compare_mkl")]
        {
            let mut mkl_coeffs = vec![0.; 4 * (len - 1)];
            let mut mkl_interp =
                mkl_spline::MklInterpF64::new(&data[..len], &mut mkl_coeffs).unwrap();
            mkl_interp.calculate_coefficients().unwrap();
            let mut mkl_out = vec![0.0; len];
            group.bench_with_input(
                BenchmarkId::new("mkl/WholeSymmetric", len),
                &access_indices,
                |b, positions| {
                    b.iter(|| {
                        mkl_interp
                            .interp_multi_generic(positions.as_slice(), &mut mkl_out)
                            .unwrap();
                        std::hint::black_box(mkl_out.as_slice());
                    });
                },
            );
        }

        group.bench_with_input(
            BenchmarkId::new("splinterpol/WholeSymmetric", len),
            &access_indices,
            |b, positions| {
                b.iter(|| {
                    interp1d.interp_multi(positions, &mut out).unwrap();
                    std::hint::black_box(out.as_slice());
                });
            },
        );
    }
}

fn bench_interpolation_2d_random_resample(c: &mut Criterion) {
    let max_len = *SIZES.iter().max().unwrap();
    let data = make_random_data(max_len * max_len, -100f64..100f64, Some(313373));
    let epsilon = 1e-12;
    let mut group = c.benchmark_group("Interpolation-2D/random-resample");

    for len in SIZES {
        let pixel_count = len * len;
        // Draw both axes together so an environment seed cannot make them identical.
        let access_indices = make_random_data(
            2 * pixel_count,
            0f64..((len - 1) as f64),
            Some(1234874623497),
        );
        let (rows, cols) = access_indices.split_at(pixel_count);
        let positions: Vec<_> = rows
            .iter()
            .zip(cols)
            .map(|(&row, &col)| Coord::row(row).col(col))
            .collect();

        let interp2d = SplineInterpolator2d::builder(BSpline3::<f64>::new())
            .extension(splinterpol::inplace::InplaceBoundaryExtension::WholeSymmetric)
            .row_blocksize::<8>()
            .col_blocksize::<8>()
            .build(
                RowMajor {
                    data: data[..pixel_count].to_vec(),
                    width: len.try_into().unwrap(),
                    height: len.try_into().unwrap(),
                },
                epsilon,
            )
            .unwrap();
        let dims = [len, len];
        let mut coeffs = vec![0.0; MultiBsplineRegular::<f64, 2>::coeff_storage_len(dims)];
        let mut scratch = vec![0.0; MultiBsplineRegular::<f64, 2>::construction_scratch_len(dims)];
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
        let mut out = vec![0.0; pixel_count];

        group.throughput(Throughput::Elements(pixel_count as u64));
        group.bench_with_input(
            BenchmarkId::new("interpn/regular", format!("{len}x{len}")),
            &[rows, cols],
            |b, positions| {
                b.iter(|| {
                    interpn.interp(positions, &mut out).unwrap();
                    std::hint::black_box(out.as_slice());
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("splinterpol/WholeSymmetric", format!("{len}x{len}")),
            &positions,
            |b, positions| {
                b.iter(|| {
                    for (&position, value) in positions.iter().zip(out.iter_mut()) {
                        *value = interp2d.interp_single(position).unwrap();
                    }
                    std::hint::black_box(out.as_slice());
                });
            },
        );
    }
}

criterion_group!(
    interpolation,
    bench_interpolation_1d_random_resample,
    bench_interpolation_2d_random_resample
);
criterion_main!(interpolation);

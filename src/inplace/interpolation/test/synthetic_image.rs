use crate::inplace::builder::RowMajor;
use num_traits::{AsPrimitive, ConstZero, Float};

pub struct SynthImage {
    nrows: usize,
    ncols: usize,
    // the shapes that this image is composed of
    gaussians: Vec<Gaussian>,
}

impl SynthImage {
    // construct a new synth image with the given width and height
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            nrows: height,
            ncols: width,
            gaussians: vec![Gaussian {
                center: (height as f64 / 2., width as f64 / 3.),
                ampl: 10.,
                sigma: (1., 2.),
            }],
        }
    }

    pub fn as_row_major<T>(&self) -> RowMajor<Vec<T>>
    where
        f64: AsPrimitive<T>,
        T: Copy + 'static + Float + ConstZero,
    {
        let mut data = vec![T::ZERO; self.nrows * self.ncols];

        for c in 0..self.ncols {
            for r in 0..self.nrows {
                //sum all the gaussians
                let mut gsum = 0.;
                for gaus in self.gaussians.iter() {
                    gsum += gaus.at(c as _, r as _);
                }

                data[c * self.nrows + r] = gsum.as_();
            }
        }
        todo!()
    }
}

struct Gaussian {
    center: (f64, f64),
    ampl: f64,
    sigma: (f64, f64),
}

impl Gaussian {
    fn at(&self, x: f64, y: f64) -> f64 {
        self.ampl
            * (((x - self.center.0) / (2. * self.sigma.0)).powi(2)
                + ((y - self.center.1) / (2. * self.sigma.1)).powi(2))
            .exp()
    }
}

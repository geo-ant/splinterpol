use crate::RowMajor;
use crate::inplace::interpolation::test::synthetic_image::shapes::{
    Blob, ExponentialPeak, Gaussian, Line, Shape, Sinusoid,
};
use core::f64;
use num_traits::{AsPrimitive, ConstZero, Float};

mod shapes;

pub struct SynthImage {
    nrows: usize,
    ncols: usize,
    // forgive me father for I have sinned, but this is test code.
    shapes: Vec<Box<dyn Shape>>,
}

impl SynthImage {
    // construct a new synth image with the given width and height
    pub fn new(width: usize, height: usize) -> Self {
        let fwidth = width as f64;
        let fheight = height as f64;
        Self {
            nrows: height,
            ncols: width,
            shapes: vec![
                Gaussian {
                    center: (height as f64 / 2., width as f64 / 3.),
                    ampl: 10.,
                    sigma: (10., 20.),
                }
                .into_boxed(),
                Gaussian {
                    center: (0., fheight),
                    ampl: 40.,
                    sigma: (fheight / 2., fwidth / 2.),
                }
                .into_boxed(),
                Gaussian {
                    center: (fwidth / 5., fheight / 4.),
                    ampl: 100.,
                    sigma: (5., 6.),
                }
                .into_boxed(),
                Sinusoid::new(10., 30., fwidth / 250.).into_boxed(),
                Sinusoid::new(5., -20., fheight / 100.).into_boxed(),
                Line::new(-10., (fwidth / 2., fheight / 2.), 0., 10., 0.3).into_boxed(),
                ExponentialPeak::new(100., (fwidth / 8., fheight / 7.), fwidth / 100.).into_boxed(),
                Blob::new(-50., (fwidth / 2., fheight), fheight / 50., 0.4).into_boxed(),
            ],
        }
    }

    pub fn as_row_major<T>(&self) -> RowMajor<Vec<T>>
    where
        f64: AsPrimitive<T>,
        T: Copy + 'static + Float + ConstZero,
    {
        let mut data = vec![T::ZERO; self.nrows * self.ncols];

        // NOTE PERF: this is slow as fuck since I'm uselessly iterating most
        // of the image. But it should be fine for testing...
        for r in 0..self.nrows {
            for c in 0..self.ncols {
                data[r * self.ncols + c] = self.at(r as _, c as _).as_();
            }
        }

        RowMajor {
            data,
            width: self.ncols.try_into().unwrap(),
            height: self.nrows.try_into().unwrap(),
        }
    }

    #[inline]
    pub fn at(&self, r: f64, c: f64) -> f64 {
        let mut sum = 0.;
        for shape in self.shapes.iter() {
            sum += shape.at(r as _, c as _);
        }
        sum
    }
}

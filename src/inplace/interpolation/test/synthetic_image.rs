use crate::inplace::builder::RowMajor;
use crate::numerics::FloatVector;
use core::f64;
use num_traits::{AsPrimitive, ConstZero, Float};

pub struct SynthImage {
    nrows: usize,
    ncols: usize,
    // all of the following shapes are additive and compose the image.
    gaussians: Vec<Gaussian>,
    sins: Vec<Sinusoid>,
}

impl SynthImage {
    // construct a new synth image with the given width and height
    pub fn new(width: usize, height: usize) -> Self {
        let fwidth = width as f64;
        let fheight = height as f64;
        Self {
            nrows: height,
            ncols: width,
            gaussians: vec![
                Gaussian {
                    center: (height as f64 / 2., width as f64 / 3.),
                    ampl: 10.,
                    sigma: (1., 2.),
                },
                Gaussian {
                    center: (0., fheight),
                    ampl: 40.,
                    sigma: (fheight / 2., fwidth / 2.),
                },
            ],
            sins: vec![],
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
        for c in 0..self.ncols {
            for r in 0..self.nrows {
                let mut gsum = 0.;
                for gauss in self.gaussians.iter() {
                    gsum += gauss.at(r as _, c as _);
                }

                let mut sinsum = 0.;
                for sin in self.sins.iter() {
                    sinsum += sin.at(r as _, c as _);
                }

                data[c * self.nrows + r] = gsum.as_() + sinsum.as_();
            }
        }
        RowMajor {
            data,
            width: self.ncols,
            height: self.nrows,
        }
    }
}

struct Gaussian {
    // use (row,col) indexing logic
    center: (f64, f64),
    ampl: f64,
    sigma: (f64, f64),
}

impl Gaussian {
    #[inline(always)]
    fn at(&self, r: f64, c: f64) -> f64 {
        self.ampl
            * (-((r - self.center.0) / (2. * self.sigma.0)).powi(2)
                - ((c - self.center.1) / (2. * self.sigma.1)).powi(2))
            .exp()
    }
}

struct Sinusoid {
    // ampl*sin (k*x), where k = (k1,k2), x = (x1,x2)
    k: (f64, f64),
    ampl: f64,
}

impl Sinusoid {
    #[inline(always)]
    fn at(&self, r: f64, c: f64) -> f64 {
        self.ampl * (self.k.0 * r + self.k.1 * c).sin()
    }
}

// this guy is just a circle with a certain value, but it also has
// a small decay at the edges using an 1+cos decay of a certain width
// delta to go do zero. If that width is longer than one pixel, it should still
// be pretty much bandlimited, but it can also be made smaller to give structures
// that violate the bandlimiting.
struct Blob {
    center: (f64, f64),
    ampl: f64,
    radius: f64,
    // in pixels, choose < 1 for non-bandlimited structure in image.
    // The total radius covered by the structure is `radius+rim_width`
    rim_width: f64,
}

impl Blob {
    #[inline(always)]
    fn at(&self, r: f64, c: f64) -> f64 {
        let rd = (r * r + c * c).sqrt();
        if rd <= self.radius {
            self.ampl
        } else if rd < self.radius + self.rim_width {
            let d = self.radius + self.rim_width - rd;
            0.5 + 0.5 * (d * 2. * core::f64::consts::PI / self.rim_width).cos()
        } else {
            0.
        }
    }
}

struct Line {
    point: [f64; 2],
    // normalized to length 1
    normal: [f64; 2],
    // thickness of line = 2*radius
    radius: f64,
    // rim width (as in blob) on each side
    rim_width: f64,
    // brightness value to be added
    ampl: f64,
}

impl Line {
    #[inline(always)]
    fn at(&self, r: f64, c: f64) -> f64 {
        let d = ([r, c].sub(self.point).dot(self.normal)).abs();
        if d <= self.radius {
            self.ampl
        } else if d <= self.radius + self.rim_width {
            let d = self.radius + self.rim_width - d;
            0.5 + 0.5 * (d * 2. * core::f64::consts::PI / self.rim_width).cos()
        } else {
            0.
        }
    }
}

// NOTE(geo): first draft my own, then had llm fix some tedious bugs in code,
// then some manual work again... etc
use crate::numerics::FloatVector;
use core::f64::consts::PI;

// yes, I know I hate it too I'm going full OOP here because I'm too lazy
// to optimize this code really. But I'll only lazily create a few images
// and not use this on the fly, so that shouldn't completely kill my test
// performance I hope.
pub trait Shape: 'static {
    /// eval the shape at the given row and column coord
    fn at(&self, r: f64, c: f64) -> f64;
    fn into_boxed(self) -> Box<dyn Shape + 'static>
    where
        Self: Sized,
    {
        Box::new(self)
    }
}

/// Constant amplitude up to `radius`, followed by a raised-cosine
/// decay to zero over `rim_width`.
#[inline(always)]
fn tapered(ampl: f64, distance: f64, radius: f64, rim_width: f64) -> f64 {
    if distance <= radius {
        ampl
    } else if distance < radius + rim_width {
        let t = (distance - radius) / rim_width;
        ampl * 0.5 * (1.0 + (PI * t).cos())
    } else {
        0.0
    }
}

pub struct Gaussian {
    // use (row,col) indexing logic
    pub center: (f64, f64),
    pub ampl: f64,
    pub sigma: (f64, f64),
}

impl Shape for Gaussian {
    fn at(&self, r: f64, c: f64) -> f64 {
        self.ampl
            * f64::exp(
                -(r - self.center.0).powi(2) / (2. * self.sigma.0.powi(2))
                    - (c - self.center.1).powi(2) / (2. * self.sigma.1.powi(2)),
            )
    }
}

pub struct Sinusoid {
    // ampl*sin(2π k·x), where x = (row,col)
    k: (f64, f64),
    ampl: f64,
}

impl Sinusoid {
    /// Angle of the wavevector relative to the x-axis (columns),
    /// amplitude, and period (in pixels).
    pub fn new(ampl: f64, angle_deg: f64, period: f64) -> Self {
        let (sin, cos) = angle_deg.to_radians().sin_cos();

        Sinusoid {
            k: (sin / period, cos / period),
            ampl,
        }
    }
}

impl Shape for Sinusoid {
    fn at(&self, r: f64, c: f64) -> f64 {
        self.ampl * ((self.k.0 * r + self.k.1 * c) * 2. * PI).sin()
    }
}

pub struct Blob {
    center: (f64, f64),
    ampl: f64,
    radius: f64,
    // The total radius covered is radius + rim_width.
    rim_width: f64,
}

impl Blob {
    pub fn new(ampl: f64, center: (f64, f64), radius: f64, rim_width: f64) -> Self {
        Self {
            center,
            ampl,
            radius,
            rim_width,
        }
    }
}

impl Shape for Blob {
    fn at(&self, r: f64, c: f64) -> f64 {
        let dr = r - self.center.0;
        let dc = c - self.center.1;
        let distance = dr.hypot(dc);

        tapered(self.ampl, distance, self.radius, self.rim_width)
    }
}

pub struct Line {
    aufpunkt: [f64; 2],
    // Unit normal in (row,col) coordinates.
    normal: [f64; 2],
    // Thickness of line = 2*radius.
    radius: f64,
    // Rim width on each side.
    rim_width: f64,
    // Brightness value to be added.
    ampl: f64,
}

impl Line {
    pub fn new(
        value: f64,
        aufpunkt: (f64, f64),
        angle_deg: f64,
        diameter: f64,
        rim_width: f64,
    ) -> Self {
        let (sin, cos) = angle_deg.to_radians().sin_cos();

        // Angle describes the direction of the line, not its normal.
        // Tangent in (row,col): (sin, cos)
        // Normal in (row,col):  (cos,-sin)
        let normal = [cos, -sin];

        Line {
            aufpunkt: [aufpunkt.0, aufpunkt.1],
            normal,
            radius: diameter / 2.,
            rim_width,
            ampl: value,
        }
    }
}

impl Shape for Line {
    fn at(&self, r: f64, c: f64) -> f64 {
        let distance = ([r, c].sub(self.aufpunkt).dot(self.normal)).abs();
        tapered(self.ampl, distance, self.radius, self.rim_width)
    }
}

pub struct ExponentialPeak {
    // row, col
    center: (f64, f64),
    // decay constant
    tau: f64,
    // amplitude at max
    ampl: f64,
}

impl ExponentialPeak {
    pub fn new(ampl: f64, center: (f64, f64), tau: f64) -> Self {
        Self {
            center,
            tau: tau.abs(),
            ampl,
        }
    }
}

impl Shape for ExponentialPeak {
    fn at(&self, r: f64, c: f64) -> f64 {
        let r = ((r - self.center.0).powi(2) + (c - self.center.1).powi(2)).sqrt();
        self.ampl * f64::exp(-r / self.tau)
    }
}

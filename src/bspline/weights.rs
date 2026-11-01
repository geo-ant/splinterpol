use crate::{bspline::*, numerics::FusedMuladd};
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float};

#[cfg(test)]
mod test;

// NOTE: why the stupid constant for the number of weights,
// even though the number of weights is defined by each degree only?
// See the comments on the other trait for a rant.
/// Every B-Spline in this crate implements this tait. It's used during the
/// interpolation step.
///
/// **Note** The number `NWEIGHTS` is the number of elements in the weight
/// vector (`interpn` calls it the _footprint_) used for evaluating the values
/// of the spline around the grid points for interpolation. The number of
/// weights is exactly `n+1` for a BSpline with polynomial degree `n`.
pub trait BSplineWeights<const NWEIGHTS: usize>: BSpline {
    /// evaluate the centered B-Spline at the given coordinate.
    fn eval(self, coord: Self::Float) -> Self::Float;
    /// The radius of the B-Spline support. This is always
    /// `(NWEIGHTS as Float)/2.0`
    fn radius() -> Self::Float;
    //TODO!
    fn weights(t: Self::Float) -> [Self::Float; NWEIGHTS];
}

impl<F> BSplineWeights<3> for BSpline2<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn eval(self, x: Self::Float) -> Self::Float {
        let x = x.abs();
        if x < 0.5.as_() {
            let y = x;
            return 1.5.as_() - 2.0.as_() * (y).powi(2);
        }
        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            return (y).powi(2);
        }
        F::ZERO
    }

    #[inline(always)]
    fn radius() -> Self::Float {
        1.5.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 3] {
        debug_assert!(t >= 0.5.as_() && t <= 1.5.as_());

        let q = t - 0.5.as_();
        let s = F::ONE - q;
        let z = t - F::ONE;

        [s * s, (z * z).fmuladd(-2.0.as_(), 1.5.as_()), q * q]
    }
}

impl<F> BSplineWeights<4> for BSpline3<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn eval(self, x: F) -> F {
        // generated this function with `python scripts/bspline-eval-function.py 3`
        // and then did optimizations and cleanups by hand
        let x = x.abs();
        if x < F::ONE {
            let y = x;
            return (y)
                .powi(2)
                .fmuladd(y.fmuladd(3.0.as_(), -6.0.as_()), 4.0.as_());
        }
        if x < 2.0.as_() {
            let y = (2.0.as_()) - x;
            return (y).powi(3);
        }
        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        Self::Float::ONE + Self::Float::ONE
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 4] {
        debug_assert!(t >= F::ONE && t <= 2.0.as_());

        let q = t - F::ONE;
        let s = F::ONE - q;
        let q2 = q * q;
        let s2 = s * s;
        let three_q = q * 3.0.as_();

        [
            s2 * s,
            q2.fmuladd(three_q - 6.0.as_(), 4.0.as_()),
            s2.fmuladd(-three_q - 3.0.as_(), 4.0.as_()),
            q2 * q,
        ]
    }
}

impl<F> BSplineWeights<5> for BSpline4<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: F) -> F {
        // NOTE(geo): also created with the script and then had an llm apply the muladd.
        // The muladd refactoring is very tedious but schematic because we
        // can always go from inside to outside. The one thing we should consider
        // is whether it makes sense to split the calculation into multiple
        // dependency chains for larger calculations, which I don't think
        // is neccessarily possible if we use the Horner rewrite for calculating
        // the polynomials, which is what Briand & Monasse also do... gotta think
        // about that.
        let x = x.abs();

        if x < 0.5.as_() {
            let x2 = x * x;
            return x2.fmuladd(x2.fmuladd(6.0.as_(), (-15.0).as_()), 14.375.as_());
        }

        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            let p = y.fmuladd((-4.0).as_(), 4.0.as_());
            let p = y.fmuladd(p, 6.0.as_());
            let p = y.fmuladd(p, 4.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 2.5.as_() {
            let y = 2.5.as_() - x;
            let y2 = y * y;
            return y2 * y2;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        2.5.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 5] {
        todo!()
    }
}

impl<F> BSplineWeights<6> for BSpline5<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // from the script, then manually refined, then let an llm handle
        // the tedious work of writing this as muladds
        let x = x.abs();

        if x < F::ONE {
            let x2 = x * x;

            // 66 + x² * (-60 + x² * (30 - 10*x))
            let p = x.fmuladd((-10.0).as_(), 30.0.as_());
            let p = x2.fmuladd(p, (-60.0).as_());
            return x2.fmuladd(p, 66.0.as_());
        }

        if x < 2.0.as_() {
            let y = 2.0.as_() - x;

            // 1 + y * (5 + y * (10 + y * (10 + y * (5 - 5*y))))
            let p = y.fmuladd((-5.0).as_(), 5.0.as_());
            let p = y.fmuladd(p, 10.0.as_());
            let p = y.fmuladd(p, 10.0.as_());
            let p = y.fmuladd(p, 5.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 3.0.as_() {
            let y = 3.0.as_() - x;
            let y2 = y * y;
            let y4 = y2 * y2;
            return y4 * y;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        3.0.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 6] {
        todo!()
    }
}

impl<F> BSplineWeights<7> for BSpline6<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE(geo): same comments as above applies..
        let x = x.abs();

        if x < 0.5.as_() {
            let y2 = x * x;
            let p = y2.fmuladd((-20.0).as_(), 105.0.as_());
            let p = y2.fmuladd(p, (-288.75).as_());
            return y2.fmuladd(p, 367.9375.as_());
        }

        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            let p = y.fmuladd(15.0.as_(), (-30.0).as_());
            let p = y.fmuladd(p, (-45.0).as_());
            let p = y.fmuladd(p, 20.0.as_());
            let p = y.fmuladd(p, 135.0.as_());
            let p = y.fmuladd(p, 150.0.as_());
            return y.fmuladd(p, 57.0.as_());
        }

        if x < 2.5.as_() {
            let y = 2.5.as_() - x;
            let p = y.fmuladd((-6.0).as_(), 6.0.as_());
            let p = y.fmuladd(p, 15.0.as_());
            let p = y.fmuladd(p, 20.0.as_());
            let p = y.fmuladd(p, 15.0.as_());
            let p = y.fmuladd(p, 6.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 3.5.as_() {
            let y = 3.5.as_() - x;
            let y2 = y * y;
            let y3 = y2 * y;
            return y3 * y3;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        3.5.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 7] {
        todo!()
    }
}

impl<F> BSplineWeights<8> for BSpline7<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE(geo): I directly had the llm translate the output of the python
        // script
        let x = x.abs();

        if x < F::ONE {
            let y2 = x * x;
            let p = x.fmuladd(35.0.as_(), (-140.0).as_());
            let p = y2.fmuladd(p, 560.0.as_());
            let p = y2.fmuladd(p, (-1680.0).as_());
            return y2.fmuladd(p, 2416.0.as_());
        }

        if x < 2.0.as_() {
            let y = 2.0.as_() - x;
            let y2 = y * y;
            let p = y.fmuladd(21.0.as_(), (-42.0).as_());
            let p = y.fmuladd(p, (-84.0).as_());
            let p = y2.fmuladd(p, 280.0.as_());
            let p = y.fmuladd(p, 504.0.as_());
            let p = y.fmuladd(p, 392.0.as_());
            return y.fmuladd(p, 120.0.as_());
        }

        if x < 3.0.as_() {
            let y = 3.0.as_() - x;
            let p = y.fmuladd((-7.0).as_(), 7.0.as_());
            let p = y.fmuladd(p, 21.0.as_());
            let p = y.fmuladd(p, 35.0.as_());
            let p = y.fmuladd(p, 35.0.as_());
            let p = y.fmuladd(p, 21.0.as_());
            let p = y.fmuladd(p, 7.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 4.0.as_() {
            let y = 4.0.as_() - x;
            let y2 = y * y;
            let y3 = y2 * y;
            let y4 = y2 * y2;
            return y3 * y4;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        4.0.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 8] {
        todo!()
    }
}

impl<F> BSplineWeights<9> for BSpline8<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE(geo): same comment as above
        let x = x.abs();

        if x < 0.5.as_() {
            let y2 = x * x;
            let p = y2.fmuladd(70.0.as_(), (-630.0).as_());
            let p = y2.fmuladd(p, 3386.25.as_());
            let p = y2.fmuladd(p, (-11379.375).as_());
            return y2.fmuladd(p, 18261.7734375.as_());
        }

        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            let p = y.fmuladd((-56.0).as_(), 168.0.as_());
            let p = y.fmuladd(p, 252.0.as_());
            let p = y.fmuladd(p, (-504.0).as_());
            let p = y.fmuladd(p, (-1890.0).as_());
            let p = y.fmuladd(p, (-504.0).as_());
            let p = y.fmuladd(p, 5292.0.as_());
            let p = y.fmuladd(p, 8568.0.as_());
            return y.fmuladd(p, 4293.0.as_());
        }

        if x < 2.5.as_() {
            let y = 2.5.as_() - x;
            let p = y.fmuladd(28.0.as_(), (-56.0).as_());
            let p = y.fmuladd(p, (-140.0).as_());
            let p = y.fmuladd(p, (-56.0).as_());
            let p = y.fmuladd(p, 490.0.as_());
            let p = y.fmuladd(p, 1288.0.as_());
            let p = y.fmuladd(p, 1540.0.as_());
            let p = y.fmuladd(p, 952.0.as_());
            return y.fmuladd(p, 247.0.as_());
        }

        if x < 3.5.as_() {
            let y = 3.5.as_() - x;
            let p = y.fmuladd((-8.0).as_(), 8.0.as_());
            let p = y.fmuladd(p, 28.0.as_());
            let p = y.fmuladd(p, 56.0.as_());
            let p = y.fmuladd(p, 70.0.as_());
            let p = y.fmuladd(p, 56.0.as_());
            let p = y.fmuladd(p, 28.0.as_());
            let p = y.fmuladd(p, 8.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 4.5.as_() {
            let y = 4.5.as_() - x;
            let y2 = y * y;
            let y4 = y2 * y2;
            return y4 * y4;
        }

        F::ZERO
    }
    #[inline]
    fn radius() -> Self::Float {
        4.5.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 9] {
        todo!()
    }
}

impl<F> BSplineWeights<10> for BSpline9<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE: same comment as before
        let x = x.abs();

        if x < F::ONE {
            let y2 = x * x;
            let p = x.fmuladd((-126.0).as_(), 630.0.as_());
            let p = y2.fmuladd(p, (-4200.0).as_());
            let p = y2.fmuladd(p, 23940.0.as_());
            let p = y2.fmuladd(p, (-88200.0).as_());
            return y2.fmuladd(p, 156190.0.as_());
        }

        if x < 2.0.as_() {
            let y = 2.0.as_() - x;
            let p = y.fmuladd((-84.0).as_(), 252.0.as_());
            let p = y.fmuladd(p, 504.0.as_());
            let p = y.fmuladd(p, (-672.0).as_());
            let p = y.fmuladd(p, (-4284.0).as_());
            let p = y.fmuladd(p, (-4032.0).as_());
            let p = y.fmuladd(p, 11256.0.as_());
            let p = y.fmuladd(p, 34272.0.as_());
            let p = y.fmuladd(p, 36414.0.as_());
            return y.fmuladd(p, 14608.0.as_());
        }

        if x < 3.0.as_() {
            let y = 3.0.as_() - x;
            let p = y.fmuladd(36.0.as_(), (-72.0).as_());
            let p = y.fmuladd(p, (-216.0).as_());
            let p = y.fmuladd(p, (-168.0).as_());
            let p = y.fmuladd(p, 756.0.as_());
            let p = y.fmuladd(p, 2772.0.as_());
            let p = y.fmuladd(p, 4536.0.as_());
            let p = y.fmuladd(p, 4248.0.as_());
            let p = y.fmuladd(p, 2214.0.as_());
            return y.fmuladd(p, 502.0.as_());
        }

        if x < 4.0.as_() {
            let y = 4.0.as_() - x;
            let p = y.fmuladd((-9.0).as_(), 9.0.as_());
            let p = y.fmuladd(p, 36.0.as_());
            let p = y.fmuladd(p, 84.0.as_());
            let p = y.fmuladd(p, 126.0.as_());
            let p = y.fmuladd(p, 126.0.as_());
            let p = y.fmuladd(p, 84.0.as_());
            let p = y.fmuladd(p, 36.0.as_());
            let p = y.fmuladd(p, 9.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 5.0.as_() {
            let y = 5.0.as_() - x;
            let y2 = y * y;
            let y4 = y2 * y2;
            let y8 = y4 * y4;
            return y8 * y;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        5.0.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 10] {
        todo!()
    }
}

impl<F> BSplineWeights<11> for BSpline10<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE(geo): same comment as before
        let x = x.abs();

        if x < 0.5.as_() {
            let y2 = x * x;
            let p = y2.fmuladd((-252.0).as_(), 3465.0.as_());
            let p = y2.fmuladd(p, (-30607.5).as_());
            let p = y2.fmuladd(p, 191585.625.as_());
            let p = y2.fmuladd(p, (-769825.546875).as_());
            return y2.fmuladd(p, 1491301.23828125.as_());
        }

        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            let p = y.fmuladd(210.0.as_(), (-840.0).as_());
            let p = y.fmuladd(p, (-1260.0).as_());
            let p = y.fmuladd(p, 5040.0.as_());
            let p = y.fmuladd(p, 16800.0.as_());
            let p = y.fmuladd(p, (-13608.0).as_());
            let p = y.fmuladd(p, (-119280.0).as_());
            let p = y.fmuladd(p, (-95760.0).as_());
            let p = y.fmuladd(p, 327600.0.as_());
            let p = y.fmuladd(p, 736260.0.as_());
            return y.fmuladd(p, 455192.0.as_());
        }

        if x < 2.5.as_() {
            let y = 2.5.as_() - x;
            let p = y.fmuladd((-120.0).as_(), 360.0.as_());
            let p = y.fmuladd(p, 900.0.as_());
            let p = y.fmuladd(p, (-720.0).as_());
            let p = y.fmuladd(p, (-8400.0).as_());
            let p = y.fmuladd(p, (-13608.0).as_());
            let p = y.fmuladd(p, 16800.0.as_());
            let p = y.fmuladd(p, 100080.0.as_());
            let p = y.fmuladd(p, 171000.0.as_());
            let p = y.fmuladd(p, 141060.0.as_());
            return y.fmuladd(p, 47840.0.as_());
        }

        if x < 3.5.as_() {
            let y = 3.5.as_() - x;
            let p = y.fmuladd(45.0.as_(), (-90.0).as_());
            let p = y.fmuladd(p, (-315.0).as_());
            let p = y.fmuladd(p, (-360.0).as_());
            let p = y.fmuladd(p, 1050.0.as_());
            let p = y.fmuladd(p, 5292.0.as_());
            let p = y.fmuladd(p, 11130.0.as_());
            let p = y.fmuladd(p, 14040.0.as_());
            let p = y.fmuladd(p, 11025.0.as_());
            let p = y.fmuladd(p, 5010.0.as_());
            return y.fmuladd(p, 1013.0.as_());
        }

        if x < 4.5.as_() {
            let y = 4.5.as_() - x;
            let p = y.fmuladd((-10.0).as_(), 10.0.as_());
            let p = y.fmuladd(p, 45.0.as_());
            let p = y.fmuladd(p, 120.0.as_());
            let p = y.fmuladd(p, 210.0.as_());
            let p = y.fmuladd(p, 252.0.as_());
            let p = y.fmuladd(p, 210.0.as_());
            let p = y.fmuladd(p, 120.0.as_());
            let p = y.fmuladd(p, 45.0.as_());
            let p = y.fmuladd(p, 10.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 5.5.as_() {
            let y = 5.5.as_() - x;
            let y2 = y * y;
            let y4 = y2 * y2;
            let y5 = y4 * y;
            return y5 * y5;
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        5.5.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 11] {
        todo!()
    }
}

impl<F> BSplineWeights<12> for BSpline11<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE(geo): same comment as before
        let x = x.abs();

        if x < F::ONE {
            let y2 = x * x;
            let p = x.fmuladd(462.0.as_(), (-2772.0).as_());
            let p = y2.fmuladd(p, 27720.0.as_());
            let p = y2.fmuladd(p, (-255024.0).as_());
            let p = y2.fmuladd(p, 1718640.0.as_());
            let p = y2.fmuladd(p, (-7475160.0).as_());
            return y2.fmuladd(p, 15724248.0.as_());
        }

        if x < 2.0.as_() {
            let y = 2.0.as_() - x;
            let y2 = y * y;
            let p = y.fmuladd(330.0.as_(), (-1320.0).as_());
            let p = y.fmuladd(p, (-2640.0).as_());
            let p = y.fmuladd(p, 7920.0.as_());
            let p = y.fmuladd(p, 39600.0.as_());
            let p = y2.fmuladd(p, (-299376.0).as_());
            let p = y.fmuladd(p, (-538560.0).as_());
            let p = y.fmuladd(p, 574200.0.as_());
            let p = y.fmuladd(p, 3273600.0.as_());
            let p = y.fmuladd(p, 4480872.0.as_());
            return y.fmuladd(p, 2203488.0.as_());
        }

        if x < 3.0.as_() {
            let y = 3.0.as_() - x;
            let p = y.fmuladd((-165.0).as_(), 495.0.as_());
            let p = y.fmuladd(p, 1485.0.as_());
            let p = y.fmuladd(p, (-495.0).as_());
            let p = y.fmuladd(p, (-14850.0).as_());
            let p = y.fmuladd(p, (-34650.0).as_());
            let p = y.fmuladd(p, 12474.0.as_());
            let p = y.fmuladd(p, 236610.0.as_());
            let p = y.fmuladd(p, 586575.0.as_());
            let p = y.fmuladd(p, 748275.0.as_());
            let p = y.fmuladd(p, 515097.0.as_());
            return y.fmuladd(p, 152637.0.as_());
        }

        if x < 4.0.as_() {
            let y = 4.0.as_() - x;
            let p = y.fmuladd(55.0.as_(), (-110.0).as_());
            let p = y.fmuladd(p, (-440.0).as_());
            let p = y.fmuladd(p, (-660.0).as_());
            let p = y.fmuladd(p, 1320.0.as_());
            let p = y.fmuladd(p, 9240.0.as_());
            let p = y.fmuladd(p, 24024.0.as_());
            let p = y.fmuladd(p, 38280.0.as_());
            let p = y.fmuladd(p, 40260.0.as_());
            let p = y.fmuladd(p, 27500.0.as_());
            let p = y.fmuladd(p, 11132.0.as_());
            return y.fmuladd(p, 2036.0.as_());
        }

        if x < 5.0.as_() {
            let y = 5.0.as_() - x;
            let p = y.fmuladd((-11.0).as_(), 11.0.as_());
            let p = y.fmuladd(p, 55.0.as_());
            let p = y.fmuladd(p, 165.0.as_());
            let p = y.fmuladd(p, 330.0.as_());
            let p = y.fmuladd(p, 462.0.as_());
            let p = y.fmuladd(p, 462.0.as_());
            let p = y.fmuladd(p, 330.0.as_());
            let p = y.fmuladd(p, 165.0.as_());
            let p = y.fmuladd(p, 55.0.as_());
            let p = y.fmuladd(p, 11.0.as_());
            return y.fmuladd(p, F::ONE);
        }

        if x < 6.0.as_() {
            let y = 6.0.as_() - x;
            let y2 = y * y;
            let y4 = y2 * y2;
            let y5 = y4 * y;
            let y10 = y5 * y5;
            return y10 * y;
        }

        F::ZERO
    }
    #[inline]
    fn radius() -> Self::Float {
        6.0.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 12] {
        todo!()
    }
}

impl<F> BSplineWeights<13> for BSpline12<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline]
    fn eval(self, x: Self::Float) -> Self::Float {
        // NOTE: generated with python script, llm assisted reformulation
        let x = x.abs();

        if x < 0.5.as_() {
            let y2 = x * x;
            return y2
                .fmuladd(924.0.as_(), -18018.0.as_())
                .fmuladd(y2, 236486.25.as_())
                .fmuladd(y2, -2358105.75.as_())
                .fmuladd(y2, 17075574.140625.as_())
                .fmuladd(y2, -79943684.1328125.as_())
                .fmuladd(y2, 181466922.6474609375.as_());
        }

        if x < 1.5.as_() {
            let y = 1.5.as_() - x;
            return y
                .fmuladd(-792.0.as_(), 3960.0.as_())
                .fmuladd(y, 5940.0.as_())
                .fmuladd(y, -38280.0.as_())
                .fmuladd(y, -121770.0.as_())
                .fmuladd(y, 261360.0.as_())
                .fmuladd(y, 1524600.0.as_())
                .fmuladd(y, -137808.0.as_())
                .fmuladd(y, -10932570.0.as_())
                .fmuladd(y, -14341800.0.as_())
                .fmuladd(y, 29555460.0.as_())
                .fmuladd(y, 90415512.0.as_())
                .fmuladd(y, 66318474.0.as_());
        }

        if x < 2.5.as_() {
            let y = 2.5.as_() - x;
            return y
                .fmuladd(495.0.as_(), -1980.0.as_())
                .fmuladd(y, -4950.0.as_())
                .fmuladd(y, 11220.0.as_())
                .fmuladd(y, 81675.0.as_())
                .fmuladd(y, 59400.0.as_())
                .fmuladd(y, -623700.0.as_())
                .fmuladd(y, -1860408.0.as_())
                .fmuladd(y, -37125.0.as_())
                .fmuladd(y, 10101300.0.as_())
                .fmuladd(y, 23794650.0.as_())
                .fmuladd(y, 24610212.0.as_())
                .fmuladd(y, 10187685.0.as_());
        }

        if x < 3.5.as_() {
            let y = 3.5.as_() - x;
            return y
                .fmuladd(-220.0.as_(), 660.0.as_())
                .fmuladd(y, 2310.0.as_())
                .fmuladd(y, 220.0.as_())
                .fmuladd(y, -24255.0.as_())
                .fmuladd(y, -75240.0.as_())
                .fmuladd(y, -23100.0.as_())
                .fmuladd(y, 475992.0.as_())
                .fmuladd(y, 1638945.0.as_())
                .fmuladd(y, 2883100.0.as_())
                .fmuladd(y, 3023790.0.as_())
                .fmuladd(y, 1807212.0.as_())
                .fmuladd(y, 478271.0.as_());
        }

        if x < 4.5.as_() {
            let y = 4.5.as_() - x;
            return y
                .fmuladd(66.0.as_(), -132.0.as_())
                .fmuladd(y, -594.0.as_())
                .fmuladd(y, -1100.0.as_())
                .fmuladd(y, 1485.0.as_())
                .fmuladd(y, 15048.0.as_())
                .fmuladd(y, 47124.0.as_())
                .fmuladd(y, 91080.0.as_())
                .fmuladd(y, 120285.0.as_())
                .fmuladd(y, 109780.0.as_())
                .fmuladd(y, 66726.0.as_())
                .fmuladd(y, 24420.0.as_())
                .fmuladd(y, 4083.0.as_());
        }

        if x < 5.5.as_() {
            let y = 5.5.as_() - x;
            return y
                .fmuladd(-12.0.as_(), 12.0.as_())
                .fmuladd(y, 66.0.as_())
                .fmuladd(y, 220.0.as_())
                .fmuladd(y, 495.0.as_())
                .fmuladd(y, 792.0.as_())
                .fmuladd(y, 924.0.as_())
                .fmuladd(y, 792.0.as_())
                .fmuladd(y, 495.0.as_())
                .fmuladd(y, 220.0.as_())
                .fmuladd(y, 66.0.as_())
                .fmuladd(y, 12.0.as_())
                .fmuladd(y, F::ONE);
        }

        if x < 6.5.as_() {
            let y = 6.5.as_() - x;
            return y.powi(12);
        }

        F::ZERO
    }

    #[inline]
    fn radius() -> Self::Float {
        6.5.as_()
    }

    #[inline(always)]
    fn weights(t: Self::Float) -> [Self::Float; 13] {
        todo!()
    }
}

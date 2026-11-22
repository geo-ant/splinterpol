use crate::{bspline::*, numerics::FusedMuladd};
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float};

#[cfg(test)]
mod test;

// NOTE: why the stupid constant for the number of weights,
// even though the number of weights is defined by each degree only?
// See the comments on the other trait for a rant.
/// **Internal Trait**. Every B-Spline in this crate implements this tait.
/// It's used during the interpolation step.
///
/// **Note** The number `NWEIGHTS` is the number of elements in the weight
/// vector (`interpn` calls it the _footprint_) used for evaluating the values
/// of the spline around the grid points for interpolation. The number of
/// weights is exactly `n+1` for a B-Spline with polynomial degree `n`.
pub trait BSplineWeights<const NWEIGHTS: usize>: BSpline {
    /// **Internal.** The radius of the B-Spline support. This is always
    /// `(NWEIGHTS as Float)/2.0`
    fn radius() -> Self::Float;
    /// **Internal.** Calculate the B-Spline weights needed for interpolation
    // NOTE: The weights and the interpolation logic need to go hand in hand.
    // This crate expects t to be calculated as `t = x - (x-R).ceil()`
    // where `x` is the coordinate (row or col index) and `R` is the radius
    // of the B-Spline. Then this function will give the weights at the
    // corresponding indices [i0,i0+1, ... i0+NWEIGHTS-1], where i0 is the
    // left hand size index of the footprint/weights `i0 = (x-R).ceil()`.
    fn weights(t: Self::Float) -> [Self::Float; NWEIGHTS];
}

impl<F> BSplineWeights<3> for BSpline2<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        1.5.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 3] {
        debug_assert!(t >= 0.5.as_() && t <= 1.5.as_());

        let q = t - 0.5.as_();
        let s = 1.5.as_() - t;
        let z = t - F::ONE;

        [s * s, z.fmuladd(-(z + z), 1.5.as_()), q * q]
    }
}

impl<F> BSplineWeights<4> for BSpline3<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        Self::Float::ONE + Self::Float::ONE
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 4] {
        debug_assert!(t >= F::ONE && t <= 2.0.as_());

        let q = t - F::ONE;
        let s = 2.0.as_() - t;
        let q2 = q * q;
        let s2 = s * s;
        let three_t = 3.0.as_() * t;

        [
            s2 * s,
            q2.fmuladd(three_t - 9.0.as_(), 4.0.as_()),
            s2.fmuladd(-three_t, 4.0.as_()),
            q2 * q,
        ]
    }
}

impl<F> BSplineWeights<5> for BSpline4<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        2.5.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 5] {
        debug_assert!(t >= 1.5.as_() && t <= 2.5.as_());

        let q = t - 1.5.as_();
        let s = 2.5.as_() - t;
        let z = t - 2.0.as_();
        let q2 = q * q;
        let s2 = s * s;
        let z2 = z * z;
        let h = (q * 4.0.as_()).fmuladd(s, 6.0.as_());

        [
            s2 * s2,
            h.fmuladd(s, 4.0.as_()).fmuladd(s, F::ONE),
            z2.fmuladd(6.0.as_(), -15.0.as_()).fmuladd(z2, 14.375.as_()),
            h.fmuladd(q, 4.0.as_()).fmuladd(q, F::ONE),
            q2 * q2,
        ]
    }
}

impl<F> BSplineWeights<6> for BSpline5<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        3.0.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 6] {
        debug_assert!(t >= 2.0.as_() && t <= 3.0.as_());

        let q = t - 2.0.as_();
        let s = 3.0.as_() - t;
        let q2 = q * q;
        let s2 = s * s;
        let ten_t = t * 10.0.as_();
        let h = (q * 5.0.as_()).fmuladd(s, 10.0.as_());

        [
            s2 * s2 * s,
            h.fmuladd(s, 10.0.as_())
                .fmuladd(s, 5.0.as_())
                .fmuladd(s, F::ONE),
            q2.fmuladd(50.0.as_() - ten_t, -60.0.as_())
                .fmuladd(q2, 66.0.as_()),
            s2.fmuladd(ten_t, -60.0.as_()).fmuladd(s2, 66.0.as_()),
            h.fmuladd(q, 10.0.as_())
                .fmuladd(q, 5.0.as_())
                .fmuladd(q, F::ONE),
            q2 * q2 * q,
        ]
    }
}

impl<F> BSplineWeights<7> for BSpline6<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        3.5.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 7] {
        debug_assert!(t >= 2.5.as_() && t <= 3.5.as_());

        let q = t - 2.5.as_();
        let s = 3.5.as_() - t;
        let z = t - 3.0.as_();
        let z2 = z * z;
        let q3 = q * q * q;
        let s3 = s * s * s;

        let qs = q * s;
        let h = qs.fmuladd(6.0.as_(), 15.0.as_());
        let k = qs.fmuladd(-15.0.as_(), -45.0.as_());

        [
            s3 * s3,
            h.fmuladd(s, 20.0.as_())
                .fmuladd(s, 15.0.as_())
                .fmuladd(s, 6.0.as_())
                .fmuladd(s, F::ONE),
            s.fmuladd(-15.0.as_(), k)
                .fmuladd(s, 20.0.as_())
                .fmuladd(s, 135.0.as_())
                .fmuladd(s, 150.0.as_())
                .fmuladd(s, 57.0.as_()),
            z2.fmuladd(-20.0.as_(), 105.0.as_())
                .fmuladd(z2, -288.75.as_())
                .fmuladd(z2, 367.9375.as_()),
            q.fmuladd(-15.0.as_(), k)
                .fmuladd(q, 20.0.as_())
                .fmuladd(q, 135.0.as_())
                .fmuladd(q, 150.0.as_())
                .fmuladd(q, 57.0.as_()),
            h.fmuladd(q, 20.0.as_())
                .fmuladd(q, 15.0.as_())
                .fmuladd(q, 6.0.as_())
                .fmuladd(q, F::ONE),
            q3 * q3,
        ]
    }
}

impl<F> BSplineWeights<8> for BSpline7<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        4.0.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 8] {
        debug_assert!(t >= 3.0.as_() && t <= 4.0.as_());

        let q = t - 3.0.as_();
        let s = 4.0.as_() - t;
        let q2 = q * q;
        let s2 = s * s;
        let qs = q * s;

        let h = qs.fmuladd(7.0.as_(), 21.0.as_());
        let k = qs.fmuladd(-21.0.as_(), -84.0.as_());
        let thirty_five_t = t * 35.0.as_();

        [
            (s2 * s2) * (s2 * s),
            h.fmuladd(s, 35.0.as_())
                .fmuladd(s, 35.0.as_())
                .fmuladd(s, 21.0.as_())
                .fmuladd(s, 7.0.as_())
                .fmuladd(s, F::ONE),
            s.fmuladd(-21.0.as_(), k)
                .fmuladd(s2, 280.0.as_())
                .fmuladd(s, 504.0.as_())
                .fmuladd(s, 392.0.as_())
                .fmuladd(s, 120.0.as_()),
            q2.fmuladd(thirty_five_t - 245.0.as_(), 560.0.as_())
                .fmuladd(q2, -1680.0.as_())
                .fmuladd(q2, 2416.0.as_()),
            s2.fmuladd(-thirty_five_t, 560.0.as_())
                .fmuladd(s2, -1680.0.as_())
                .fmuladd(s2, 2416.0.as_()),
            q.fmuladd(-21.0.as_(), k)
                .fmuladd(q2, 280.0.as_())
                .fmuladd(q, 504.0.as_())
                .fmuladd(q, 392.0.as_())
                .fmuladd(q, 120.0.as_()),
            h.fmuladd(q, 35.0.as_())
                .fmuladd(q, 35.0.as_())
                .fmuladd(q, 21.0.as_())
                .fmuladd(q, 7.0.as_())
                .fmuladd(q, F::ONE),
            (q2 * q2) * (q2 * q),
        ]
    }
}

impl<F> BSplineWeights<9> for BSpline8<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        4.5.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 9] {
        debug_assert!(t >= 3.5.as_() && t <= 4.5.as_());

        let q = t - 3.5.as_();
        let s = 4.5.as_() - t;
        let z = t - 4.0.as_();
        let z2 = z * z;
        let q2 = q * q;
        let s2 = s * s;
        let q4 = q2 * q2;
        let s4 = s2 * s2;
        let qs = q * s;

        let h = qs.fmuladd(8.0.as_(), 28.0.as_());
        let k = qs.fmuladd(-28.0.as_(), -140.0.as_());
        let l = qs.fmuladd(56.0.as_(), 252.0.as_());

        // Each pair uses the same Horner chain: u = s on the left, q on the right.
        let pair1 = |u: F| {
            h.fmuladd(u, 56.0.as_())
                .fmuladd(u, 70.0.as_())
                .fmuladd(u, 56.0.as_())
                .fmuladd(u, 28.0.as_())
                .fmuladd(u, 8.0.as_())
                .fmuladd(u, F::ONE)
        };

        let pair2 = |u: F| {
            u.fmuladd(-28.0.as_(), k)
                .fmuladd(u, -56.0.as_())
                .fmuladd(u, 490.0.as_())
                .fmuladd(u, 1288.0.as_())
                .fmuladd(u, 1540.0.as_())
                .fmuladd(u, 952.0.as_())
                .fmuladd(u, 247.0.as_())
        };

        let pair3 = |u: F| {
            u.fmuladd(112.0.as_(), l)
                .fmuladd(u, -504.0.as_())
                .fmuladd(u, -1890.0.as_())
                .fmuladd(u, -504.0.as_())
                .fmuladd(u, 5292.0.as_())
                .fmuladd(u, 8568.0.as_())
                .fmuladd(u, 4293.0.as_())
        };

        let center = z2
            .fmuladd(70.0.as_(), -630.0.as_())
            .fmuladd(z2, 3386.25.as_())
            .fmuladd(z2, -11379.375.as_())
            .fmuladd(z2, 18261.7734375.as_());

        [
            s4 * s4,
            pair1(s),
            pair2(s),
            pair3(s),
            center,
            pair3(q),
            pair2(q),
            pair1(q),
            q4 * q4,
        ]
    }
}

impl<F> BSplineWeights<10> for BSpline9<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        5.0.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 10] {
        debug_assert!(t >= 4.0.as_() && t <= 5.0.as_());

        let q = t - 4.0.as_();
        let s = 5.0.as_() - t;
        let q2 = q * q;
        let s2 = s * s;
        let q4 = q2 * q2;
        let s4 = s2 * s2;
        let qs = q * s;

        let h = qs.fmuladd(9.0.as_(), 36.0.as_());
        let k = qs.fmuladd(-36.0.as_(), -216.0.as_());
        let l = qs.fmuladd(84.0.as_(), 504.0.as_());

        let pair1 = |u: F| {
            h.fmuladd(u, 84.0.as_())
                .fmuladd(u, 126.0.as_())
                .fmuladd(u, 126.0.as_())
                .fmuladd(u, 84.0.as_())
                .fmuladd(u, 36.0.as_())
                .fmuladd(u, 9.0.as_())
                .fmuladd(u, F::ONE)
        };

        let pair2 = |u: F| {
            u.fmuladd(-36.0.as_(), k)
                .fmuladd(u, -168.0.as_())
                .fmuladd(u, 756.0.as_())
                .fmuladd(u, 2772.0.as_())
                .fmuladd(u, 4536.0.as_())
                .fmuladd(u, 4248.0.as_())
                .fmuladd(u, 2214.0.as_())
                .fmuladd(u, 502.0.as_())
        };

        let pair3 = |u: F| {
            u.fmuladd(168.0.as_(), l)
                .fmuladd(u, -672.0.as_())
                .fmuladd(u, -4284.0.as_())
                .fmuladd(u, -4032.0.as_())
                .fmuladd(u, 11256.0.as_())
                .fmuladd(u, 34272.0.as_())
                .fmuladd(u, 36414.0.as_())
                .fmuladd(u, 14608.0.as_())
        };

        let center = |u: F, u2: F| {
            u.fmuladd(-126.0.as_(), 630.0.as_())
                .fmuladd(u2, -4200.0.as_())
                .fmuladd(u2, 23940.0.as_())
                .fmuladd(u2, -88200.0.as_())
                .fmuladd(u2, 156190.0.as_())
        };

        [
            s4 * s4 * s,
            pair1(s),
            pair2(s),
            pair3(s),
            center(q, q2),
            center(s, s2),
            pair3(q),
            pair2(q),
            pair1(q),
            q4 * q4 * q,
        ]
    }
}

impl<F> BSplineWeights<11> for BSpline10<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        5.5.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 11] {
        debug_assert!(t >= 4.5.as_() && t <= 5.5.as_());

        let q = t - 4.5.as_();
        let s = 5.5.as_() - t;
        let z = t - 5.0.as_();
        let z2 = z * z;
        let q2 = q * q;
        let s2 = s * s;
        let q5 = q2 * q2 * q;
        let s5 = s2 * s2 * s;
        let qs = q * s;

        let h = qs.fmuladd(10.0.as_(), 45.0.as_());
        let k = qs.fmuladd(-45.0.as_(), -315.0.as_());
        let l = qs.fmuladd(120.0.as_(), 900.0.as_());
        let m = qs.fmuladd(-210.0.as_(), -1260.0.as_());

        let pair1 = |u: F| {
            h.fmuladd(u, 120.0.as_())
                .fmuladd(u, 210.0.as_())
                .fmuladd(u, 252.0.as_())
                .fmuladd(u, 210.0.as_())
                .fmuladd(u, 120.0.as_())
                .fmuladd(u, 45.0.as_())
                .fmuladd(u, 10.0.as_())
                .fmuladd(u, F::ONE)
        };

        let pair2 = |u: F| {
            u.fmuladd(-45.0.as_(), k)
                .fmuladd(u, -360.0.as_())
                .fmuladd(u, 1050.0.as_())
                .fmuladd(u, 5292.0.as_())
                .fmuladd(u, 11130.0.as_())
                .fmuladd(u, 14040.0.as_())
                .fmuladd(u, 11025.0.as_())
                .fmuladd(u, 5010.0.as_())
                .fmuladd(u, 1013.0.as_())
        };

        let pair3 = |u: F| {
            u.fmuladd(240.0.as_(), l)
                .fmuladd(u, -720.0.as_())
                .fmuladd(u, -8400.0.as_())
                .fmuladd(u, -13608.0.as_())
                .fmuladd(u, 16800.0.as_())
                .fmuladd(u, 100080.0.as_())
                .fmuladd(u, 171000.0.as_())
                .fmuladd(u, 141060.0.as_())
                .fmuladd(u, 47840.0.as_())
        };

        let pair4 = |u: F| {
            u.fmuladd(-630.0.as_(), m)
                .fmuladd(u, 5040.0.as_())
                .fmuladd(u, 16800.0.as_())
                .fmuladd(u, -13608.0.as_())
                .fmuladd(u, -119280.0.as_())
                .fmuladd(u, -95760.0.as_())
                .fmuladd(u, 327600.0.as_())
                .fmuladd(u, 736260.0.as_())
                .fmuladd(u, 455192.0.as_())
        };

        let center = z2
            .fmuladd(-252.0.as_(), 3465.0.as_())
            .fmuladd(z2, -30607.5.as_())
            .fmuladd(z2, 191585.625.as_())
            .fmuladd(z2, -769825.546875.as_())
            .fmuladd(z2, 1491301.23828125.as_());

        [
            s5 * s5,
            pair1(s),
            pair2(s),
            pair3(s),
            pair4(s),
            center,
            pair4(q),
            pair3(q),
            pair2(q),
            pair1(q),
            q5 * q5,
        ]
    }
}

impl<F> BSplineWeights<12> for BSpline11<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        6.0.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 12] {
        debug_assert!(t >= 5.0.as_() && t <= 6.0.as_());

        let q = t - 5.0.as_();
        let s = 6.0.as_() - t;
        let q2 = q * q;
        let s2 = s * s;
        let q5 = q2 * q2 * q;
        let s5 = s2 * s2 * s;
        let qs = q * s;

        let h = qs.fmuladd(11.0.as_(), 55.0.as_());
        let k = qs.fmuladd(-55.0.as_(), -440.0.as_());
        let l = qs.fmuladd(165.0.as_(), 1485.0.as_());
        let m = qs.fmuladd(-330.0.as_(), -2640.0.as_());

        let pair1 = |u: F| {
            h.fmuladd(u, 165.0.as_())
                .fmuladd(u, 330.0.as_())
                .fmuladd(u, 462.0.as_())
                .fmuladd(u, 462.0.as_())
                .fmuladd(u, 330.0.as_())
                .fmuladd(u, 165.0.as_())
                .fmuladd(u, 55.0.as_())
                .fmuladd(u, 11.0.as_())
                .fmuladd(u, F::ONE)
        };

        let pair2 = |u: F| {
            u.fmuladd(-55.0.as_(), k)
                .fmuladd(u, -660.0.as_())
                .fmuladd(u, 1320.0.as_())
                .fmuladd(u, 9240.0.as_())
                .fmuladd(u, 24024.0.as_())
                .fmuladd(u, 38280.0.as_())
                .fmuladd(u, 40260.0.as_())
                .fmuladd(u, 27500.0.as_())
                .fmuladd(u, 11132.0.as_())
                .fmuladd(u, 2036.0.as_())
        };

        let pair3 = |u: F| {
            u.fmuladd(330.0.as_(), l)
                .fmuladd(u, -495.0.as_())
                .fmuladd(u, -14850.0.as_())
                .fmuladd(u, -34650.0.as_())
                .fmuladd(u, 12474.0.as_())
                .fmuladd(u, 236610.0.as_())
                .fmuladd(u, 586575.0.as_())
                .fmuladd(u, 748275.0.as_())
                .fmuladd(u, 515097.0.as_())
                .fmuladd(u, 152637.0.as_())
        };

        let pair4 = |u: F, u2: F| {
            u.fmuladd(-990.0.as_(), m)
                .fmuladd(u, 7920.0.as_())
                .fmuladd(u, 39600.0.as_())
                .fmuladd(u2, -299376.0.as_()) // Skips the missing power.
                .fmuladd(u, -538560.0.as_())
                .fmuladd(u, 574200.0.as_())
                .fmuladd(u, 3273600.0.as_())
                .fmuladd(u, 4480872.0.as_())
                .fmuladd(u, 2203488.0.as_())
        };

        let center = |u: F, u2: F| {
            u.fmuladd(462.0.as_(), -2772.0.as_())
                .fmuladd(u2, 27720.0.as_())
                .fmuladd(u2, -255024.0.as_())
                .fmuladd(u2, 1718640.0.as_())
                .fmuladd(u2, -7475160.0.as_())
                .fmuladd(u2, 15724248.0.as_())
        };

        [
            s5 * s5 * s,
            pair1(s),
            pair2(s),
            pair3(s),
            pair4(s, s2),
            center(q, q2),
            center(s, s2),
            pair4(q, q2),
            pair3(q),
            pair2(q),
            pair1(q),
            q5 * q5 * q,
        ]
    }
}

impl<F> BSplineWeights<13> for BSpline12<F>
where
    F: Float + 'static + ConstOne + ConstZero,
    f64: AsPrimitive<F>,
{
    #[inline(always)]
    fn radius() -> Self::Float {
        6.5.as_()
    }

    #[inline(always)]
    fn weights(t: F) -> [F; 13] {
        debug_assert!(t >= 5.5.as_() && t <= 6.5.as_());

        let q = t - 5.5.as_();
        let s = 6.5.as_() - t;
        let z = t - 6.0.as_();
        let z2 = z * z;
        let q3 = q * q * q;
        let s3 = s * s * s;
        let q6 = q3 * q3;
        let s6 = s3 * s3;
        let qs = q * s;

        let h1 = qs.fmuladd(12.0.as_(), 66.0.as_());
        let h2 = qs.fmuladd(-66.0.as_(), -594.0.as_());
        let h3 = qs.fmuladd(220.0.as_(), 2310.0.as_());
        let h4 = qs.fmuladd(-495.0.as_(), -4950.0.as_());
        let h5 = qs.fmuladd(792.0.as_(), 5940.0.as_());

        let pair1 = |u: F| {
            h1.fmuladd(u, 220.0.as_())
                .fmuladd(u, 495.0.as_())
                .fmuladd(u, 792.0.as_())
                .fmuladd(u, 924.0.as_())
                .fmuladd(u, 792.0.as_())
                .fmuladd(u, 495.0.as_())
                .fmuladd(u, 220.0.as_())
                .fmuladd(u, 66.0.as_())
                .fmuladd(u, 12.0.as_())
                .fmuladd(u, F::ONE)
        };

        let pair2 = |u: F| {
            u.fmuladd(-66.0.as_(), h2)
                .fmuladd(u, -1100.0.as_())
                .fmuladd(u, 1485.0.as_())
                .fmuladd(u, 15048.0.as_())
                .fmuladd(u, 47124.0.as_())
                .fmuladd(u, 91080.0.as_())
                .fmuladd(u, 120285.0.as_())
                .fmuladd(u, 109780.0.as_())
                .fmuladd(u, 66726.0.as_())
                .fmuladd(u, 24420.0.as_())
                .fmuladd(u, 4083.0.as_())
        };

        let pair3 = |u: F| {
            u.fmuladd(440.0.as_(), h3)
                .fmuladd(u, 220.0.as_())
                .fmuladd(u, -24255.0.as_())
                .fmuladd(u, -75240.0.as_())
                .fmuladd(u, -23100.0.as_())
                .fmuladd(u, 475992.0.as_())
                .fmuladd(u, 1638945.0.as_())
                .fmuladd(u, 2883100.0.as_())
                .fmuladd(u, 3023790.0.as_())
                .fmuladd(u, 1807212.0.as_())
                .fmuladd(u, 478271.0.as_())
        };

        let pair4 = |u: F| {
            u.fmuladd(-1485.0.as_(), h4)
                .fmuladd(u, 11220.0.as_())
                .fmuladd(u, 81675.0.as_())
                .fmuladd(u, 59400.0.as_())
                .fmuladd(u, -623700.0.as_())
                .fmuladd(u, -1860408.0.as_())
                .fmuladd(u, -37125.0.as_())
                .fmuladd(u, 10101300.0.as_())
                .fmuladd(u, 23794650.0.as_())
                .fmuladd(u, 24610212.0.as_())
                .fmuladd(u, 10187685.0.as_())
        };

        let pair5 = |u: F| {
            u.fmuladd(3168.0.as_(), h5)
                .fmuladd(u, -38280.0.as_())
                .fmuladd(u, -121770.0.as_())
                .fmuladd(u, 261360.0.as_())
                .fmuladd(u, 1524600.0.as_())
                .fmuladd(u, -137808.0.as_())
                .fmuladd(u, -10932570.0.as_())
                .fmuladd(u, -14341800.0.as_())
                .fmuladd(u, 29555460.0.as_())
                .fmuladd(u, 90415512.0.as_())
                .fmuladd(u, 66318474.0.as_())
        };

        let center = z2
            .fmuladd(924.0.as_(), -18018.0.as_())
            .fmuladd(z2, 236486.25.as_())
            .fmuladd(z2, -2358105.75.as_())
            .fmuladd(z2, 17075574.140625.as_())
            .fmuladd(z2, -79943684.1328125.as_())
            .fmuladd(z2, 181466922.6474609375.as_());

        [
            s6 * s6,
            pair1(s),
            pair2(s),
            pair3(s),
            pair4(s),
            pair5(s),
            center,
            pair5(q),
            pair4(q),
            pair3(q),
            pair2(q),
            pair1(q),
            q6 * q6,
        ]
    }
}

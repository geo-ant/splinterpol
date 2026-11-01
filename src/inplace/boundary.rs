#[cfg(test)]
mod test;

/// Describes a way to extend the signal beyond the measured domain. In
/// particular these particular conditions is what B&M call _transmitted_
/// c.f section 3.2.2 in B&M and Algorithm 5 therein. These kinds of boundary
/// conditions enable calculations to be performed completely _in-place_
/// in the orignal image.
///
/// Note, for example, that constant extension can _not_ be written as a
/// transmittable boundary condition. Only transmittable boundary conditions
/// can be used for in-place coefficient calculation and interpolation.
///
/// ## References
///
/// * [B&M](https://www.ipol.im/pub/art/2018/221/) T. Briand, P. Monasse: "Theory
/// and Practice of Image B-Spline Interpolation"
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum InplaceBoundaryExtension {
    /// Signal `abcde` is extended to `...cba|abcde|edc` and is called
    /// _half-symmetric_ in the Briand & Monasse paper. In `scipy`
    /// this is also called
    /// [`reflected`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter).
    Reflected = 1,
    /// Signal `abcde` is extended to `...dcb|abcde|dcb...`. In `scipy` this is
    /// is called
    /// [`mirrored`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter)
    /// and is the default mode for [`spline_filter`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter).
    WholeSymmetric = 2,
    /// Signal `abcde` is extended to `...cde|abcde|abc...`. In `scipy` this
    /// is called [`grid-wrap`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter).
    Periodic = 3,
}

impl Default for InplaceBoundaryExtension {
    /// We chose the same default here as the `scipy` function
    /// [`spline_filter`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.spline_filter.html#scipy.ndimage.spline_filter).
    /// This doesn't mean that this is the best continuation, your mileage may vary.
    fn default() -> Self {
        Self::WholeSymmetric
    }
}

impl InplaceBoundaryExtension {
    #[inline]
    pub(crate) fn extended_index(self, idx: isize, n: usize) -> usize {
        debug_assert!(n > 0, "logic error: index extension for zero length row");

        //TODO PERF: the reference code uses repeated subtraction, which might
        // be the smarter choice compared to modulo, since likely there's only
        // one subtraction to be performed, but we'd have to measure that.

        let n = n as isize;

        let mapped = match self {
            Self::Reflected => {
                let period = 2 * n;
                let i = idx.rem_euclid(period);

                if i < n { i } else { period - 1 - i }
            }

            Self::WholeSymmetric => {
                if n == 1 {
                    return 0;
                }

                let period = 2 * (n - 1);
                let i = idx.rem_euclid(period);

                if i < n { i } else { period - i }
            }

            Self::Periodic => idx.rem_euclid(n),
        };

        mapped as usize
    }
}

#[cfg(test)]
mod test;

/// Describes a way to extend the signal beyond the measured domain in a way
/// that can be _transmitted_ like described in section 3.2.2 of the
/// Briand & Monasse paper (cf also Algorithm 5 therein). Note that constant
/// extension can _not_ be written as a transmittable boundary condition.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum TransmittableBoundaryExtension {
    /// Signal `abcde` is extended to `...cba|abcde|edc` and is called
    /// _half-symmetric_ in the Briand & Monasse paper.
    Mirrored = 1,
    /// Signal `abcde` is extended to `...dcb|abcde|dcb...`
    WholeSymmetric = 2,
    /// Signal `abcde` is extended to `...cde|abcde|abc...`
    Periodic = 3,
}

impl TransmittableBoundaryExtension {
    #[inline]
    pub(crate) fn extended_index(self, idx: isize, n: usize) -> usize {
        debug_assert!(n > 0, "logic error: index extension for zero length row");

        //TODO PERF: the reference code uses repeated subtraction, which might
        // be the smarter choice compared to modulo, since likely there's only
        // one subtraction to be performed, but we'd have to measure that.

        let n = n as isize;

        let mapped = match self {
            Self::Mirrored => {
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

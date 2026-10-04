use crate::{Error, TransmittableBoundaryExtension, bspline::BSpline, numerics::FusedMuladd};
use core::iter::Sum;
use num_traits::{ConstOne, ConstZero, Float};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
/// This is purely a convenience type that acts as a builder to make the
/// calls to [`crate::inplace::coefficients2d`] and [`crate::inplace::coefficients2d_blocked`]
/// less awkward to write. It performs _no setup at all_ and only acts as pure
/// syntax sugar to the method calls.
///
/// ## Examples
///
/// ### Without Blocking
///
/// To replace a call to [`crate::inplace::coefficients2d`] use it as
///
///
/// ```rust
/// # use spleen::*;
/// # use spleen::bspline::*;
/// # let mut some_data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = 9;
/// # let height = 10;
/// # let extension = spleen::TransmittableBoundaryExtension::Periodic;
/// let mut data1 = some_data;
/// let mut data2 = some_data;
///
/// // using the rather awkward function call
/// inplace::coefficients2d::<_,_>(BSpline3::<f64>::new(),
///                              &mut data1,
///                              width,
///                              height,
///                              epsilon,
///                              extension).unwrap();
///
///
/// // using the builder
/// inplace::Coefficients2d::new(BSpline3::<f64>::new())
///     .compute(
///         &mut data2,
///         width,
///         height,
///         epsilon,
///         extension,
///     ).unwrap();
///
/// // produces the same coefficients in data1 and data2
/// assert_eq!(data1,data2);
/// ```
///
/// ### With Blocking
///
/// The example above uses block size of `1` for both rows and columns. To explicily
/// set a non-unit block size for rows and/or columns, use the builder like so:
///
/// ```rust
/// # use spleen::*;
/// # use spleen::bspline::*;
/// # let mut some_data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = 9;
/// # let height = 10;
/// # let extension = spleen::TransmittableBoundaryExtension::Periodic;
/// let mut data1 = some_data;
/// let mut data2 = some_data;
///
/// // using the rather awkward function call
/// inplace::coefficients2d_blocked::<_,_,2,3>(BSpline3::<f64>::new(),
///                              &mut data1,
///                              width,
///                              height,
///                              epsilon,
///                              extension).unwrap();
///
///
/// // using the builder
/// inplace::Coefficients2d::new(BSpline3::<f64>::new())
///     .with_row_blocksize::<2>()
///     .with_column_blocksize::<3>()
///     .compute(
///         &mut data2,
///         width,
///         height,
///         epsilon,
///         extension,
///     ).unwrap();
///
/// // produces the same coefficients in data1 and data2
/// assert_eq!(data1,data2);
/// ```
pub struct Coefficients2d<
    // number of poles of the spline
    const NPOLES: usize,
    // bspline to use for interpolationj
    B: BSpline<NPOLES>,
    // row-blocksize for calculations
    const RB: usize,
    // column-blocksize for calculations
    const CB: usize,
> {
    bspline: B,
}

impl<const NPOLES: usize, B: BSpline<NPOLES>> Coefficients2d<NPOLES, B, 1, 1> {
    #[inline(always)]
    /// create a new builder for computing the coefficients of the given spline
    /// in-place
    pub fn new(bspline: B) -> Self {
        Self { bspline }
    }
}

impl<const NPOLES: usize, B: BSpline<NPOLES>, const RB: usize> Coefficients2d<NPOLES, B, RB, 1> {
    #[inline(always)]
    /// set the column blocksize for calculation to this nonzero compile time value
    pub fn with_column_blocksize<const CB: usize>(self) -> Coefficients2d<NPOLES, B, RB, CB> {
        const { assert!(CB > 0) };
        let Self { bspline } = self;
        Coefficients2d { bspline }
    }
}

impl<const NPOLES: usize, B: BSpline<NPOLES>, const CB: usize> Coefficients2d<NPOLES, B, 1, CB> {
    #[inline(always)]
    /// set the row blocksize for calculation to this nonzero compile time value
    pub fn with_row_blocksize<const RB: usize>(self) -> Coefficients2d<NPOLES, B, RB, CB> {
        const { assert!(RB > 0) };
        let Self { bspline } = self;
        Coefficients2d { bspline }
    }
}

impl<const NPOLES: usize, B: BSpline<NPOLES>, const RB: usize, const CB: usize>
    Coefficients2d<NPOLES, B, RB, CB>
{
    #[inline(always)]
    /// perform the computation. Again, this is just syntax sugar for
    /// [`crate::inplace::coefficients2d_blocked`] calls with the same arguments.
    /// There is no advantage or disadvantage in terms of computation in using
    /// this builder pattern. The only advantage is nicer syntax.
    pub fn compute(
        self,
        image_data: &mut [B::Float],
        width: usize,
        height: usize,
        epsilon: B::Float,
        boundary_extension: TransmittableBoundaryExtension,
    ) -> Result<(), Error>
    where
        B::Float: Float + ConstOne + ConstZero + Sum + FusedMuladd,
        B: BSpline<NPOLES>,
    {
        super::coefficients2d_blocked::<_, _, RB, CB>(
            self.bspline,
            image_data,
            width,
            height,
            epsilon,
            boundary_extension,
        )
    }
}

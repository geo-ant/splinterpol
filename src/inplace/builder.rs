use crate::{
    Error, RowMajor,
    bspline::{BSplinePrefilter, BSplineWeights},
    inplace::{TransmittableBoundaryExtension, interp2d_single},
    numerics::FusedMuladd,
};
use core::{iter::Sum, num::NonZero};
use num_traits::{ConstOne, ConstZero, Float, FromPrimitive, NumCast, ToPrimitive};

#[derive(Clone, PartialEq, Debug)]
pub struct SplineInterpolator2d<B, C> {
    // NOTE(geo) this guy uses owned containers because forcing a slice here
    // creates a mildly (but not totally, but enough for me to change it)
    // API when creating an instance of an interpolator, see
    // https://doc.rust-lang.org/nomicon/lifetime-mismatch.html
    // That would mean that we'd have to accept a mutable borrow to the input
    // but we want to return a shared borrow to the input on exit. That works
    // and we could create shared references to the data through this type,
    // but it would make the mutable borrow to the input stay alive as long
    // as this type lives which could create surprises to downstream users.
    coeffs: RowMajor<C>,
    bspline: B,
    boundary_extension: TransmittableBoundaryExtension,
}

impl<B, C> SplineInterpolator2d<B, C> {
    pub fn from_image<const NPOLES: usize>(
        bspline: B,
        image_data: RowMajor<C>,
        boundary_extension: TransmittableBoundaryExtension,
        epsilon: B::Float,
    ) -> Result<Self, Error>
    where
        C: AsMut<[B::Float]>,
        B: BSplinePrefilter<NPOLES>,
        B::Float: Float + ConstOne + ConstZero + Sum,
    {
        let RowMajor {
            mut data,
            width,
            height,
        } = image_data;
        crate::inplace::coefficients2d_blocked::<_, _, 1, 1>(
            bspline,
            RowMajor {
                data: data.as_mut(),
                width,
                height,
            },
            epsilon,
            boundary_extension,
        )?;
        Ok(Self {
            coeffs: RowMajor {
                data,
                width,
                height,
            },
            bspline,
            boundary_extension,
        })
    }

    pub fn interp_single<const NWEIGHTS: usize>(
        &self,
        row: B::Float,
        col: B::Float,
    ) -> Result<B::Float, Error>
    where
        C: AsRef<[B::Float]>,
        B: BSplineWeights<NWEIGHTS>,
        B::Float: Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast,
    {
        let width = self.coeffs.width;
        let height = self.coeffs.height;
        interp2d_single(
            self.bspline,
            RowMajor {
                data: self.coeffs.as_ref(),
                width,
                height,
            },
            self.boundary_extension,
            (row, col),
        )
    }
}

impl<C, B> SplineInterpolator2d<B, C> {
    /// consumes the interpolator and returns the calculated
    /// B-Spline coefficients in row major fashion
    pub fn coefficients(self) -> RowMajor<C> {
        self.coeffs
    }
}

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
/// # use spleen::RowMajor;
/// # use core::num::NonZero;
/// # let mut some_data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = NonZero::new(9).unwrap();
/// # let height = NonZero::new(10).unwrap();
/// # let extension = spleen::inplace::TransmittableBoundaryExtension::Periodic;
/// let mut data1 = some_data;
/// let mut data2 = some_data;
///
/// // using the rather awkward function call
/// inplace::coefficients2d::<_,_>(BSpline3::new(),
///                              RowMajor {
///                                data: &mut data1,
///                                width,
///                                height,
///                              },
///                              epsilon,
///                              extension).unwrap();
///
///
/// // using the builder
/// inplace::Coefficients2d::new(BSpline3::new())
///     .compute(
///         RowMajor {
///           data: &mut data2,
///           width,
///           height,
///         },
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
/// # fn foo() {
/// # use spleen::*;
/// # use spleen::bspline::*;
/// # let mut some_data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = core::num::NonZero::new(9).unwrap();
/// # let height = core::num::NonZero::new(10).unwrap();
/// # let extension = spleen::inplace::TransmittableBoundaryExtension::Periodic;
/// let mut data1 = some_data;
/// let mut data2 = some_data;
///
/// // using the rather awkward function call
/// inplace::coefficients2d_blocked::<_,_,2,3>(BSpline3::<f64>::new(),
///                              RowMajor {
///                                data: &mut data1,
///                                width,
///                                height,
///                              },
///                              epsilon,
///                              extension).unwrap();
///
///
/// // using the builder
/// inplace::Coefficients2d::new(BSpline3::<f64>::new())
///     .with_row_blocksize::<2>()
///     .with_column_blocksize::<3>()
///     .compute(
///         RowMajor {
///           data: &mut data2,
///           width,
///           height,
///         },
///         epsilon,
///         extension,
///     ).unwrap();
///
/// // produces the same coefficients in data1 and data2
/// assert_eq!(data1,data2);
/// # }
/// ```
pub struct Coefficients2d<
    // number of poles of the spline
    const NPOLES: usize,
    // bspline to use for interpolationj
    B: BSplinePrefilter<NPOLES>,
    // row-blocksize for calculations
    const RB: usize,
    // column-blocksize for calculations
    const CB: usize,
> {
    bspline: B,
}

impl<const NPOLES: usize, B: BSplinePrefilter<NPOLES>> Coefficients2d<NPOLES, B, 1, 1> {
    #[inline(always)]
    /// create a new builder for computing the coefficients of the given spline
    /// in-place
    pub fn new(bspline: B) -> Self {
        Self { bspline }
    }
}

impl<const NPOLES: usize, B: BSplinePrefilter<NPOLES>, const RB: usize>
    Coefficients2d<NPOLES, B, RB, 1>
{
    #[inline(always)]
    /// set the column blocksize for calculation to this nonzero compile time value
    pub fn with_column_blocksize<const CB: usize>(self) -> Coefficients2d<NPOLES, B, RB, CB> {
        const { assert!(CB > 0) };
        let Self { bspline } = self;
        Coefficients2d { bspline }
    }
}

impl<const NPOLES: usize, B: BSplinePrefilter<NPOLES>, const CB: usize>
    Coefficients2d<NPOLES, B, 1, CB>
{
    #[inline(always)]
    /// set the row blocksize for calculation to this nonzero compile time value
    pub fn with_row_blocksize<const RB: usize>(self) -> Coefficients2d<NPOLES, B, RB, CB> {
        const { assert!(RB > 0) };
        let Self { bspline } = self;
        Coefficients2d { bspline }
    }
}

impl<const NPOLES: usize, B: BSplinePrefilter<NPOLES>, const RB: usize, const CB: usize>
    Coefficients2d<NPOLES, B, RB, CB>
{
    #[inline(always)]
    /// perform the computation. Again, this is just syntax sugar for
    /// [`crate::inplace::coefficients2d_blocked`] calls with the same arguments.
    /// There is no advantage or disadvantage in terms of computation in using
    /// this builder pattern. The only advantage is nicer syntax.
    pub fn compute(
        self,
        image_data: RowMajor<&mut [B::Float]>,
        epsilon: B::Float,
        boundary_extension: TransmittableBoundaryExtension,
    ) -> Result<(), Error>
    where
        B::Float: Float + ConstOne + ConstZero + Sum + FusedMuladd,
        B: BSplinePrefilter<NPOLES>,
    {
        super::coefficients2d_blocked::<_, _, RB, CB>(
            self.bspline,
            image_data,
            epsilon,
            boundary_extension,
        )
    }
}

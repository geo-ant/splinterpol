use crate::{
    Error, RowMajor,
    bspline::{BSpline, BSplinePrefilter, BSplineWeights},
    inplace::{InplaceBoundaryExtension, interp2d_single},
    numerics::FusedMuladd,
};
use core::iter::Sum;
use num_traits::{ConstOne, ConstZero, Float, FromPrimitive, NumCast, ToPrimitive};

#[cfg(test)]
mod test;

#[derive(Copy, Clone, PartialEq, Eq)]
/// Higher-level interface to spline-interpolation of 2D images. Construct
/// an instance of this by using the builder pattern in the example below.
/// The interpolator takes ownership of the image data and the builder
/// pattern guides you through the 2-step process of coefficient calculation
/// and interpolation.
///
/// ## Example
///
/// ```rust
/// # fn foo() {
/// # use splinterpol::*;
/// use splinterpol::bspline::*;
/// use splinterpol::inplace::SplineInterpolator2d;
/// use splinterpol::inplace::InplaceBoundaryExtension::*;
/// # let mut data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = core::num::NonZero::new(9).unwrap();
/// # let height = core::num::NonZero::new(10).unwrap();
/// # let extension = splinterpol::inplace::InplaceBoundaryExtension::Periodic;
/// # let image_data = RowMajor {data, width, height};
///
/// let interp = SplineInterpolator2d::new(BSpline3::new())
///                // this is optional and defaults
///                // to whole-symmetric extension
///                .extension(WholeSymmetric)
///                // giving the block sizes for
///                // computation is optional and
///                // defaults to `1`. Block sizes
///                // are in number of elements, not
///                // bytes.
///                .row_blocksize::<8>()
///                .col_blocksize::<1>()
///                // this is mandatory and finalizes
///                // the computation of the interpolator
///                // after which it can be used for
///                // interpolation.
///                .compute(image_data, epsilon)
///                .unwrap();
///
/// // now use the interpolator to find the values at
/// // arbitraty points inside the original index grid
/// let value = interp.interp_single(1.5, 2.4);
/// # }
/// ```
pub struct SplineInterpolator2d<B, C> {
    bspline: B,
    // NOTE(geo) this guy uses owned containers because forcing a slice here
    // creates a mildly (but not totally, but enough for me to change it)
    // API when creating an instance of an interpolator, see
    // https://doc.rust-lang.org/nomicon/lifetime-mismatch.html
    // That would mean that we'd have to accept a mutable borrow to the input
    // but we want to return a shared borrow to the input on exit. That works
    // and we could create shared references to the data through this type,
    // but it would make the mutable borrow to the input stay alive as long
    // as this type lives which could create surprises to downstream users.
    //
    // NOTE: The reason we use () for the container initially is that we
    // make the construction of this guy follow a builder pattern, which
    // makes for some pretty ugly generic, with very nice usage.
    coeffs: RowMajor<C>,
    boundary_extension: InplaceBoundaryExtension,
}

/// This is just an intermediate structure on the way to configuring and
/// building a spline interpolator.
pub struct SplineInterpBuilder2d<B, C, const RB: usize, const CB: usize> {
    bspline: B,
    coeffs: C,
    boundary_extension: InplaceBoundaryExtension,
}

impl SplineInterpolator2d<(), ()> {
    #[inline]
    /// Starts the builder process for a spline interpolator, which
    /// is finalized with the call to [`compute`](SplineInterpBuilder2d::compute)
    pub fn new<B: BSpline>(bspline: B) -> SplineInterpBuilder2d<B, (), 1, 1> {
        SplineInterpBuilder2d {
            bspline,
            coeffs: (),
            boundary_extension: InplaceBoundaryExtension::WholeSymmetric,
        }
    }
}

impl<B, C, const CB: usize> SplineInterpBuilder2d<B, C, 1, CB> {
    #[inline]
    /// Set the block size for column calculations at compile time, see
    /// the comments on [`crate::inplace::coefficients2d_blocked`].
    pub fn col_blocksize<const RB: usize>(self) -> SplineInterpBuilder2d<B, C, RB, CB> {
        let SplineInterpBuilder2d {
            bspline,
            coeffs,
            boundary_extension,
        } = self;
        SplineInterpBuilder2d {
            bspline,
            coeffs,
            boundary_extension,
        }
    }
}

impl<B, C, const RB: usize> SplineInterpBuilder2d<B, C, RB, 1> {
    #[inline]
    /// Set the block size for row calculations at compile time, see
    /// the comments on [`crate::inplace::coefficients2d_blocked`].
    pub fn row_blocksize<const CB: usize>(self) -> SplineInterpBuilder2d<B, C, RB, CB> {
        let SplineInterpBuilder2d {
            bspline,
            coeffs,
            boundary_extension,
        } = self;
        SplineInterpBuilder2d {
            bspline,
            coeffs,
            boundary_extension,
        }
    }
}

impl<B, const RB: usize, const CB: usize> SplineInterpBuilder2d<B, (), RB, CB> {
    /// Set the boundary extension method used by the interpolator. If this
    /// method is not invoked, the whole-symmetric extension is used. Obviously,
    /// this must be called before the call to [`compute`](Self::compute).
    pub fn extension(self, boundary_extension: InplaceBoundaryExtension) -> Self {
        Self {
            boundary_extension,
            ..self
        }
    }
}

impl<B, const RB: usize, const CB: usize> SplineInterpBuilder2d<B, (), RB, CB> {
    /// Try and perform the necessary precomputations (i.e. coefficient calculations)
    /// to obtain a usable interpolator which we can use to interpolate the
    /// given `image_data`. To understand the meaning of `epsilon` and its
    /// recommended values, see the documentation of [`crate::inplace::coefficients2d_blocked`].
    /// In short, `epsilon` bounds the _relative_ accuracy of the interpolation
    /// _on the original grid points_. A conservative rule of thumb is to chose
    /// `1e-12` for `f64` data and `1e-6` for `f32`.
    pub fn compute<const NPOLES: usize, C>(
        self,
        image_data: RowMajor<C>,
        epsilon: B::Float,
    ) -> Result<SplineInterpolator2d<B, C>, Error>
    where
        C: AsMut<[B::Float]>,
        B: BSplinePrefilter<NPOLES>,
        B::Float: Float + ConstOne + ConstZero + Sum,
    {
        let Self {
            bspline,
            coeffs: _coeffs,
            boundary_extension,
        } = self;
        let RowMajor {
            mut data,
            width,
            height,
        } = image_data;
        crate::inplace::coefficients2d_blocked::<_, _, RB, CB>(
            bspline,
            RowMajor {
                data: data.as_mut(),
                width,
                height,
            },
            epsilon,
            boundary_extension,
        )?;
        Ok(SplineInterpolator2d {
            coeffs: RowMajor {
                data,
                width,
                height,
            },
            bspline,
            boundary_extension,
        })
    }
}

impl<B, C> SplineInterpolator2d<B, C> {
    /// Interpolate the value of a single point at the given floating point
    /// indices for the `row` and `column` of the image.
    ///
    /// **Attention**: this function indexes `(row, col)`, which
    /// is `(y,x)` in the common `x,y` coordinate system.
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

impl<B, C, const RB: usize, const CB: usize> SplineInterpBuilder2d<B, RowMajor<C>, RB, CB>
where
    B: BSpline,
{
    /// Consume the interpolator and returns the calculated
    /// B-Spline coefficients in row major fashion
    pub fn coefficients(self) -> RowMajor<C> {
        self.coeffs
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
/// This is purely a convenience type that acts as a builder to make the
/// calls to [`crate::inplace::coefficients2d`] and [`crate::inplace::coefficients2d_blocked`]
/// less awkward to write. It performs _no setup at all_ and only acts as pure
/// syntax sugar to the method calls. There is no advantage or disadvantage
/// in terms of computation to using this pattern for coefficient calculation.
/// The only advantage is nicer syntax.
///
/// ## Examples
///
/// ### Without Blocking
///
/// To replace a call to [`crate::inplace::coefficients2d`] use it as
///
///
/// ```rust
/// # use splinterpol::*;
/// # use splinterpol::bspline::*;
/// # use splinterpol::RowMajor;
/// # use core::num::NonZero;
/// # let mut some_data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = NonZero::new(9).unwrap();
/// # let height = NonZero::new(10).unwrap();
/// # let extension = splinterpol::inplace::InplaceBoundaryExtension::Periodic;
/// let mut data1 = some_data;
/// let mut data2 = some_data;
///
/// // using the rather awkward function call with unit block sizes
/// inplace::coefficients2d_blocked::<_,_,1,1>(BSpline3::new(),
///                              RowMajor {
///                                data: &mut data1,
///                                width,
///                                height,
///                              },
///                              epsilon,
///                              extension).unwrap();
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
/// # use splinterpol::*;
/// # use splinterpol::bspline::*;
/// # let mut some_data = [10.;90];
/// # let epsilon = 1e-1;
/// # let width = core::num::NonZero::new(9).unwrap();
/// # let height = core::num::NonZero::new(10).unwrap();
/// # let extension = splinterpol::inplace::InplaceBoundaryExtension::Periodic;
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
    /// Set the column blocksize for calculation to this nonzero compile time value.
    /// Defaults to `1` if not given.
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
    /// Set the row blocksize for calculation to this nonzero compile time value
    /// Defaults to `1` if not given.
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
    /// Please consult its documentation for the meaning of the parameters
    /// and some guidance.
    pub fn compute(
        self,
        image_data: RowMajor<&mut [B::Float]>,
        epsilon: B::Float,
        boundary_extension: InplaceBoundaryExtension,
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

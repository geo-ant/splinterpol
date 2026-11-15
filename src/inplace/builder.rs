use crate::{
    Error, RowMajor,
    bspline::{BSpline, BSplinePrefilter, BSplineWeights},
    inplace::{InplaceBoundaryExtension, interp2d_single, interpolation::interp1d_single},
    shared::{ColBlockSize, RowBlockSize},
};
use core::{iter::Sum, marker::PhantomData};
use num_traits::{ConstOne, ConstZero, Float, FromPrimitive, NumCast, ToPrimitive};

#[cfg(test)]
mod test;

#[derive(Copy, Clone, PartialEq, Eq)]
/// High-level interface to spline-interpolation of 2D images. Construct
/// an instance of this by using the builder pattern in the example below.
/// The interpolator takes ownership of the image data and the builder
/// pattern guides you through the 2-step process of coefficient calculation
/// and interpolation.
///
/// ## Example: General Purpose Interpolation
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
/// let interp = SplineInterpolator2d::builder(BSpline3::new())
///                // this is optional and defaults to whole-symmetric extension
///                .extension(WholeSymmetric)
///                // giving the block sizes for computation is optional and
///                // defaults to `1`. Block sizes are in number of elements,
///                // not bytes.
///                .row_blocksize::<8>()
///                .col_blocksize::<1>()
///                // this is mandatory and finalizes the builder by the
///                // computing the spline coefficients. After this, we can
///                // perform interpolation
///                .build(image_data, epsilon)
///                .unwrap();
///
/// // now use the interpolator to find the values at arbitraty points in
/// // the image
///
/// let value = interp.interp_single(1.5, 2.4);
/// # }
/// ```
pub struct SplineInterpolator2d<B, C> {
    phantom: PhantomData<B>,
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

/// High-level interface to spline-interpolation of 1D data. It's the 1D
/// analog to [`SplineInterpolator2d`] and functions very much the same way.
///
/// ## Example: General Purpose Interpolation
///
/// ```rust
/// # fn foo() {
/// # use splinterpol::*;
/// use splinterpol::bspline::*;
/// use splinterpol::inplace::SplineInterpolator1d;
/// use splinterpol::inplace::InplaceBoundaryExtension::*;
/// # let mut data = [10.;90];
/// # let epsilon = 1e-1;
/// # let extension = splinterpol::inplace::InplaceBoundaryExtension::Periodic;
///
/// let interp = SplineInterpolator1d::builder(BSpline3::new())
///                // this is optional and defaults to whole-symmetric extension
///                .extension(WholeSymmetric)
///                // giving the block size for computation is optional and
///                // defaults to `1`. Block size is in number of elements,
///                // not bytes.
///                .row_blocksize::<8>()
///                // this is mandatory and finalizes the builder by the
///                // computing the spline coefficients. After this, we can
///                // perform interpolation
///                .build(data, epsilon)
///                .unwrap();
///
/// // now use the interpolator to find the values at an arbitraty position
///
/// let value = interp.interp_single(1.5);
/// # }
/// ```
pub struct SplineInterpolator1d<B, C> {
    phantom: PhantomData<B>,
    coeffs: C,
    boundary_extension: InplaceBoundaryExtension,
}

/// This is just an intermediate structure on the way to configuring and
/// building a spline interpolator in 1d or 2d.
// NOTE: we're using `DIM` here to differentiate between 1d and 2d code where
// necessary, and share the other code otherwise.
pub struct SplineInterpBuilder<const DIM: usize, B, C, const RB: usize, const CB: usize> {
    bspline: B,
    coeffs: C,
    boundary_extension: InplaceBoundaryExtension,
}

impl SplineInterpolator2d<(), ()> {
    #[inline]
    /// Starts the builder process for a spline interpolator, which
    /// is finalized with the call to the 2D specialization of the `build` method.
    pub fn builder<B: BSpline>(bspline: B) -> SplineInterpBuilder<2, B, (), 1, 1> {
        SplineInterpBuilder {
            bspline,
            coeffs: (),
            boundary_extension: InplaceBoundaryExtension::WholeSymmetric,
        }
    }
}

impl SplineInterpolator1d<(), ()> {
    #[inline]
    /// Starts the builder process for a spline interpolator, which
    /// is finalized with the call to the 1D specialization of the `build` method.
    pub fn builder<B: BSpline>(bspline: B) -> SplineInterpBuilder<1, B, (), 1, 1> {
        SplineInterpBuilder {
            bspline,
            coeffs: (),
            boundary_extension: InplaceBoundaryExtension::WholeSymmetric,
        }
    }
}

impl<B, C, const CB: usize> SplineInterpBuilder<2, B, C, 1, CB> {
    #[inline]
    /// Set the block size for column calculations at compile time, see
    /// the comments on [`crate::inplace::coefficients2d_blocked`].
    pub fn col_blocksize<const RB: usize>(self) -> SplineInterpBuilder<2, B, C, RB, CB> {
        let SplineInterpBuilder {
            bspline,
            coeffs,
            boundary_extension,
        } = self;
        SplineInterpBuilder {
            bspline,
            coeffs,
            boundary_extension,
        }
    }
}

impl<const DIM: usize, B, C, const RB: usize> SplineInterpBuilder<DIM, B, C, RB, 1> {
    #[inline]
    /// Set the block size for row calculations at compile time, see
    /// the comments on [`crate::inplace::coefficients2d_blocked`].
    pub fn row_blocksize<const CB: usize>(self) -> SplineInterpBuilder<DIM, B, C, RB, CB> {
        let SplineInterpBuilder {
            bspline,
            coeffs,
            boundary_extension,
        } = self;
        SplineInterpBuilder {
            bspline,
            coeffs,
            boundary_extension,
        }
    }
}

impl<const DIM: usize, B, const RB: usize, const CB: usize>
    SplineInterpBuilder<DIM, B, (), RB, CB>
{
    /// Set the boundary extension method used by the interpolator. If this
    /// method is not invoked, the whole-symmetric extension is used. Obviously,
    /// this must be called before the call to [`build`](Self::build).
    pub fn extension(self, boundary_extension: InplaceBoundaryExtension) -> Self {
        Self {
            boundary_extension,
            ..self
        }
    }
}

// specialization for 2d interpolation
impl<B, const RB: usize, const CB: usize> SplineInterpBuilder<2, B, (), RB, CB> {
    /// **2D Interpolation Specialization**
    /// Try and perform the necessary precomputations (i.e. coefficient calculations)
    /// to obtain a usable interpolator which we can use to interpolate the
    /// given `image_data`.
    ///
    /// To understand the meaning of `epsilon` and its recommended values, see
    /// the documentation of [`crate::inplace::coefficients2d_blocked`].
    /// In short, `epsilon` bounds the _relative_ accuracy of the interpolation
    /// _on the original grid points_. A conservative rule of thumb is to chose
    /// `1e-12` for `f64` data and `1e-6` for `f32`.
    pub fn build<const NPOLES: usize, C>(
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
        crate::inplace::coefficients2d(
            bspline,
            RowMajor {
                data: data.as_mut(),
                width,
                height,
            },
            epsilon,
            boundary_extension,
            RowBlockSize::<RB>,
            ColBlockSize::<CB>,
        )?;
        Ok(SplineInterpolator2d {
            coeffs: RowMajor {
                data,
                width,
                height,
            },
            phantom: Default::default(),
            boundary_extension,
        })
    }
}

// specialization for 1D interpolation

impl<B, const RB: usize, const CB: usize> SplineInterpBuilder<1, B, (), RB, CB> {
    /// **1D Interpolation Specialization**. Build the 1D interpolator by
    /// performing the necessary precomputations. The whole content of `data` is
    /// assumed as 1D data for interpolation.
    ///
    /// To understand the meaning of `epsilon` and its recommended values, see
    /// the documentation of [`crate::inplace::coefficients2d_blocked`].
    /// In short, `epsilon` bounds the _relative_ accuracy of the interpolation
    /// _on the original grid points_. A conservative rule of thumb is to chose
    /// `1e-12` for `f64` data and `1e-6` for `f32`.
    pub fn build<const NPOLES: usize, C>(
        self,
        mut data: C,
        epsilon: B::Float,
    ) -> Result<SplineInterpolator1d<B, C>, Error>
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
        crate::inplace::coefficients1d_blocked(
            bspline,
            data.as_mut(),
            epsilon,
            boundary_extension,
            RowBlockSize::<RB>,
        )?;
        Ok(SplineInterpolator1d {
            coeffs: data,
            phantom: Default::default(),
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
        interp2d_single::<B, _>(
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

impl<B, C> SplineInterpolator2d<B, C>
where
    B: BSpline,
{
    /// Consume the interpolator and return the calculated
    /// B-Spline coefficients in row major fashion
    pub fn coefficients(self) -> RowMajor<C> {
        self.coeffs
    }
}

impl<B, C> SplineInterpolator1d<B, C> {
    /// Interpolate the value of a single point at the given floating point
    /// indices for the `row` and `column` of the image.
    ///
    /// **Attention**: this function indexes `(row, col)`, which
    /// is `(y,x)` in the common `x,y` coordinate system.
    pub fn interp_single<const NWEIGHTS: usize>(&self, pos: B::Float) -> Result<B::Float, Error>
    where
        C: AsRef<[B::Float]>,
        B: BSplineWeights<NWEIGHTS>,
        B::Float: Float + FromPrimitive + ToPrimitive + ConstOne + ConstZero + NumCast,
    {
        interp1d_single::<B, _>(self.coeffs.as_ref(), self.boundary_extension, pos)
    }
}

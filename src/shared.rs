//! Shared logic and types for the [`inplace`](crate::inplace) module and
//! future extended domain calculation modules.

#[cfg(test)]
mod test;

/// Type-safe way of passing the row-blocksize for calculations
#[derive(Debug)]
pub struct RowBlockSize<const N: usize>;

/// Type-safe way of passing the column-blocksize for calculations
#[derive(Debug)]
pub struct ColBlockSize<const N: usize>;

impl Default for RowBlockSize<1> {
    fn default() -> Self {
        Self {}
    }
}

impl Default for ColBlockSize<1> {
    fn default() -> Self {
        Self {}
    }
}

/// Since image coordinates can be a little confusing, this library makes it
/// explicit which logic is used to specify a location in the image.
/// We can specify coordinates by either giving their (floating point) row and
/// column index (where the index starts at `0`) _or_ we can give x,y coordinates
/// where `(x,y) = (0,0)` is the top left corner of an image as is common in
/// image processing and `x` runs along the width of the image (meaning `x` indexes
/// the _columns_).
///
/// Both ways of specifying the location below are equivalent, but one might be
/// more natural to you than the other:
///
/// ```rust
/// # use splinterpol::Coord;
/// let p1 = Coord::x(-3.7).y(10.5);
/// let p2 = Coord::row(10.5).col(-3.7);
/// assert_eq!(p1, p2);
/// ```
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Coord<T>
where
    T: Copy,
{
    pub(crate) row: T,
    pub(crate) col: T,
}

impl<T: Copy> Coord<T> {
    /// Construct a coordinate with x,y coordinate logic
    ///
    /// ```rust
    /// # use splinterpol::Coord;
    /// let coord = Coord::x(-3.7).y(10.5);
    /// ```
    #[inline]
    pub fn x(x: T) -> UnfinishedCoord<(), T> {
        UnfinishedCoord { row: (), col: x }
    }

    /// Construct a coordinate with row and column logic:
    ///
    /// ```rust
    /// # use splinterpol::Coord;
    /// let coord = Coord::row(10.5).col(-3.7);
    /// ```
    #[inline]
    pub fn row(row: T) -> UnfinishedCoord<T, ()> {
        UnfinishedCoord { row, col: () }
    }
}

/// Just an intermediate struct in the build process with no use on its
/// own.
pub struct UnfinishedCoord<U, V> {
    row: U,
    col: V,
}

// NOTE: since I only allow the sequence to be either
// `Location::row(...).col(...)` or `Location::x(...).y(...)`,
// I know that I don't have to put another x method here.
// In fact, by restricting the order of construction I save myself another
// generic and railroad the user to do the correct thing.
impl<T: Copy> UnfinishedCoord<T, ()> {
    #[inline]
    /// Floating point column index
    pub fn col(self, col: T) -> Coord<T> {
        let UnfinishedCoord { row, col: _col } = self;
        Coord { row, col }
    }
}

impl<T: Copy> UnfinishedCoord<(), T> {
    #[inline]
    /// Floating point y coordinate.
    pub fn y(self, y: T) -> Coord<T> {
        let UnfinishedCoord { row: _row, col } = self;
        Coord { col, row: y }
    }
}
use core::num::NonZero;

#[derive(Copy, Clone, PartialEq, Eq)]
/// A thin wrapper structure around row-major 2D image data stored in a container
/// `C` that might be owning (like a `Vec`) or non-owning (like a slice).
/// This wrapper alone does not enforce that the data is indeed row major nor
/// does it enforce that the number of elements in the data structure match
/// the product of `width` and `height`.
pub struct RowMajor<C> {
    /// Data structure containing the row major 2D image
    pub data: C,
    /// The width of the image (i.e. the _number of columns_,
    /// i.e. the _length of a row_).
    pub width: NonZero<usize>,
    /// The height of the image (i.e. the _number of rows_,
    /// i.e. the _length of a column_).
    pub height: NonZero<usize>,
}

impl<C> core::fmt::Debug for RowMajor<C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("RowMajor")
            .field("width", &self.width)
            .field("height", &self.height)
            .finish_non_exhaustive()
    }
}

impl<T, C> AsRef<[T]> for RowMajor<C>
where
    C: AsRef<[T]>,
{
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl<T, C> AsMut<[T]> for RowMajor<C>
where
    C: AsMut<[T]>,
{
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.data.as_mut()
    }
}

impl<C> RowMajor<C> {
    /// just a convenience implementation to get a value the given 2D index,
    /// don't rely on this for high performance computations
    #[inline]
    pub fn at<T>(&self, row: usize, col: usize) -> &T
    where
        C: AsRef<[T]>,
    {
        &self.data.as_ref()[row * self.width.get() + col]
    }
}

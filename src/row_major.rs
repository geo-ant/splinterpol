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

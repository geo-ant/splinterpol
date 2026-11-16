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
/// # use splinterpol::shared::Coord;
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
    /// # use splinterpol::shared::Coord;
    /// let coord = Coord::x(-3.7).y(10.5);
    /// ```
    #[inline]
    pub fn x(x: T) -> UnfinishedCoord<(), T> {
        UnfinishedCoord { row: (), col: x }
    }

    /// Construct a coordinate with row and column logic:
    ///
    /// ```rust
    /// # use splinterpol::shared::Coord;
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

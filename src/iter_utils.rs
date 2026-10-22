use core::num::NonZero;

mod blocked_strided_slice;
mod strided_slice;

pub use blocked_strided_slice::BlockedStridedSliceMut;
pub use strided_slice::StridedSliceMut;

#[cfg(test)]
mod test;

/// Utilitly trait for iterating over a slice with a given stride. For strides
/// of `1` use regular slices rather than this iterator.
pub trait StridedIter {
    type Item: Copy;
    /// an iterator over the elements in the sequence
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item>;
    /// a reverse iterator over the elements in the sequence
    // NOTE: adding the explicit ..._rev implementations turned out I could
    // make this more efficient than using the ExactSizeDoubleEndedIterator bound on
    // the normal iter/iter_mut methods. I could do a couple things that helped
    // the compiler generate better code. I think. Or I chased ghosts in the
    // benchmarks and assembly but here we are
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item>;
    /// iterator with mut access to the elements in the sequence
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item>;
    /// reverse iterator with mut access to the elements in the sequence
    fn strided_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item>;
    /// is the sequence empty
    fn strided_is_empty(&self) -> bool;
    /// number of elements in the sequence
    fn strided_len(&self) -> usize;
}

/// abstracts over sequence access in a strided or non-strided slice.
/// Has some extra methods to get the first, last, and second to last element
/// because the recursive filter function requires that.
///
/// This trait is implemented for regular slices, as well as for the special
/// `[StridedSliceMut]` type. For strides equal to `1` just use regular slices.
// NOTE: the `strided_` prefixes are because otherwise it's just annoying
// because the compiler will get confused with the methods on slice of the
// same name...
pub trait StridedRandomAccess: StridedIter {
    /// get a copy of the first element (if it exists)
    fn strided_first(&self) -> Option<Self::Item>;
    /// get copy of the last element (if it exists)
    fn strided_last(&self) -> Option<Self::Item>;
    /// and the second to last elements... because the filter needs that
    /// for one of the initialization passes for one of the boundary conditions
    fn strided_2nd_last(&self) -> Option<Self::Item>;
    /// write to the last element. Fails silently if there's no such element
    /// (at least in release mode)
    fn strided_write_last(&mut self, value: Self::Item);
    /// write to the first element. Fails silently if there's no such element
    /// (at least in release mode)
    fn strided_write_first(&mut self, value: Self::Item);
}
/// just a way of naming the arguments to make the `new` function less
/// confusing
pub struct StridedAccessArgs {
    pub offset: usize,
    pub stride: NonZero<usize>,
    pub count: NonZero<usize>,
}

impl<T> StridedIter for &mut [T]
where
    T: Copy,
{
    type Item = T;

    #[inline]
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.iter()
    }

    #[inline]
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item> {
        self.iter().rev()
    }

    #[inline]
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.iter_mut()
    }

    #[inline]
    fn strided_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.iter_mut().rev()
    }

    #[inline]
    fn strided_len(&self) -> usize {
        self.len()
    }

    #[inline]
    fn strided_is_empty(&self) -> bool {
        self.is_empty()
    }
}

impl<T> StridedRandomAccess for &mut [T]
where
    T: Copy,
{
    #[inline]
    fn strided_first(&self) -> Option<Self::Item> {
        self.first().copied()
    }

    #[inline]
    fn strided_last(&self) -> Option<Self::Item> {
        self.last().copied()
    }

    #[inline]
    fn strided_write_last(&mut self, value: Self::Item) {
        if let Some(last) = self.last_mut() {
            *last = value;
        }
    }

    #[inline]
    fn strided_2nd_last(&self) -> Option<Self::Item> {
        let len = self.strided_len();
        if len >= 2 { Some(self[len - 2]) } else { None }
    }

    #[inline]
    fn strided_write_first(&mut self, value: Self::Item) {
        if let Some(first) = self.first_mut() {
            *first = value;
        }
    }
}

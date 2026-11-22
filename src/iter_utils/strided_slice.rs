use super::StridedIter;
use super::StridedRandomAccess;
use crate::iter_utils::StridedAccessArgs;
use core::num::NonZero;

#[cfg(test)]
mod test;

/// Access a sequence in a slice in strided fashion. Useful for iterating over the
/// columns of an image that's stored as a contiguous slice in row major
/// layout.
pub struct StridedSliceMut<'a, T> {
    // this actually doesn't contain the full slice but starts at the first
    // element of the elements we want to access and ends at the last.
    // The stride indicates the number of elements in between and on construction
    // we enforce that the layout actually works and this is a multiple of
    // stride
    data: &'a mut [T],
    stride: NonZero<usize>,
}

impl<'a, T> StridedSliceMut<'a, T> {
    #[inline(always)]
    /// create a new instance, make sure that the dimensions are okay.
    pub fn new(data: &'a mut [T], access: StridedAccessArgs) -> Option<Self> {
        let count = access.count.get();
        let stride = access.stride.get();

        let last = (count - 1)
            .checked_mul(stride)?
            .checked_add(access.offset)?;

        if last >= data.len() {
            return None;
        }

        Some(Self {
            data: &mut data[access.offset..=last],
            stride: access.stride,
        })
    }
}

impl<'a, T> StridedIter for StridedSliceMut<'a, T>
where
    T: Copy,
{
    type Item = T;

    #[inline(always)]
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item> {
        self.data.iter().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item> {
        self.data.iter().rev().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.data.iter_mut().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.data.iter_mut().rev().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_is_empty(&self) -> bool {
        // this is actually always false because we enforce nonzero count,
        // but that's alright.
        self.data.is_empty()
    }

    #[inline(always)]
    fn strided_len(&self) -> usize {
        self.data.len().div_ceil(self.stride.get())
    }
}

impl<'a, T> StridedRandomAccess for StridedSliceMut<'a, T>
where
    T: Copy,
{
    #[inline(always)]
    fn strided_first(&self) -> Option<Self::Item> {
        self.data.first().copied()
    }

    #[inline(always)]
    fn strided_last(&self) -> Option<Self::Item> {
        self.data.last().copied()
    }

    #[inline(always)]
    fn strided_2nd_last(&self) -> Option<Self::Item> {
        self.strided_iter_rev().nth(1).copied()
    }

    #[inline(always)]
    fn strided_write_last(&mut self, value: Self::Item) {
        // NOTE: this works because the last element of the strided elements
        // is actually the last element that we want to access
        if let Some(last) = self.data.last_mut() {
            *last = value;
        }
    }

    #[inline(always)]
    fn strided_write_first(&mut self, value: Self::Item) {
        // NOTE: same as above, the first element is the actual first
        // element in the sequence
        if let Some(first) = self.data.get_mut(0) {
            *first = value;
        }
    }
}

// TODO(geo) feature gate this...
#![allow(unsafe_code)]

use core::marker::PhantomData;
use core::mem::ManuallyDrop;
use core::num::NonZero;

// TODO: bench again, but I think we can delete this...
pub(crate) struct StridedPointerAccess<T> {
    ptr: *mut T,
    stride: NonZero<usize>,
    count: NonZero<usize>,
}

#[cfg_attr(not(test), expect(dead_code))]
impl<T> StridedPointerAccess<T> {
    pub unsafe fn new(ptr: *mut T, stride: NonZero<usize>, count: NonZero<usize>) -> Option<Self> {
        if ptr.is_null() {
            return None;
        }
        if count.get().overflowing_mul(stride.get()).1 {
            return None;
        }
        Some(Self { ptr, stride, count })
    }

    #[inline(always)]
    pub fn at(&self, index: usize) -> T {
        // TODO re enable
        // if index < self.count.get() {
        unsafe { self.ptr.add(index * self.stride.get()).read() }
        // } else {
        //     panic!(
        //         "out of bounds access of index {} in [0, {})",
        //         index, self.count
        //     );
        // }
    }

    #[inline(always)]
    pub fn at_mut(&'_ mut self, index: usize) -> WriteAccessAt<'_, T> {
        // TODO re-enable
        // if index < self.count.get() {
        unsafe {
            WriteAccessAt {
                ptr: self.ptr.add(index * self.stride.get()),
                phantom: Default::default(),
            }
        }
        // } else {
        //     panic!(
        //         "out of bounds access of index {} in [0, {})",
        //         index, self.count
        //     );
        // }
    }

    #[inline(always)]
    pub fn iter(&'_ self) -> StridedPointerIter<'_, T> {
        unsafe {
            StridedPointerIter {
                cur: self.ptr,
                stride: self.stride,
                end: self.ptr.add(self.count.get() * self.stride.get()),
                phantom: Default::default(),
            }
        }
    }

    #[inline(always)]
    #[expect(dead_code)]
    pub fn iter_mut(&'_ mut self) -> StridedPointerIterMut<'_, T> {
        unsafe {
            StridedPointerIterMut {
                cur: self.ptr,
                stride: self.stride,
                end: self.ptr.add(self.count.get() * self.stride.get()),
                phantom: Default::default(),
            }
        }
    }

    #[inline(always)]
    pub fn first(&self) -> T {
        unsafe { self.ptr.read() }
    }

    #[inline(always)]
    pub fn last(&self) -> T {
        unsafe {
            self.ptr
                .add((self.count.get() - 1) * self.stride.get())
                .read()
        }
    }

    pub fn len(&self) -> usize {
        self.count.get()
    }
}

#[cfg_attr(not(test), expect(dead_code))]
pub(crate) struct WriteAccessAt<'a, T> {
    ptr: *mut T,
    phantom: PhantomData<&'a T>,
}

#[cfg_attr(not(test), expect(dead_code))]
impl<'a, T> WriteAccessAt<'a, T> {
    #[inline(always)]
    pub fn write(self, val: T) {
        let me = ManuallyDrop::new(self);
        unsafe {
            me.ptr.write(val);
        }
    }
}

impl<'a, T> Drop for WriteAccessAt<'a, T> {
    fn drop(&mut self) {
        const {
            panic!(
                "you're using this wrong. Call either read or write after requesting index access"
            )
        }
    }
}

pub(crate) struct StridedPointerIter<'a, T> {
    cur: *const T,
    stride: NonZero<usize>,
    end: *const T,
    phantom: PhantomData<&'a T>,
}

impl<'a, T: 'a> Iterator for StridedPointerIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur == self.end {
            return None;
        }
        unsafe {
            let item = self.cur.as_ref_unchecked();
            self.cur = self.cur.add(self.stride.get());
            Some(item)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }
}

impl<T> ExactSizeIterator for StridedPointerIter<'_, T> {
    fn len(&self) -> usize {
        unsafe { self.end.sub(self.cur as _) as usize / self.stride }
    }
}

impl<T> DoubleEndedIterator for StridedPointerIter<'_, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.end == self.cur {
            return None;
        }
        let item = unsafe {
            self.end = self.end.sub(self.stride.get());
            self.end.as_ref_unchecked()
        };
        Some(item)
    }
}

pub(crate) struct StridedPointerIterMut<'a, T> {
    cur: *mut T,
    stride: NonZero<usize>,
    end: *mut T,
    phantom: PhantomData<&'a T>,
}

impl<'a, T: 'a> Iterator for StridedPointerIterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur == self.end {
            return None;
        }
        unsafe {
            let item = self.cur.as_mut_unchecked();
            self.cur = self.cur.add(self.stride.get());
            Some(item)
        }
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }
}

impl<T> ExactSizeIterator for StridedPointerIterMut<'_, T> {
    fn len(&self) -> usize {
        unsafe { self.end.sub(self.cur as _) as usize / self.stride }
    }
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
pub trait StridedSequence {
    type Item: Copy;
    /// an iterator over the elements in the sequence
    fn strided_iter(&self) -> impl ExactSizeIterator<Item = &Self::Item>;
    /// a reverse iterator over the elements in the sequence
    // NOTE: adding the explicit ..._rev implementations turned out I could
    // make this more efficient than using the ExactSizeDoubleEndedIterator bound on
    // the normal iter/iter_mut methods. I could do a couple things that helped
    // the compiler generate better code. I think. Or I chased ghosts in the
    // benchmarks and assembly but here we are
    fn strided_iter_rev(&self) -> impl ExactSizeIterator<Item = &Self::Item>;
    /// iterator with mut access to the elements in the sequence
    fn strided_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut Self::Item>;
    /// reverse iterator with mut access to the elements in the sequence
    fn strided_iter_mut_rev(&mut self) -> impl ExactSizeIterator<Item = &mut Self::Item>;
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
    /// is the sequence empty
    fn strided_is_empty(&self) -> bool;
    /// number of elements in the sequence
    fn strided_len(&self) -> usize;
}

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

/// just a way of naming the arguments to make the `new` function less
/// confusing
pub struct StridedAccessArgs {
    pub offset: usize,
    pub stride: NonZero<usize>,
    pub count: NonZero<usize>,
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

impl<'a, T> StridedSequence for StridedSliceMut<'a, T>
where
    T: Copy,
{
    type Item = T;

    #[inline(always)]
    fn strided_iter(&self) -> impl ExactSizeIterator<Item = &Self::Item> {
        self.data.iter().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_rev(&self) -> impl ExactSizeIterator<Item = &Self::Item> {
        self.data.iter().rev().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut Self::Item> {
        self.data.iter_mut().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_mut_rev(&mut self) -> impl ExactSizeIterator<Item = &mut Self::Item> {
        self.data.iter_mut().rev().step_by(self.stride.get())
    }

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
        self.strided_iter_rev().skip(1).next().copied()
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

    #[inline(always)]
    fn strided_is_empty(&self) -> bool {
        // this is actually always false because we enforce nonzero count,
        // but that's alright.
        self.data.is_empty()
    }

    #[inline(always)]
    fn strided_len(&self) -> usize {
        self.data.len() / self.stride.get()
    }
}

impl<'a, T> StridedSequence for &'a mut [T]
where
    T: Copy,
{
    type Item = T;

    #[inline]
    fn strided_iter(&self) -> impl ExactSizeIterator<Item = &Self::Item> {
        self.iter()
    }

    #[inline]
    fn strided_iter_rev(&self) -> impl ExactSizeIterator<Item = &Self::Item> {
        self.iter().rev()
    }

    #[inline]
    fn strided_iter_mut(&mut self) -> impl ExactSizeIterator<Item = &mut Self::Item> {
        self.iter_mut()
    }

    #[inline]
    fn strided_iter_mut_rev(&mut self) -> impl ExactSizeIterator<Item = &mut Self::Item> {
        self.iter_mut().rev()
    }

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

    #[inline]
    fn strided_is_empty(&self) -> bool {
        self.is_empty()
    }

    #[inline]
    fn strided_len(&self) -> usize {
        self.len()
    }
}

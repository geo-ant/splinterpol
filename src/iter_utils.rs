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

pub(crate) struct WriteAccessAt<'a, T> {
    ptr: *mut T,
    phantom: PhantomData<&'a T>,
}

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

pub struct StridedMut<'a, T> {
    // start at first element of the column, end at last
    data: &'a mut [T],
    stride: NonZero<usize>,
}

pub struct StridedAccessArgs {
    pub offset: usize,
    pub stride: NonZero<usize>,
    pub count: NonZero<usize>,
}

impl<'a, T> StridedMut<'a, T> {
    #[inline(always)]
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

impl<'a, T> StridedSlice for StridedMut<'a, T>
where
    T: Copy,
{
    type Item = T;

    #[inline(always)]
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        self.data.iter().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        self.data.iter().rev().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        self.data.iter_mut().step_by(self.stride.get())
    }

    #[inline(always)]
    fn strided_iter_mut_rev(
        &mut self,
    ) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
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
        if let Some(last) = self.data.last_mut() {
            *last = value;
        }
    }

    #[inline(always)]
    fn strided_write_first(&mut self, value: Self::Item) {
        if let Some(first) = self.data.get_mut(0) {
            *first = value;
        }
    }

    #[inline(always)]
    fn strided_is_empty(&self) -> bool {
        self.data.is_empty()
    }

    #[inline(always)]
    fn strided_len(&self) -> usize {
        self.data.len() / self.stride.get()
    }
}

pub trait StridedSlice {
    type Item: Copy;
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator;
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator;
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator;
    fn strided_iter_mut_rev(&mut self)
    -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator;
    fn strided_first(&self) -> Option<Self::Item>;
    fn strided_last(&self) -> Option<Self::Item>;
    fn strided_2nd_last(&self) -> Option<Self::Item>;
    fn strided_write_last(&mut self, value: Self::Item);
    fn strided_write_first(&mut self, value: Self::Item);
    fn strided_is_empty(&self) -> bool;
    fn strided_len(&self) -> usize;
}

impl<'a, T> StridedSlice for &'a mut [T]
where
    T: Copy,
{
    type Item = T;

    #[inline]
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        (&**self).iter()
    }

    #[inline]
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        (&**self).iter().rev()
    }

    #[inline]
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        (&mut **self).iter_mut()
    }

    #[inline]
    fn strided_iter_mut_rev(
        &mut self,
    ) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        (&mut **self).iter_mut().rev()
    }

    #[inline]
    fn strided_first(&self) -> Option<Self::Item> {
        (&**self).first().copied()
    }

    #[inline]
    fn strided_last(&self) -> Option<Self::Item> {
        (&**self).last().copied()
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
        (&**self).is_empty()
    }

    #[inline]
    fn strided_len(&self) -> usize {
        (&**self).len()
    }
}

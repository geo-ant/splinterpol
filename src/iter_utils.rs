// TODO(geo) feature gate this...
#![allow(unsafe_code)]

use core::marker::PhantomData;
use core::mem::ManuallyDrop;
use core::num::NonZero;

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

pub struct ColumnAccessMut<'a, T> {
    slice: &'a mut [T],
    row_len: usize,
}

impl<'a, T> ColumnAccessMut<'a, T> {
    #[inline(always)]
    pub fn new(row_major_data: &'a mut [T], row_len: usize, num_rows: usize) -> Option<Self> {
        let (len, overflow) = row_len.overflowing_mul(num_rows);
        if overflow {
            return None;
        }
        if len != row_major_data.len() {
            return None;
        }
        Some(Self {
            slice: row_major_data,
            row_len,
        })
    }
}

pub struct ColumnIteratorMut<'a, T> {
    slice: &'a mut [T],
    stride: usize,
    count: usize,
}

// impl<'a, T> Iterator for ColumnIteratorMut<'a, T> {
//     type Item = ColRefMut<'a, T>;

//     fn next(&mut self) -> Option<Self::Item> {
//         if self.slice.len() >= self.stride {
//             let returned = Some(ColRefMut {
//                 start: self.slice,
//                 stride: self.stride,
//                 count: self.count,
//             });
//             // self.slice = &mut self.slice[self.stride..];
//             returned
//         } else {
//             None
//         }
//     }
// }

pub struct StridedMut<'a, T> {
    data: &'a mut [T],
    offset: usize,
    stride: NonZero<usize>,
    count: NonZero<usize>,
}

pub struct StridedAccessArgs {
    pub offset: usize,
    pub stride: NonZero<usize>,
    pub count: NonZero<usize>,
}

impl<'a, T> StridedMut<'a, T> {
    #[inline(always)]
    pub fn new(data: &'a mut [T], access: StridedAccessArgs) -> Option<Self> {
        if data.len() <= access.count.get() * access.stride.get() + access.offset {
            Some(Self {
                data,
                offset: access.offset,
                stride: access.stride,
                count: access.count,
            })
        } else {
            None
        }
    }
}

impl<'a, T> StridedSlice for StridedMut<'a, T>
where
    T: Copy,
{
    type Item = T;

    #[inline(always)]
    fn ss_iter(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        self.data
            .iter()
            .skip(self.offset)
            .step_by(self.stride.get())
            .take(self.count.get())
    }
    #[inline(always)]
    fn ss_iter_rev(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        let last = self.offset + (self.count.get() - 1) * self.stride.get();
        self.data[..=last]
            .iter()
            .rev()
            .step_by(self.stride.get())
            .take(self.count.get())
    }

    #[inline(always)]
    fn ss_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        self.data
            .iter_mut()
            .skip(self.offset)
            .step_by(self.stride.get())
            .take(self.count.get())
    }

    #[inline(always)]
    fn ss_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        let last = self.offset + (self.count.get() - 1) * self.stride.get();

        self.data[..=last]
            .iter_mut()
            .rev()
            .step_by(self.stride.get())
            .take(self.count.get())
    }

    #[inline(always)]
    fn ss_first(&self) -> Option<Self::Item> {
        self.data.get(self.offset).copied()
    }

    #[inline(always)]
    fn ss_last(&self) -> Option<Self::Item> {
        self.data
            .get(self.offset + self.stride.get() * (self.count.get() - 1))
            .copied()
    }

    #[inline(always)]
    fn ss_2nd_to_last(&self) -> Option<Self::Item> {
        if self.count.get() < 2 {
            return None;
        }

        self.data
            .get((self.count.get() - 2) * self.stride.get() + self.offset)
            .copied()
    }

    #[inline(always)]
    fn ss_write_last(&mut self, value: Self::Item) {
        if let Some(last) = self
            .data
            .get_mut(self.offset + (self.count.get() - 1) * self.stride.get())
        {
            *last = value;
        }
    }

    #[inline(always)]
    fn ss_write_first(&mut self, value: Self::Item) {
        if let Some(first) = self.data.get_mut(self.offset) {
            *first = value;
        }
    }

    #[inline(always)]
    fn ss_is_empty(&self) -> bool {
        // NOTE is always false
        self.count.get() == 0
    }

    #[inline(always)]
    fn ss_len(&self) -> usize {
        self.count.get()
    }
}

pub trait StridedSlice {
    type Item: Copy;
    fn ss_iter(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator;
    fn ss_iter_rev(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator;
    fn ss_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator;
    fn ss_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator;
    fn ss_first(&self) -> Option<Self::Item>;
    fn ss_last(&self) -> Option<Self::Item>;
    fn ss_2nd_to_last(&self) -> Option<Self::Item>;
    fn ss_write_last(&mut self, value: Self::Item);
    fn ss_write_first(&mut self, value: Self::Item);
    fn ss_is_empty(&self) -> bool;
    fn ss_len(&self) -> usize;
}

impl<'a, T> StridedSlice for &'a mut [T]
where
    T: Copy,
{
    type Item = T;

    #[inline]
    fn ss_iter(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        (&**self).iter()
    }

    #[inline]
    fn ss_iter_rev(&self) -> impl Iterator<Item = &Self::Item> + ExactSizeIterator {
        (&**self).iter().rev()
    }

    #[inline]
    fn ss_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        (&mut **self).iter_mut()
    }

    #[inline]
    fn ss_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item> + ExactSizeIterator {
        (&mut **self).iter_mut().rev()
    }

    #[inline]
    fn ss_first(&self) -> Option<Self::Item> {
        (&**self).first().copied()
    }

    #[inline]
    fn ss_last(&self) -> Option<Self::Item> {
        (&**self).last().copied()
    }

    #[inline]
    fn ss_write_last(&mut self, value: Self::Item) {
        if let Some(last) = self.last_mut() {
            *last = value;
        }
    }

    #[inline]
    fn ss_2nd_to_last(&self) -> Option<Self::Item> {
        let len = self.ss_len();
        if len >= 2 { Some(self[len - 2]) } else { None }
    }

    #[inline]
    fn ss_write_first(&mut self, value: Self::Item) {
        if let Some(first) = self.first_mut() {
            *first = value;
        }
    }

    #[inline]
    fn ss_is_empty(&self) -> bool {
        (&**self).is_empty()
    }

    #[inline]
    fn ss_len(&self) -> usize {
        (&**self).len()
    }
}

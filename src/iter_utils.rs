use core::marker::PhantomData;
use core::num::NonZero;
use core::ops::Index;
use core::ops::IndexMut;
use core::ptr::NonNull;

pub(crate) struct StridedSliceMut<'a, T> {
    data: &'a mut [T],
    first: usize,
    stride: NonZero<usize>,
    count: NonZero<usize>,
}

pub(crate) struct StridedSliceArgs {
    pub first: usize,
    pub stride: NonZero<usize>,
    pub count: NonZero<usize>,
}

#[expect(dead_code)]
impl<'a, T> StridedSliceMut<'a, T> {
    pub fn new(args: StridedSliceArgs, data: &'a mut [T]) -> Option<Self> {
        let StridedSliceArgs {
            first,
            stride,
            count,
        } = args;
        debug_assert!(first < data.len());
        debug_assert!((data.len() - first) / stride >= count.get());
        debug_assert!(!data.is_empty());
        if data.is_empty() {
            return None;
        }

        if data.len() <= first {
            return None;
        }

        if (data.len() - first) / stride < count.get() {
            return None;
        }

        Some(Self {
            data,
            first,
            stride,
            count,
        })
    }

    #[inline(always)]
    pub fn iter_mut(
        &'a mut self,
    ) -> impl ExactSizeIterator<Item = &'a mut T> + DoubleEndedIterator {
        self.data
            .iter_mut()
            .skip(self.first)
            .step_by(self.stride.get())
            .take(self.count.get())
    }

    #[inline(always)]
    pub fn iter(&'a self) -> impl ExactSizeIterator<Item = &'a T> + DoubleEndedIterator {
        self.data
            .iter()
            .skip(self.first)
            .step_by(self.stride.get())
            .take(self.count.get())
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        let len = (self.data.len() - self.first) / self.stride.get();
        debug_assert_eq!(self.iter().count(), self.iter().len());
        debug_assert_eq!(self.iter().len(), len);
        debug_assert_eq!(
            len,
            self.count.get(),
            "unexpected slice len: count is {}, but calculated length is {}",
            self.count,
            len
        );
        len
    }

    #[inline(always)]
    pub fn first(&self) -> Option<&'_ T> {
        self.data.first()
    }

    #[inline(always)]
    pub fn last(&self) -> Option<&'_ T> {
        self.data.last()
    }

    #[inline(always)]
    pub fn get(&self, index: usize) -> Option<&T> {
        self.data.get(self.to_index(index))
    }

    #[inline(always)]
    pub fn get_mut(&mut self, index: usize) -> Option<&mut T> {
        self.data.get_mut(self.to_index(index))
    }

    #[inline(always)]
    pub fn as_slice(&self) -> &[T] {
        self.data
    }

    #[inline(always)]
    pub fn as_slice_mut(&mut self) -> &mut [T] {
        self.data
    }

    #[inline(always)]
    fn to_index(&self, element_index: usize) -> usize {
        self.first + element_index * self.stride.get()
    }
}

impl<'a, T> Index<usize> for StridedSliceMut<'a, T> {
    type Output = T;

    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        self.data.index(self.to_index(index))
    }
}

impl<'a, T> IndexMut<usize> for StridedSliceMut<'a, T> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.data.index_mut(self.to_index(index))
    }
}

pub(crate) struct StridedPointerAccess<T> {
    ptr: *mut T,
    stride: NonZero<usize>,
    count: usize,
}

impl<T> StridedPointerAccess<T> {
    pub unsafe fn new(ptr: *mut T, stride: NonZero<usize>, count: usize) -> Option<Self> {
        if ptr.is_null() {
            return None;
        }
        Some(Self { ptr, stride, count })
    }

    pub fn at(&self, index: usize) -> T {
        if index < self.count {
            unsafe { self.ptr.add(index * self.stride.get()).read() }
        } else {
            panic!(
                "out of bounds access of index {} in [0, {})",
                index, self.count
            );
        }
    }

    pub fn at_mut(&'_ mut self, index: usize) -> WriteAccessAt<'_, T> {
        if index < self.count {
            unsafe {
                WriteAccessAt {
                    ptr: self.ptr.add(index * self.stride.get()),
                    phantom: Default::default(),
                }
            }
        } else {
            panic!(
                "out of bounds access of index {} in [0, {})",
                index, self.count
            );
        }
    }
}

pub(crate) struct WriteAccessAt<'a, T> {
    ptr: *mut T,
    phantom: PhantomData<&'a T>,
}

impl<'a, T> WriteAccessAt<'a, T> {
    pub fn write(self, val: T) {
        unsafe {
            self.ptr.write(val);
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
    cur: *mut T,
    stride: NonZero<usize>,
    end: *mut T,
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
}

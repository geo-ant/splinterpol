use core::num::NonZero;
use core::ops::Index;
use core::ops::IndexMut;

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
    ) -> impl Iterator<Item = &'_ mut T> + ExactSizeIterator + DoubleEndedIterator {
        self.data
            .iter_mut()
            .skip(self.first)
            .step_by(self.stride.get())
            .take(self.count.get())
    }

    #[inline(always)]
    pub fn iter(&'a self) -> impl Iterator<Item = &'_ T> + ExactSizeIterator + DoubleEndedIterator {
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
    fn to_index(&self, element_index: usize) -> usize {
        self.first + element_index * self.stride.get()
    }
}

impl<'a, T> Index<usize> for StridedSliceMut<'a, T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        self.data.index(self.to_index(index))
    }
}

impl<'a, T> IndexMut<usize> for StridedSliceMut<'a, T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.data.index_mut(self.to_index(index))
    }
}

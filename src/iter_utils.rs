// TODO(geo) feature gate this...
#![allow(unsafe_code)]

use core::num::NonZero;

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

// TODO: document
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
        self.data.len() / self.stride.get()
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

// TODO:doc
pub struct BlockedStridedSliceMut<'a, const BLOCKSIZE: usize, T> {
    // this actually doesn't contain the full slice but starts at the first
    // element of the elements we want to access and ends at the last.
    // The stride indicates the number of elements in between and on construction
    // we enforce that the layout actually works and this is a multiple of
    // stride
    data: &'a mut [T],
    stride: NonZero<usize>,
}

impl<'a, const BLOCKSIZE: usize, T> BlockedStridedSliceMut<'a, BLOCKSIZE, T> {
    #[inline(always)]
    /// create a new instance, make sure that the dimensions are okay.
    pub fn new(data: &'a mut [T], access: StridedAccessArgs) -> Option<Self> {
        const {
            assert!(BLOCKSIZE > 0);
        }
        let count = access.count.get();
        let stride = access.stride.get();

        // NOTE(geo): we're using a trick here to make returning blocks
        // work well with iterators. For that, we're dividing the
        // strided slice into chunks of size `stride` (not block-size)
        // and then returning the first block-size elements from the chunk.
        //
        // Example:
        // say offset = 1
        // stride = 4
        // blocksize = 2
        // count = 3
        //
        // |0 |1 |2 |3 |4 |5 |6 |7 |8 |9 |10|11|12|13|14|15...
        //    |           |           |           |  <-- stride
        //    |block|     |block|     |block|        <-- blocks
        //                                         ^
        //                                         | end = offset+count*stride
        //                                           (past the end element)
        //    |-----------------------------------|
        //     this is what we slice to do iteration
        //
        // note htis requires that the block size is <= stride,
        // which makes sense since stride is the number of columns
        // (=length of a row) in a row major storage format.
        let past_the_end = count.checked_mul(stride)?.checked_add(access.offset)?;

        if stride < BLOCKSIZE {
            return None;
        }

        if past_the_end > data.len() {
            return None;
        }

        Some(Self {
            data: &mut data[access.offset..past_the_end],
            stride: access.stride,
        })
    }
}

// NOTE(geo): I've verified that the generated assembly in the exponential
// filter function doesn't contain panic checks. Everything is inlined as
// it should be.
#[allow(clippy::unwrap_used)]
impl<'a, const BLOCKSIZE: usize, T> StridedIter for BlockedStridedSliceMut<'a, BLOCKSIZE, T>
where
    T: Copy,
{
    type Item = [T; BLOCKSIZE];

    #[inline(always)]
    fn strided_iter(&self) -> impl Iterator<Item = &Self::Item> {
        // NOTE: the trick is to chunk the data according to the stride and
        // then extract the blocks from those stride-sized chunks. This
        // will give use the multi-columns
        self.data
            .chunks_exact(self.stride.get())
            // NOTE: this works without panic checks in the tight loop of
            // the exponential filter but it's ugly due to the unwrap
            .map(move |block| block.first_chunk::<BLOCKSIZE>().unwrap())
        // .flat_map(move |block| unsafe { block.first_chunk::<BLOCKSIZE>() })
    }

    #[inline(always)]
    fn strided_iter_rev(&self) -> impl Iterator<Item = &Self::Item> {
        self.data
            .chunks_exact(self.stride.get())
            .rev()
            // TODO:check whether unwrap is optimized out
            .map(move |block| block.first_chunk::<BLOCKSIZE>().unwrap())
        // // NOTE this doesn't get optimized to the same assembly
        // .flat_map(move |block| unsafe { block.first_chunk::<BLOCKSIZE>() })
    }

    #[inline(always)]
    fn strided_iter_mut(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.data
            .chunks_exact_mut(self.stride.get())
            // TODO:check whether unwrap is optimized out
            .map(move |block| block.first_chunk_mut::<BLOCKSIZE>().unwrap())
        // // NOTE this doesn't get optimized to the same assembly
        // .flat_map(move |block| unsafe { block.first_chunk_mut::<BLOCKSIZE>() })
    }

    #[inline(always)]
    fn strided_iter_mut_rev(&mut self) -> impl Iterator<Item = &mut Self::Item> {
        self.data
            .chunks_exact_mut(self.stride.get())
            .rev()
            // TODO check whether unwrap is optimized out, otherwise use unwrap unchecked
            .map(move |block| block.first_chunk_mut::<BLOCKSIZE>().unwrap())
        // .flat_map(move |block| unsafe { block.first_chunk_mut::<BLOCKSIZE>() })
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

impl<'a, const BLOCKSIZE: usize, T> StridedRandomAccess for BlockedStridedSliceMut<'a, BLOCKSIZE, T>
where
    T: Copy,
{
    fn strided_first(&self) -> Option<Self::Item> {
        self.data.first_chunk().copied()
    }

    fn strided_last(&self) -> Option<Self::Item> {
        self.data
            .rchunks_exact(self.stride.get())
            .next()?
            .first_chunk()
            .copied()
    }

    fn strided_2nd_last(&self) -> Option<Self::Item> {
        self.data
            .rchunks_exact(self.stride.get())
            .nth(1)?
            .first_chunk()
            .copied()
    }

    fn strided_write_last(&mut self, value: Self::Item) {
        if let Some(last) = self
            .data
            .rchunks_exact_mut(self.stride.get())
            .next()
            .and_then(|slice| slice.first_chunk_mut::<BLOCKSIZE>())
        {
            *last = value
        }
    }

    fn strided_write_first(&mut self, value: Self::Item) {
        if let Some(first) = self.data.first_chunk_mut::<BLOCKSIZE>() {
            *first = value
        }
    }
}

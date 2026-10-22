use super::StridedAccessArgs;
use super::StridedIter;
use super::StridedRandomAccess;
use core::num::NonZero;

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

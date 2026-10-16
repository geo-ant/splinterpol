// NOTE(geo): this is directly taken from the `hashbrown` crate here:
// https://github.com/rust-lang/hashbrown/blob/7109e3a6388b0a6a0e843befc8876864e1d44d2e/src/util.rs
//
// That crate is dual-licensed MIT/Apache, to which this crate's license is
// compatible
//
// TODO: Replace with `core::hint::{likely, unlikely}` once they are stable.
//

use num_traits::AsPrimitive;

#[inline(always)]
#[cold]
#[expect(dead_code)]
fn cold_path() {}

#[inline(always)]
#[expect(dead_code)]
pub(crate) fn likely(b: bool) -> bool {
    if b {
        true
    } else {
        cold_path();
        false
    }
}

#[inline(always)]
#[expect(dead_code)]
pub(crate) fn unlikely(b: bool) -> bool {
    if b {
        cold_path();
        true
    } else {
        false
    }
}

/// utility type to convert arrays of castable numbers (via num_traits)
/// as arryas of other numbers. I'm writing this as a struct for more
/// convenient usage because the compiler will conveniently infer the
/// necessary generic, which wouldn't be the case when written as a function.
/// In the latter case we're forced to put `_` markers in there... at least
/// I think so.
pub(crate) struct ArrayConvert<T, const N: usize> {
    data: [T; N],
}

impl<T, const N: usize> ArrayConvert<T, N> {
    #[inline(always)]
    pub const fn new(data: [T; N]) -> Self {
        Self { data }
    }
}

impl<T, const N: usize> ArrayConvert<T, N> {
    #[inline(always)]
    pub fn as_<U>(self) -> [U; N]
    where
        T: AsPrimitive<U> + 'static + Copy,
        U: 'static + Copy,
    {
        // we could have just written this call at every callsite,
        // but the wrapper makes sure that only arrays of same size will
        // be converted which adds an extra layer of compile time safety
        core::array::from_fn(|idx| self.data[idx].as_())
    }
}

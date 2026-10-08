// NOTE(geo): this is directly taken from the `hashbrown` crate here:
// https://github.com/rust-lang/hashbrown/blob/7109e3a6388b0a6a0e843befc8876864e1d44d2e/src/util.rs
//
// That crate is dual-licensed MIT/Apache, to which this crate's license is
// compatible
//
// TODO: Replace with `core::hint::{likely, unlikely}` once they are stable.
//

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

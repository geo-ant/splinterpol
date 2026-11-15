//! Shared logic and types for the [`inplace`](crate::inplace) module and
//! future extended domain calculation modules.

/// Type-safe way of passing the row-blocksize for calculations
#[derive(Debug)]
pub struct RowBlockSize<const N: usize>;

/// Type-safe way of passing the column-blocksize for calculations
#[derive(Debug)]
pub struct ColBlockSize<const N: usize>;

impl Default for RowBlockSize<1> {
    fn default() -> Self {
        Self {}
    }
}

impl Default for ColBlockSize<1> {
    fn default() -> Self {
        Self {}
    }
}

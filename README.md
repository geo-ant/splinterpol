# splinterpol

Rust-native 1D and 2D spline interpolation up to degree 12. No-alloc and
no-std compatible, fast and accurate.

!! TODO BADGES

## Overview



## Features and Performance

### `no_std` vs `std`

By default, this crate has the `no_std` feature enabled, which makes the crate
itself `#![no_std]`. That implies that floating point calculations use the `libm`
crate. If you're in an `std`-compatible environment you _might_ consider enabling
the `std` feature of this crate to use the floating point operations of the
standard library. To enable the `std` feature, make sure to disable the default
features and enable all desired features individually:

```toml
splinterpol = {version = "0.1", default-features=false, features=["std", "fma"]}
```

The author of [`interpn`](https://docs.rs/interpn/latest/interpn/) reports
[more consistent](https://jlogan.dev/blog/2025/11/10/2025-11-10-interpn-fast-interpolation/#a-brief-aside-about-libm)
behavior of `libm` across platforms, though I've (anecdotally) observed that
the generated assembly with the standard library _looks_ better (which neither
proves or disproves that it _performs better_).

### `fma`

Following [the advice](https://jlogan.dev/blog/2025/11/10/2025-11-10-interpn-fast-interpolation/#fused-multiply-add)
of the `interpn` crate author, `splinterpol` makes fused multiply-add operations 
_opt-in_. On consumer PCs, it's recommended to enable the `fma` feature, while
on embedded systems it should be a case-by-case decision.

## Comparison to Other Crates

This is a non-exhaustive list of other libraries in the Rust ecosystem that serve
similar purposes:

* [`interpn`](https://docs.rs/interpn/latest/interpn/): an excellent crate and
  in many ways an inspiration for `splinterpol`. While `splinterpol` restricts
  itself to 1D and 2D interpolation of spline orders up to 12!!!TODO MAKE SURE THIS IS RIGHT!!!
  on regular grids, `interpn` offers n-dimensional interpolation on regular and
  rectilinear grids but restricts itself to interpolation order 3 or less. Both
  crates are no-alloc and no-std compatible and should be well suited even for
  resource constrained environments.
* [`ninterp`](https://docs.rs/ninterp/latest/ninterp/): offers n-dimensional
  interpolation up to 3rd order on rectilinear grids. Is no-std compatible, but
  isn't no-alloc and introduces `ndarray` as a dependency.
* [ndarray-interp]: TODO
* [gridr]: TODO


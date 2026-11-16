# splinterpol: Spline Interpolation

Rust-native 1D and 2D spline interpolation up to degree 12. Fast, accurate, 
no-std, and no-alloc.

!! TODO BADGES

## Quick Start

This crate exposes both high-level and low-level interfaces to its functionality.
The high-level interface for 2D image interpolation is shown below and it
should be preferred because it offers good usability without runtime penalties.

```rust
// actual image data in row-major layout.
// Anything that you can convert into // a (mutable) slice will do.
let mut data = vec![1.; 1024 * 768];

// 1. Describe the image layout
let width = core::num::NonZero::new(1024).unwrap();
let height = core::num::NonZero::new(768).unwrap();
let image_data = RowMajor {
    data,
    width,
    height,
};

// 2. Create an interpolator

// epsilon the relative accuracy of the interpolation on the
// grid points. See the documentation for guidance.
let epsilon = 1e-12;

let interp = SplineInterpolator2d::builder(BSpline3::new())
    // this is optional and defaults to whole-symmetric extension
    .extension(WholeSymmetric)
    // giving the block sizes for computation is optional and
    // defaults to `1`. Block sizes are in number of elements,
    // not bytes.
    .row_blocksize::<8>()
    .col_blocksize::<1>()
    // this is mandatory and finalizes the builder by the
    // precomputing the spline coefficients. After this, we can
    // perform interpolation
    .build(image_data, epsilon)
    .unwrap();

// 3. use the interpolator to interpolate values at arbitraty coordinates.
let value = interp.interp_single(Coord::x(1.5).y(2.4)).unwrap();
```

Interpolation in `1D` happens analogously and uses the `SplineInterpolator1d`
structure. There is also a low-level interface documented in the library for
the rare usecase where you cannot move ownership of the image data into the
interpolator. There is no reason to use the low-level interface for performance,
since the high-level interface has zero runtime-cost.

## Design Principles

* **No allocations**.
* **Safe by default**: no `unsafe` unless explicitly enabled.
* **Friendly**: zero rutime-cost high-level interfaces for all. Low-level
  interfaces for the few usecases that require it.
* **Explicit**: heavy use of builder patters to make interfaces easy to
  use and hard to misuse.

Later versions of the library will add unsafe code for extra functionality
(such as parallel coefficient computations), but will make unsafe code _opt-in_
via feature gates.

## Features and Performance

### `no_std` vs `std`

By default, this crate has the `no_std` feature enabled, which makes the crate
itself `#![no_std]`. That implies that floating point calculations use the `libm`
crate. If you're in an `std`-compatible environment you _might_ consider enabling
the `std` feature of this crate to use the floating point operations of the
standard library. To enable the `std` feature, make sure to disable the default
features and enable all desired features individually, for example:

```toml
splinterpol = {version = "0.1", default-features=false, features=["std", "fma"]}
```

The author of [`interpn`](https://docs.rs/interpn/latest/interpn/) reports
[more consistent](https://jlogan.dev/blog/2025/11/10/2025-11-10-interpn-fast-interpolation/#a-brief-aside-about-libm)
behavior of `libm` across platforms, though I've (anecdotally) observed that
the generated assembly with the standard library _looks_ better (which neither
proves nor disproves that it _performs better_).

### `fma`

Following [the advice](https://jlogan.dev/blog/2025/11/10/2025-11-10-interpn-fast-interpolation/#fused-multiply-add)
of the `interpn` crate author, `splinterpol` makes fused multiply-add operations 
_opt-in_. On consumer PCs, it's recommended to enable the `fma` feature, while
on embedded systems it should be a case-by-case decision.

## Alternatives

This is a non-exhaustive list of other libraries in the Rust ecosystem that serve
a similar purpose:

* [`interpn v0.11`](https://docs.rs/interpn/latest/interpn/): an excellent crate and
  in many ways an inspiration for `splinterpol`. While `splinterpol` restricts
  itself to 1D and 2D interpolation of spline orders up to 12 on regular grids,
  `interpn` offers n-dimensional interpolation on regular and rectilinear grids,
  but restricts itself to degree 3 or less for interpolation. Both crates are
  `no-std` and `no-alloc` compatible and should be well suited even for resource
  constrained environments.
* [`ninterp v0.11`](https://docs.rs/ninterp/latest/ninterp/): offers n-dimensional
  interpolation up to 3rd order on rectilinear grids. It's no-std compatible, but
  isn't no-alloc and introduces `ndarray` as a dependency.
* [`gridr v0.7`](https://crates.io/crates/gridr): Standalone crate which is the
  backend to a larger [python library](https://github.com/CNES/gridr). It's
  currently neither `no-std` nor `no-alloc` compatible. Uses the same underlying
  algorithm as `splinterpol`, but currently supports only _odd_ B-Spline degrees
  up to 11.
* [`ndarray-interp v0.6`](https://docs.rs/ndarray-interp/latest/ndarray_interp/):
  offers 1D spline interpolation up to degree 3 and only bilinear interpolation
  for 2D images. Introduces `ndarray` as a dependency.

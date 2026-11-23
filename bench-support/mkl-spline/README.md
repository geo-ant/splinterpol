# MKL benchmark support

* Utility crate for exposing MKLs 1D spline interpolation in my benchmarks.
* Without the `mkl` feature, the crate is empty and its build script does nothing.
  It's enabled via `splinterpol`s `__compare_mkl` feature.
* We also feature-gate the code using this crate in the `splinterpol` benchmarks
  so that when `__compare_mkl` isn't selected, the crate isn't used at all.
  It just needs to compile, which is what the `mkl` flag in this crate is for.
* The build script first asks `pkg-config` for `mkl-static-lp64-seq`. If that
  fails, set the `MKLROOT` environment variable when running the benchmarks from
  the crate above. On my machine (TM), that's:

  ```sh
  export MKLROOT=/opt/intel/oneapi/mkl/2025.2
  ```

* The build script uses `intel-mkl-tool` to find both the headers and libraries
  under the prefix found by `pkg-config`, or under `MKLROOT` as a fallback.
  It doesn't search other installations or download MKL if both fail.
  This keeps the C headers and the linked libraries from the same installation.
* Linking is static, LP64, and sequential. Don't define `MKL_ILP64` when
  compiling the C code, since that would no longer match the libraries.

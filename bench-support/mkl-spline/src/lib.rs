//! minimal crate that binds specialized C-functions that call into MKL for
//! 1D spline interpolation.
#![cfg(feature = "mkl")]

use std::{marker::PhantomData, ptr::null_mut};

use crate::bindings::MklInterpolator_1d_f64;

#[allow(
    dead_code,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals
)]
mod bindings {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Error(&'static str);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Error: {}", self.0)
    }
}

impl std::error::Error for Error {}

/// Prepared MKL task for cubic interpolation with whole sym boundary conditions
pub struct MklInterpF64<'a> {
    handle: bindings::MklInterpolator_1d_f64,
    // internal coefficient storage, allocated once and never resized
    // I use to have Box<[f64]> here which sounded better than vec, but see
    // https://doc.rust-lang.org/std/boxed/index.html#considerations-for-unsafe-code
    // and https://doc.rust-lang.org/std/vec/struct.Vec.html#method.as_mut_ptr
    coeffs: &'a mut [f64],
    // MKL keeps a pointer into this slice so we'll keep it borrowed
    _data: PhantomData<&'a [f64]>,
}

impl<'a> MklInterpF64<'a> {
    pub fn new(data: &'a [f64], coeffs_out: &'a mut [f64]) -> Result<Self, Error> {
        if data.len() < 4 {
            return Err(Error("at least four samples are required"));
        }
        let ndata = data
            .len()
            .try_into()
            .map_err(|_| Error("sample count exceeds MKL limits"))?;
        if coeffs_out.len() != 4 * (data.len() - 1) {
            return Err(Error("wrong len for coefficient output array"));
        }
        let mut me = Self {
            handle: MklInterpolator_1d_f64 {
                task: null_mut(),
                range: null_mut(),
                bc: null_mut(),
            },
            coeffs: coeffs_out,
            _data: PhantomData,
        };

        if unsafe {
            bindings::MklInterpolator_1d_f64_init(
                &mut me.handle,
                data.as_ptr(),
                ndata,
                me.coeffs.as_mut_ptr(),
                me.coeffs.len(),
            )
        } {
            Ok(me)
        } else {
            Err(Error("creating mkl interpolator"))
        }
    }

    pub fn calculate_coefficients(&mut self) -> Result<(), Error> {
        if unsafe { bindings::MklInterpolator_1d_f64_calculate_coefficients(self.handle) } {
            Ok(())
        } else {
            Err(Error("calculating coeffs"))
        }
    }

    pub fn coefficients(&self) -> &[f64] {
        &self.coeffs
    }

    pub fn interp_multi_generic(&self, sites: &[f64], out: &mut [f64]) -> Result<(), Error> {
        if unsafe {
            bindings::MklInterpolator_1d_f64_interp_multi_generic(
                self.handle,
                sites.as_ptr(),
                sites.len().try_into().unwrap(),
                out.as_mut_ptr(),
            )
        } {
            Ok(())
        } else {
            Err(Error("calculating coeffs"))
        }
    }
}

impl Drop for MklInterpF64<'_> {
    fn drop(&mut self) {
        unsafe { bindings::MklInterpolator_1d_f64_free(&mut self.handle) };
    }
}

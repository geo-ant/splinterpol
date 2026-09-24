#[cfg(test)]
mod test;

pub trait FloatVector: Copy {
    type Scalar: num_traits::Float + num_traits::ConstOne + num_traits::ConstZero;

    /// BLAS inspired syntax a*x+y where x = self,y = addend are vectors and
    /// a = factor is the scalar factor. Translates to fused multiply add (hopefully)
    /// when FMA is enabled.
    fn axpy(self, factor: Self::Scalar, addend: Self) -> Self;
    /// Also blas inspired naming, where this is just scaling the vector
    /// value with a scalar
    fn scal(self, factor: Self::Scalar) -> Self;
    /// vector subtraction: `self`-`subtrahend`
    fn sub(self, subtrahend: Self) -> Self;
}

impl FloatVector for f64 {
    type Scalar = Self;

    #[inline(always)]
    fn axpy(self, factor: Self::Scalar, addend: Self) -> Self {
        // NOTE(geo) when inspecting the assembly it looks like using
        // the std::mul_add emits better instructions compared to the
        // libm-based Float::mul_add implementation. However, to be
        // honest I haven't seen that change the actual performance,
        // so this might be overengineering. But it's cheap so what
        // the heck...
        #[cfg(all(feature = "fma", feature = "std"))]
        {
            f64::mul_add(self, factor, addend)
        }
        #[cfg(all(feature = "fma", not(feature = "std")))]
        {
            num_traits::Float::mul_add(self, factor, addend)
        }
        #[cfg(not(feature = "fma"))]
        {
            self * factor + addend
        }
    }

    #[inline(always)]
    fn scal(self, factor: Self::Scalar) -> Self {
        self * factor
    }

    #[inline(always)]
    fn sub(self, subtrahend: Self) -> Self {
        self - subtrahend
    }
}

impl FloatVector for f32 {
    type Scalar = Self;

    #[inline(always)]
    fn axpy(self, factor: Self::Scalar, addend: Self) -> Self {
        // NOTE(geo) when inspecting the assembly it looks like using
        // the std::mul_add emits better instructions compared to the
        // libm-based Float::mul_add implementation. However, to be
        // honest I haven't seen that change the actual performance,
        // so this might be overengineering. But it's cheap so what
        // the heck...
        #[cfg(all(feature = "fma", feature = "std"))]
        {
            f32::mul_add(self, factor, addend)
        }
        #[cfg(all(feature = "fma", not(feature = "std")))]
        {
            num_traits::Float::mul_add(self, factor, addend)
        }
        #[cfg(not(feature = "fma"))]
        {
            self * factor + addend
        }
    }

    #[inline(always)]
    fn scal(self, factor: Self::Scalar) -> Self {
        self * factor
    }

    #[inline(always)]
    fn sub(self, subtrahend: Self) -> Self {
        self - subtrahend
    }
}

impl<const N: usize> FloatVector for [f32; N] {
    type Scalar = f32;

    fn axpy(mut self, factor: Self::Scalar, addend: Self) -> Self {
        for j in 0..N {
            // self[j] = factor * self[j] + addend[j];
            self[j] = self[j].axpy(factor, addend[j]);
        }
        self
    }

    fn scal(mut self, factor: Self::Scalar) -> Self {
        for j in 0..N {
            self[j] *= factor;
        }
        self
    }

    fn sub(mut self, subtrahend: Self) -> Self {
        for j in 0..N {
            self[j] -= subtrahend[j];
        }
        self
    }
}

impl<const N: usize> FloatVector for [f64; N] {
    type Scalar = f64;

    fn axpy(mut self, factor: Self::Scalar, addend: Self) -> Self {
        for j in 0..N {
            // self[j] = factor * self[j] + addend[j];
            self[j] = self[j].axpy(factor, addend[j]);
        }
        self
    }

    fn scal(mut self, factor: Self::Scalar) -> Self {
        for j in 0..N {
            self[j] *= factor;
        }
        self
    }

    fn sub(mut self, subtrahend: Self) -> Self {
        for j in 0..N {
            self[j] -= subtrahend[j];
        }
        self
    }
}

/// Trait exposing a fused mul-add or a fallback if the instruction doesn't exist
/// This idea is taken from the `interpN` crate.
pub trait FusedMuladd {
    fn fmuladd(self, factor: Self, addend: Self) -> Self;
}

impl FusedMuladd for f64 {
    #[inline(always)]
    /// try and make this a fused mul-add operation `self*factor+addend`
    /// and fall back on non-fma when the FMA feature is not enabled.
    fn fmuladd(self, factor: Self, addend: Self) -> Self {
        // NOTE(geo) when inspecting the assembly it looks like using
        // the std::mul_add emits better instructions compared to the
        // libm-based Float::mul_add implementation. However, to be
        // honest I haven't seen that change the actual performance,
        // so this might be overengineering. But it's cheap so what
        // the heck...
        #[cfg(all(feature = "fma", feature = "std"))]
        {
            f64::mul_add(self, factor, addend)
        }
        #[cfg(all(feature = "fma", not(feature = "std")))]
        {
            num_traits::Float::mul_add(self, factor, addend)
        }
        #[cfg(not(feature = "fma"))]
        {
            self * factor + addend
        }
    }
}

impl FusedMuladd for f32 {
    #[inline(always)]
    /// try and make this a fused mul-add operation `self*factor+addend`
    /// and fall back on non-fma when the FMA feature is not enabled.
    fn fmuladd(self, factor: Self, addend: Self) -> Self {
        // NOTE(geo) when inspecting the assembly it looks like using
        // the std::mul_add emits better instructions compared to the
        // libm-based Float::mul_add implementation. However, to be
        // honest I haven't seen that change the actual performance,
        // so this might be overengineering. But it's cheap so what
        // the heck...
        #[cfg(all(feature = "fma", feature = "std"))]
        {
            f32::mul_add(self, factor, addend)
        }
        #[cfg(all(feature = "fma", not(feature = "std")))]
        {
            num_traits::Float::mul_add(self, factor, addend)
        }
        #[cfg(not(feature = "fma"))]
        {
            self * factor + addend
        }
    }
}

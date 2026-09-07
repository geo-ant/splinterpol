#[cfg(test)]
mod test;

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

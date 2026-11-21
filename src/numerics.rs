#[cfg(test)]
mod test;

pub trait FloatVector: Copy + IntoIterator<Item = Self::Scalar>
where
    Self::Scalar: core::iter::Sum,
{
    type Scalar: num_traits::Float;

    /// BLAS inspired syntax a*x+y where x = self,y = addend are vectors and
    /// a = factor is the scalar factor. Translates to fused multiply add (hopefully)
    /// when FMA is enabled.
    fn axpy(self, factor: Self::Scalar, addend: Self) -> Self;
    /// Also blas inspired naming, where this is just scaling the vector
    /// value with a scalar
    fn scal(self, factor: Self::Scalar) -> Self;
    /// vector subtraction: `self`-`subtrahend`
    fn sub(self, subtrahend: Self) -> Self;
    /// calculate the dot product
    fn dot(self, other: Self) -> Self::Scalar;
}

impl<const N: usize, T> FloatVector for [T; N]
where
    T: num_traits::Float + num_traits::ConstZero + core::iter::Sum,
{
    type Scalar = T;

    #[inline(always)]
    fn axpy(mut self, factor: Self::Scalar, addend: Self) -> Self {
        for j in 0..N {
            self[j] = self[j].fmuladd(factor, addend[j]);
        }
        self
    }

    #[inline(always)]
    fn scal(mut self, factor: Self::Scalar) -> Self {
        for j in 0..N {
            self[j] = self[j] * factor;
        }
        self
    }

    #[inline(always)]
    fn sub(mut self, subtrahend: Self) -> Self {
        for j in 0..N {
            self[j] = self[j] - subtrahend[j];
        }
        self
    }

    #[inline(always)]
    fn dot(self, other: Self) -> Self::Scalar {
        const { assert!(N > 0) };
        // TODO PERF: see if there's a better way to group the terms to use
        // muladds like so:
        // v^T w = v[0] w[0] + v[1] w[1] + v[2] w[2] + v[3] w[3] + ...
        // could be written as
        // v[0].fmuladd(w[0], v[1] w[1]) + v[2].fmuladd(w[2], v[3] w[3]) + ...
        // This is for example how interpn does it for dot4:
        // see https://github.com/jlogan03/interpn/blob/main/src/interp_math.rs
        // but some care has to be taken whether the terms are even or odd

        // NOTE(perf): in my tests, special casing this loop for N==4 (cubic
        // bspline weights) made the benchmarks for intepolation slightly
        // but noticeably slower than just leaving the loop here...
        // we special case the N==4 path here (where the compiler should eliminate
        // // the branching) because N==4 is important for cubic B-Splines
        // if N == 4 {
        //     let first = self[1].fmuladd(other[1], self[0] * other[0]);
        //     let second = self[3].fmuladd(other[3], self[2] * other[2]);
        //     first + second
        // } else {
        //   TODO PERF: this is NOT OPTIMAL DESPITE the muladd because we
        //   still have a dependency from one element to the next
        // loop based (same as zip and fold)
        let mut sum = T::ZERO;
        for j in 0..N {
            sum = self[j].fmuladd(other[j], sum);
        }
        sum
        // }

        // // NOTE(perf): this is indeed slower when FMA is enabled
        // self.into_iter()
        //     .zip(other.into_iter())
        //     .map(|(v, w)| v * w)
        //     .sum()
    }
}

/// Trait exposing a fused mul-add or a fallback if the instruction doesn't exist
/// This idea is taken from the `interpN` crate.
pub trait FusedMuladd {
    fn fmuladd(self, factor: Self, addend: Self) -> Self;
}

impl<T> FusedMuladd for T
where
    T: num_traits::Float,
{
    #[inline(always)]
    /// try and make this a fused mul-add operation `self*factor+addend`
    /// and fall back on non-fma when the FMA feature is not enabled.
    fn fmuladd(self, factor: Self, addend: Self) -> Self {
        // the idea of feature-gating the fused-muladd is from the
        // interpN crate, see:
        // https://jlogan.dev/blog/2025/11/10/2025-11-10-interpn-fast-interpolation/

        // NOTE(geo) when inspecting the assembly it looks like using
        // the std dependency for the num_traits crate emits better
        // assembly for mul-adds. However, to be honest I haven't seen
        // that change the actual performance. But still we have the ability
        // to change the dependency of the num_traits crate through the
        // std or no-std features.
        #[cfg(feature = "fma")]
        {
            num_traits::Float::mul_add(self, factor, addend)
        }
        #[cfg(not(feature = "fma"))]
        {
            self * factor + addend
        }
    }
}

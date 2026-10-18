use std::path::Path;

/// compare a f32 and a f64 array for equality in a way that can
/// be evaluated at compile time.
pub const fn array_f32_f64_eq<const N: usize>(lhs: [f32; N], rhs: [f64; N]) -> bool {
    // can't use for loops in const functions... cause of reasons I guess.
    let mut j = 0;
    loop {
        if j >= N {
            break;
        }
        if lhs[j] != rhs[j] as f32 {
            return false;
        }
        j += 1;
    }
    true
}

#[macro_export]
// utility macro for parametrizing a test with different instantiations
// of a single const-generic parameter.
macro_rules! parametrized_test {
    (
        #[test]
        #[parametrize($($val:literal),+ $(,)?)]
        fn $fn:ident<const $generic:ident : $type:ty>() $block:block
    ) => {

        fn $fn<const $generic: $type>() $block

        paste::paste! {
            $(
                #[test]
                fn [<$fn _ $generic:lower _ $val>]() {
                    $fn::<$val>();
                }
            )+
        }
    };
}

/// read little endian f64 values from a file to a Vec
pub fn read_f64le<P: AsRef<Path>>(path: P) -> Vec<f64> {
    let path = path.as_ref().to_owned();
    let bytes = std::fs::read(&path).expect(&format!("couldn't read file '{}'", path.display()));
    assert!(
        bytes.len().is_multiple_of(size_of::<f64>()),
        "file size not a multiple of f64 size"
    );
    bytes
        .chunks_exact(size_of::<f64>())
        .map(|bytes| f64::from_le_bytes(bytes.try_into().unwrap()))
        .collect()
}

/// read little endian f32 values to a vec but convert them to f64
pub fn read_f32le_asf64<P: AsRef<Path>>(path: P) -> Vec<f64> {
    let path = path.as_ref().to_owned();
    let bytes = std::fs::read(&path).expect(&format!("couldn't read file '{}'", path.display()));
    assert!(
        bytes.len().is_multiple_of(size_of::<f32>()),
        "file size not a multiple of f64 size"
    );
    bytes
        .chunks_exact(size_of::<f32>())
        .map(|bytes| f32::from_le_bytes(bytes.try_into().unwrap()) as f64)
        .collect()
}

use std::path::Path;

/// utility macro for parametrizing a test with different instantiations
/// of a single const-generic parameter.
#[macro_export]
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

/// utility function to calculate the "cross product iterator"
pub fn iter_cross<T, const N: usize, const M: usize>(
    first: [T; N],
    second: [T; M],
) -> impl Iterator<Item = (T, T)>
where
    T: Copy + 'static,
{
    const { assert!(N > 0 && M > 0) };
    first
        .into_iter()
        .flat_map(move |first| second.into_iter().map(move |second| (first, second)))
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

#[test]
fn test_iter_cross() {
    let collected: Vec<_> = iter_cross([1, 2], [3]).collect();
    assert_eq!(collected, [(1, 3), (2, 3)]);

    let collected: Vec<_> = iter_cross([1], [2, 3]).collect();
    assert_eq!(collected, [(1, 2), (1, 3)]);

    let collected: Vec<_> = iter_cross([1, 2], [3, 4]).collect();
    assert_eq!(collected, [(1, 3), (1, 4), (2, 3), (2, 4)]);

    let collected: Vec<_> = iter_cross([1, 2, 10], [3, 4]).collect();
    assert_eq!(
        collected,
        [(1, 3), (1, 4), (2, 3), (2, 4), (10, 3), (10, 4)]
    );
}

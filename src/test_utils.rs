use std::path::Path;

#[macro_export]
macro_rules! parametrized_test {
    // instantiate a test taking one const generic arguments
    // with many values of that argument
    (
        #[test]
        #[parametrize($($val:literal),+ $(,)?)]
        fn $fnname:ident<const $generic:ident : $type:ty>() $block:block
    ) => {

        // #[track_caller]
        fn $fnname<const $generic: $type>() $block

        paste::paste! {
            $(
                #[test]
                fn [<$fnname _ $generic:lower _ $val>]() {
                    $fnname::<$val>();
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

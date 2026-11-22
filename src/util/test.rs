use super::ArrayConvert;

#[test]
fn array_convert() {
    assert_eq!([1., 2., 3., 4.], ArrayConvert::new([1, 2, 3, 4]).as_());
    assert_eq!([-1., -2., 3., 4.], ArrayConvert::new([-1, -2, 3, 4]).as_());
}

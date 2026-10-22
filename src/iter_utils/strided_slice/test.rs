use crate::iter_utils::{StridedAccessArgs, StridedIter, StridedRandomAccess, StridedSliceMut};

#[test]
fn empty_strided_slice_mut_cannot_be_constructed() {
    assert!(
        StridedSliceMut::<f32>::new(
            &mut [],
            StridedAccessArgs {
                offset: 0,
                stride: 1.try_into().unwrap(),
                count: 1.try_into().unwrap()
            }
        )
        .is_none()
    );
    assert!(
        StridedSliceMut::<f32>::new(
            &mut [],
            StridedAccessArgs {
                offset: 1,
                stride: 1.try_into().unwrap(),
                count: 2.try_into().unwrap()
            }
        )
        .is_none()
    );
}

#[test]
// just to make sure this is correct, even though we don't use stride == 1
fn strided_iter_for_strided_slice_behaves_like_slice_for_unit_stride() {
    fn gen_unit_strided<'a, T>(data: &'a mut [T]) -> StridedSliceMut<'a, T> {
        let len = data.len();
        StridedSliceMut::new(
            data,
            StridedAccessArgs {
                offset: 0,
                stride: 1.try_into().unwrap(),
                count: len.try_into().unwrap(),
            },
        )
        .unwrap()
    }

    let mut original = [1isize, 2, 3, 4, 5];
    let len = original.len();

    let mut unit_strided = gen_unit_strided(&mut original);
    assert_eq!(unit_strided.strided_len(), len);
    assert!(!unit_strided.strided_is_empty());

    let actual: Vec<_> = unit_strided.strided_iter().copied().collect();
    assert_eq!([1, 2, 3, 4, 5], actual.as_slice());

    let actual: Vec<_> = unit_strided.strided_iter_rev().copied().collect();
    assert_eq!([5, 4, 3, 2, 1], actual.as_slice());

    unit_strided
        .strided_iter_mut()
        .enumerate()
        .for_each(|(idx, val)| {
            *val += idx as isize + 1;
        });
    drop(unit_strided);

    assert_eq!(original, [2, 4, 6, 8, 10]);

    let mut unit_strided = gen_unit_strided(&mut original);
    unit_strided
        .strided_iter_mut_rev()
        .enumerate()
        .for_each(|(idx, val)| *val -= idx as isize + 1);
    assert_eq!(original, [-3, 0, 3, 6, 9]);
}

#[test]
fn strided_random_access_for_strided_slice_mut_with_unit_stride_behaves_like_slice() {
    fn gen_unit_strided<'a, T>(data: &'a mut [T]) -> StridedSliceMut<'a, T> {
        let len = data.len();
        StridedSliceMut::new(
            data,
            StridedAccessArgs {
                offset: 0,
                stride: 1.try_into().unwrap(),
                count: len.try_into().unwrap(),
            },
        )
        .unwrap()
    }

    let mut original = [1isize, 2, 3, 4, 5];
    let mut unit_strided = gen_unit_strided(&mut original);

    assert_eq!(unit_strided.strided_first(), Some(1));
    assert_eq!(unit_strided.strided_last(), Some(5));
    assert_eq!(unit_strided.strided_2nd_last(), Some(4));
    unit_strided.strided_write_first(123);
    drop(unit_strided);
    assert_eq!(original, [123, 2, 3, 4, 5]);
    let mut unit_strided = gen_unit_strided(&mut original);
    unit_strided.strided_write_last(-7);
    assert_eq!(original, [123, 2, 3, 4, -7]);
}

#[test]
fn strided_slice_with_invalid_arguments() {
    let mut data = [1, 2, 3, 4];

    assert!(
        StridedSliceMut::new(
            &mut data,
            StridedAccessArgs {
                offset: 4,
                stride: 1.try_into().unwrap(),
                count: 1.try_into().unwrap()
            }
        )
        .is_none()
    );

    assert!(
        StridedSliceMut::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 3.try_into().unwrap(),
                count: 3.try_into().unwrap()
            }
        )
        .is_none()
    );

    assert!(
        StridedSliceMut::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 5.try_into().unwrap(),
                count: 2.try_into().unwrap()
            }
        )
        .is_none()
    );
}

#[test]
fn strided_slice_with_valid_arguments_must_pass() {
    let mut data = [1, 2, 3, 4, 5];

    assert!(
        StridedSliceMut::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 1.try_into().unwrap(),
                count: 2.try_into().unwrap()
            }
        )
        .is_some()
    );

    // looks weird but we've only requested one element, the element
    // at the offset!
    assert!(
        StridedSliceMut::new(
            &mut data,
            StridedAccessArgs {
                offset: 3,
                stride: 10.try_into().unwrap(),
                count: 1.try_into().unwrap()
            }
        )
        .is_some()
    );

    assert!(
        StridedSliceMut::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 2.try_into().unwrap(),
                count: 3.try_into().unwrap()
            }
        )
        .is_some()
    );
}

#[test]
fn strided_iter_for_strided_slice_mut() {
    let mut original = [
        0isize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    ];

    let mut strided = StridedSliceMut::new(
        &mut original,
        StridedAccessArgs {
            offset: 2,
            stride: 3.try_into().unwrap(),
            count: 5.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_len(), 5);
    assert!(!strided.strided_is_empty());

    let actual: Vec<_> = strided.strided_iter().copied().collect();
    assert_eq!(actual, [2, 5, 8, 11, 14]);

    let actual: Vec<_> = strided.strided_iter_rev().copied().collect();
    assert_eq!(actual, [14, 11, 8, 5, 2]);

    strided
        .strided_iter_mut()
        .enumerate()
        .for_each(|(idx, val)| *val += idx as isize + 1);

    let actual: Vec<_> = strided.strided_iter().copied().collect();
    assert_eq!(actual, [3, 7, 11, 15, 19]);

    strided
        .strided_iter_mut_rev()
        .enumerate()
        .for_each(|(idx, val)| *val += idx as isize);

    let actual: Vec<_> = strided.strided_iter().copied().collect();
    assert_eq!(actual, [7, 10, 13, 16, 19]);

    drop(strided);

    assert_eq!(
        original,
        [
            0, 1, 7, 3, 4, 10, 6, 7, 13, 9, 10, 16, 12, 13, 19, 15, 16, 17
        ]
    );
}

#[test]
fn strided_random_access_for_strided_slice_mut() {
    let mut original = [0isize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

    let mut strided = StridedSliceMut::new(
        &mut original,
        StridedAccessArgs {
            offset: 1,
            stride: 3.try_into().unwrap(),
            count: 4.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_first(), Some(1));
    assert_eq!(strided.strided_last(), Some(10));
    assert_eq!(strided.strided_2nd_last(), Some(7));

    strided.strided_write_first(99);
    strided.strided_write_last(-10);

    assert_eq!(strided.strided_first(), Some(99));
    assert_eq!(strided.strided_last(), Some(-10));
    assert_eq!(strided.strided_2nd_last(), Some(7));

    drop(strided);

    // another test, single element stride larger than sequence,
    // but it's okay because it's only one element
    assert_eq!(original, [0, 99, 2, 3, 4, 5, 6, 7, 8, 9, -10, 11]);

    let mut strided = StridedSliceMut::new(
        &mut original,
        StridedAccessArgs {
            offset: 8,
            stride: 100.try_into().unwrap(),
            count: 1.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_first(), Some(8));
    assert_eq!(strided.strided_last(), Some(8));
    assert_eq!(strided.strided_2nd_last(), None);

    strided.strided_write_first(123);
    assert_eq!(strided.strided_last(), Some(123));

    strided.strided_write_last(456);
    assert_eq!(strided.strided_first(), Some(456));

    drop(strided);

    assert_eq!(original, [0, 99, 2, 3, 4, 5, 6, 7, 456, 9, -10, 11]);
}

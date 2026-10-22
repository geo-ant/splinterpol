use crate::iter_utils::{
    BlockedStridedSliceMut, StridedAccessArgs, StridedIter, StridedRandomAccess,
};

#[test]
fn empty_strided_slice_mut_cannot_be_constructed() {
    assert!(
        BlockedStridedSliceMut::<1, f32>::new(
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
        BlockedStridedSliceMut::<1, f32>::new(
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
fn strided_iter_for_strided_slice_for_unit_stride_and_blocksize() {
    fn gen_unit_strided<'a, T>(data: &'a mut [T]) -> BlockedStridedSliceMut<'a, 1, T> {
        let len = data.len();
        BlockedStridedSliceMut::new(
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
    assert_eq!([[1], [2], [3], [4], [5]], actual.as_slice());

    let actual: Vec<_> = unit_strided.strided_iter_rev().copied().collect();
    assert_eq!([[5], [4], [3], [2], [1]], actual.as_slice());

    unit_strided
        .strided_iter_mut()
        .enumerate()
        .for_each(|(idx, val)| {
            val[0] += idx as isize + 1;
        });
    drop(unit_strided);

    assert_eq!(original, [2, 4, 6, 8, 10]);

    let mut unit_strided = gen_unit_strided(&mut original);
    unit_strided
        .strided_iter_mut_rev()
        .enumerate()
        .for_each(|(idx, val)| val[0] -= idx as isize + 1);
    drop(unit_strided);

    assert_eq!(original, [-3, 0, 3, 6, 9]);
}

#[test]
fn strided_random_access_for_strided_slice_mut_with_unit_stride_behaves_like_slice() {
    fn gen_unit_strided<'a, T>(data: &'a mut [T]) -> BlockedStridedSliceMut<'a, 1, T> {
        let len = data.len();
        BlockedStridedSliceMut::new(
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

    assert_eq!(unit_strided.strided_first(), Some([1]));
    assert_eq!(unit_strided.strided_last(), Some([5]));
    assert_eq!(unit_strided.strided_2nd_last(), Some([4]));

    unit_strided.strided_write_first([123]);
    drop(unit_strided);

    assert_eq!(original, [123, 2, 3, 4, 5]);

    let mut unit_strided = gen_unit_strided(&mut original);
    unit_strided.strided_write_last([-7]);
    drop(unit_strided);

    assert_eq!(original, [123, 2, 3, 4, -7]);
}

#[test]
fn strided_slice_with_invalid_arguments() {
    let mut data = [1, 2, 3, 4];

    assert!(
        BlockedStridedSliceMut::<1, _>::new(
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
        BlockedStridedSliceMut::<2, _>::new(
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
        BlockedStridedSliceMut::<2, _>::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 5.try_into().unwrap(),
                count: 2.try_into().unwrap()
            }
        )
        .is_none()
    );

    assert!(
        BlockedStridedSliceMut::<3, _>::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 2.try_into().unwrap(),
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
        BlockedStridedSliceMut::<1, _>::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 1.try_into().unwrap(),
                count: 2.try_into().unwrap()
            }
        )
        .is_some()
    );

    // looks weird but we've only requested one block, the block
    // at the offset!
    assert!(
        BlockedStridedSliceMut::<2, _>::new(
            &mut data,
            StridedAccessArgs {
                offset: 3,
                stride: 2.try_into().unwrap(),
                count: 1.try_into().unwrap()
            }
        )
        .is_some()
    );

    assert!(
        BlockedStridedSliceMut::<2, _>::new(
            &mut data,
            StridedAccessArgs {
                offset: 0,
                stride: 2.try_into().unwrap(),
                count: 2.try_into().unwrap()
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

    let access = StridedAccessArgs {
        offset: 2,
        stride: 3.try_into().unwrap(),
        count: 5.try_into().unwrap(),
    };

    let mut strided = BlockedStridedSliceMut::<2, _>::new(&mut original, access).unwrap();

    assert_eq!(strided.strided_len(), 5);
    assert!(!strided.strided_is_empty());

    let actual: Vec<_> = strided.strided_iter().copied().collect();
    assert_eq!(actual, [[2, 3], [5, 6], [8, 9], [11, 12], [14, 15]]);

    let actual: Vec<_> = strided.strided_iter_rev().copied().collect();
    assert_eq!(actual, [[14, 15], [11, 12], [8, 9], [5, 6], [2, 3]]);

    strided
        .strided_iter_mut()
        .enumerate()
        .for_each(|(idx, block)| {
            for val in block {
                *val += idx as isize + 1;
            }
        });

    drop(strided);

    assert_eq!(
        original,
        [
            0, 1, 3, 4, 4, 7, 8, 7, 11, 12, 10, 15, 16, 13, 19, 20, 16, 17
        ]
    );

    let mut strided = BlockedStridedSliceMut::<2, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 2,
            stride: 3.try_into().unwrap(),
            count: 5.try_into().unwrap(),
        },
    )
    .unwrap();

    strided
        .strided_iter_mut_rev()
        .enumerate()
        .for_each(|(idx, block)| {
            for val in block {
                *val += idx as isize;
            }
        });

    drop(strided);

    assert_eq!(
        original,
        [
            0, 1, 7, 8, 4, 10, 11, 7, 13, 14, 10, 16, 17, 13, 19, 20, 16, 17
        ]
    );
}

#[test]
fn strided_random_access_for_strided_slice_mut() {
    let mut original = [0isize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

    let mut strided = BlockedStridedSliceMut::<2, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 1,
            stride: 3.try_into().unwrap(),
            count: 4.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_first(), Some([1, 2]));
    assert_eq!(strided.strided_last(), Some([10, 11]));
    assert_eq!(strided.strided_2nd_last(), Some([7, 8]));

    strided.strided_write_first([99, 98]);
    drop(strided);

    assert_eq!(original, [0, 99, 98, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]);

    let mut strided = BlockedStridedSliceMut::<2, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 1,
            stride: 3.try_into().unwrap(),
            count: 4.try_into().unwrap(),
        },
    )
    .unwrap();

    strided.strided_write_last([-10, -11]);
    drop(strided);

    assert_eq!(original, [0, 99, 98, 3, 4, 5, 6, 7, 8, 9, -10, -11, 12, 13]);

    // another test, single block with stride larger than block size
    let mut strided = BlockedStridedSliceMut::<2, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 8,
            stride: 6.try_into().unwrap(),
            count: 1.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_first(), Some([8, 9]));
    assert_eq!(strided.strided_last(), Some([8, 9]));
    assert_eq!(strided.strided_2nd_last(), None);

    strided.strided_write_first([123, 124]);
    drop(strided);

    assert_eq!(
        original,
        [0, 99, 98, 3, 4, 5, 6, 7, 123, 124, -10, -11, 12, 13]
    );

    let mut strided = BlockedStridedSliceMut::<2, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 8,
            stride: 6.try_into().unwrap(),
            count: 1.try_into().unwrap(),
        },
    )
    .unwrap();

    strided.strided_write_last([456, 457]);
    drop(strided);

    assert_eq!(
        original,
        [0, 99, 98, 3, 4, 5, 6, 7, 456, 457, -10, -11, 12, 13]
    );
}

#[test]
fn strided_iter_with_blocksize_equal_to_stride() {
    let mut original = [0isize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

    let mut strided = BlockedStridedSliceMut::<3, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 1,
            stride: 3.try_into().unwrap(),
            count: 3.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_len(), 3);
    assert!(!strided.strided_is_empty());

    let actual: Vec<_> = strided.strided_iter().copied().collect();
    assert_eq!(actual, [[1, 2, 3], [4, 5, 6], [7, 8, 9]]);

    let actual: Vec<_> = strided.strided_iter_rev().copied().collect();
    assert_eq!(actual, [[7, 8, 9], [4, 5, 6], [1, 2, 3]]);

    strided
        .strided_iter_mut()
        .enumerate()
        .for_each(|(idx, block)| {
            for val in block {
                *val += idx as isize + 1;
            }
        });

    drop(strided);

    assert_eq!(original, [0, 2, 3, 4, 6, 7, 8, 10, 11, 12, 10]);

    let mut strided = BlockedStridedSliceMut::<3, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 1,
            stride: 3.try_into().unwrap(),
            count: 3.try_into().unwrap(),
        },
    )
    .unwrap();

    strided
        .strided_iter_mut_rev()
        .enumerate()
        .for_each(|(idx, block)| {
            for val in block {
                *val -= idx as isize + 1;
            }
        });

    drop(strided);

    assert_eq!(original, [0, -1, 0, 1, 4, 5, 6, 9, 10, 11, 10]);
}

#[test]
fn strided_iter_with_blocksize_not_dividing_stride() {
    let mut original = [
        0isize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    ];

    let mut strided = BlockedStridedSliceMut::<3, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 2,
            stride: 5.try_into().unwrap(),
            count: 3.try_into().unwrap(),
        },
    )
    .unwrap();

    assert_eq!(strided.strided_len(), 3);
    assert!(!strided.strided_is_empty());

    let actual: Vec<_> = strided.strided_iter().copied().collect();
    assert_eq!(actual, [[2, 3, 4], [7, 8, 9], [12, 13, 14]]);

    let actual: Vec<_> = strided.strided_iter_rev().copied().collect();
    assert_eq!(actual, [[12, 13, 14], [7, 8, 9], [2, 3, 4]]);

    strided
        .strided_iter_mut()
        .enumerate()
        .for_each(|(idx, block)| {
            for val in block {
                *val += idx as isize + 1;
            }
        });

    drop(strided);

    assert_eq!(
        original,
        [
            0, 1, 3, 4, 5, 5, 6, 9, 10, 11, 10, 11, 15, 16, 17, 15, 16, 17
        ]
    );

    let mut strided = BlockedStridedSliceMut::<3, _>::new(
        &mut original,
        StridedAccessArgs {
            offset: 2,
            stride: 5.try_into().unwrap(),
            count: 3.try_into().unwrap(),
        },
    )
    .unwrap();

    strided
        .strided_iter_mut_rev()
        .enumerate()
        .for_each(|(idx, block)| {
            for val in block {
                *val -= idx as isize + 1;
            }
        });

    drop(strided);

    assert_eq!(
        original,
        [0, 1, 0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 15, 16, 17]
    );
}

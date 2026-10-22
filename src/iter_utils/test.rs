use super::StridedIter;
use super::StridedRandomAccess;
#[test]
fn strided_iter_for_empty_slice() {
    let mut data: &mut [f32] = &mut [];

    assert!(data.strided_iter().next().is_none());
    assert!(data.strided_iter_rev().next().is_none());
    assert!(data.strided_iter_mut().next().is_none());
    assert!(data.strided_iter_mut_rev().next().is_none());
    assert!(data.strided_is_empty());
    assert_eq!(data.strided_len(), 0);
}

#[test]
fn strided_random_access_for_empty_slice() {
    let mut data: &mut [f32] = &mut [];

    assert!(data.strided_first().is_none());
    assert!(data.strided_last().is_none());
    assert!(data.strided_2nd_last().is_none());
    // no-ops
    data.strided_write_first(1.);
    data.strided_write_last(1.);
}

#[test]
fn strided_iter_for_one_element_slice() {
    let mut array = [1f64];
    let mut data = array.as_mut_slice();

    let expected = data.to_vec();

    assert_eq!(data.strided_iter().copied().collect::<Vec<_>>(), expected);
    assert_eq!(
        data.strided_iter_rev().copied().collect::<Vec<_>>(),
        expected
    );

    data.strided_iter_mut().for_each(|v| *v += 1.);
    assert_eq!(data.strided_iter().copied().collect::<Vec<_>>(), vec![2.]);
    data.strided_iter_mut_rev().for_each(|v| *v *= 2.);
    assert_eq!(data.strided_iter().copied().collect::<Vec<_>>(), vec![4.]);

    assert_eq!(data.strided_len(), 1);
    assert!(!data.strided_is_empty());
}

#[test]
fn strided_random_access_for_one_element_slice() {
    // NOTE: scoped to pass borrow checking
    let mut array = [1f64];
    {
        let mut data = array.as_mut_slice();
        assert_eq!(data.strided_first(), Some(1.));
        assert_eq!(data.strided_last(), Some(1.));
        assert_eq!(data.strided_2nd_last(), None);

        data.strided_write_first(123.);
    }
    assert_eq!(array, [123.]);
    let mut data = array.as_mut_slice();
    data.strided_write_last(-123.);
    assert_eq!(array, [-123.]);
}

#[test]
fn strided_iter_for_multi_element_slice() {
    // NOTE: scoped to pass borrow checking
    let mut array = [1f64, 3., 5.];
    {
        let mut data = array.as_mut_slice();
        assert_eq!(data.strided_len(), 3);
        assert!(!data.strided_is_empty());
        assert_eq!(
            data.strided_iter().copied().collect::<Vec<_>>(),
            vec![1., 3., 5.]
        );
        assert_eq!(
            data.strided_iter_rev().copied().collect::<Vec<_>>(),
            vec![5., 3., 1.]
        );
        data.strided_iter_mut().enumerate().for_each(|(idx, val)| {
            *val += 1. + idx as f64;
        });
    }

    assert_eq!(array, [2., 5., 8.]);
    {
        let mut data = array.as_mut_slice();
        data.strided_iter_mut_rev()
            .enumerate()
            .for_each(|(idx, v)| *v *= idx as f64);
    }
    assert_eq!(array, [4., 5., 0.]);
}

#[test]
fn strided_random_access_for_multi_element_slice() {
    let mut array = [1f64, 3., 5.];
    {
        let mut data = array.as_mut_slice();
        assert_eq!(data.strided_first(), Some(1.));
        assert_eq!(data.strided_last(), Some(5.));
        assert_eq!(data.strided_2nd_last(), Some(3.));

        data.strided_write_first(123.);
    }
    assert_eq!(array, [123., 3., 5.]);
    let mut data = array.as_mut_slice();
    data.strided_write_last(-123.);
    assert_eq!(array, [123., 3., -123.]);
}

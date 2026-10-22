#![allow(unsafe_code)]
use core::num::NonZero;

use crate::iter_utils::{StridedAccessArgs, StridedSliceMut};

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

    unsafe {
        assert!(
            StridedSliceMut::<f32>::new(
                &mut [],
                StridedAccessArgs {
                    offset: 1,
                    stride: 1.try_into().unwrap(),
                    count: NonZero::new_unchecked(0)
                }
            )
            .is_none()
        )
    }

    unsafe {
        assert!(
            StridedSliceMut::<f32>::new(
                &mut [],
                StridedAccessArgs {
                    offset: 1,
                    stride: NonZero::new_unchecked(0),
                    count: NonZero::new_unchecked(0)
                }
            )
            .is_none()
        )
    }
    todo!("waah")
}

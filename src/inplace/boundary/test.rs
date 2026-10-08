use crate::inplace::TransmittableBoundaryExtension;

#[test]
fn extensions_produce_reference_results() {
    // we have to use n=2 as the start because the reference impl has a bug for
    // n = 1 in the whole sym extension
    for n in [2, 1000, 100000] {
        for idx in -(4 * n as isize)..(4 * n as isize) {
            let actual = TransmittableBoundaryExtension::Mirrored.extended_index(idx, n);
            let expected = reference::half_sym_ext(n as _, idx as _);
            assert_eq!(
                expected,
                actual.try_into().unwrap(),
                "half sym extension does not match"
            );

            let actual = TransmittableBoundaryExtension::WholeSymmetric.extended_index(idx, n);
            let expected = reference::whole_sym_ext(n as _, idx as _);
            assert_eq!(
                expected,
                actual.try_into().unwrap(),
                "half sym extension does not match"
            );

            let actual = TransmittableBoundaryExtension::Periodic.extended_index(idx, n);
            let expected = reference::periodic_ext(n as _, idx as _);
            assert_eq!(
                expected,
                actual.try_into().unwrap(),
                "half sym extension does not match"
            );
        }
    }
}

#[test]
// just to make sure I don't have the same bug as the reference code
fn extensions_dont_stall_for_n1() {}

// semantic 1:1 translations from the reference implementation
mod reference {
    // WARN: this implementation has a bug for n==1 where it loops indefinitely
    pub fn half_sym_ext(n: i64, mut i: i64) -> i64 {
        loop {
            if i < 0 {
                i = -1 - i;
            } else if i >= n {
                i = (2 * n - 1) - i;
            } else {
                return i;
            }
        }
    }

    // WARN: this implementation has a bug for n==1 where it loops indefinitely
    pub fn whole_sym_ext(n: i64, mut i: i64) -> i64 {
        loop {
            if i < 0 {
                i = -i;
            } else if i >= n {
                i = (2 * n - 2) - i;
            } else {
                return i;
            }
        }
    }

    // WARN: this implementation has a bug for n==1 where it loops indefinitely
    pub fn periodic_ext(n: i64, mut i: i64) -> i64 {
        loop {
            if i < 0 {
                i = n + i;
            } else if i >= n {
                i = i - n;
            } else {
                return i;
            }
        }
    }
}

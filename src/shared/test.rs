use crate::shared::Coord;

#[test]
fn coordinates_are_equivalent_in_xy_or_row_col_notation() {
    assert_eq!(Coord::row(-19.4).col(22.5), Coord::x(22.5).y(-19.4));
}

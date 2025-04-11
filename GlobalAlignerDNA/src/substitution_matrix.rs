use lazy_static::lazy_static;
use ndarray;

lazy_static! {
    static ref SUBMAT: ndarray::Array2<i32> = ndarray::Array::from_shape_vec(
        (15, 15),
        vec![
            5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, -4, 5, -4, -4, -4, 1, -4, 1, 1,
            -4, -1, -4, -1, -1, -2, -4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2, -4, -4,
            -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, -4, -4, 1, 1, -1, -4, -2, -2, -2, -2,
            -1, -1, -3, -3, -1, 1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, 1, -4, 1,
            -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -4, 1, -4, 1, -2, -2, -4, -1, -2, -2,
            -1, -3, -1, -3, -1, -4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, 1, -4,
            -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -4, -1, -1, -1, -1, -3, -3, -1, -1,
            -3, -1, -2, -2, -2, -1, -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, -1,
            -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, -1, -1, -1, -4, -3, -1, -1, -3,
            -1, -3, -2, -2, -2, -1, -1, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        ]
    )
    .expect("Unable to load EDNAFULL");
}

fn query(key: u8) -> usize {
    match key {
        b'A' => 0,
        b'T' => 1,
        b'G' => 2,
        b'C' => 3,
        b'S' => 4,
        b'W' => 5,
        b'R' => 6,
        b'Y' => 7,
        b'K' => 8,
        b'M' => 9,
        b'B' => 10,
        b'V' => 11,
        b'H' => 12,
        b'D' => 13,
        b'N' => 14,
        _ => unreachable!(),
    }
}

pub fn edna_scoring(symbol_1: u8, symbol_2: u8) -> i32 {
    let index1 = query(symbol_1);
    let index2 = query(symbol_2);

    SUBMAT[(index1, index2)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edna() {
        let score_1 = edna_scoring(b'A', b'A');
        assert_eq!(score_1, 5);
        let score_1 = edna_scoring(b'T', b'T');
        assert_eq!(score_1, 5);
        let score_1 = edna_scoring(b'G', b'G');
        assert_eq!(score_1, 5);
        let score_1 = edna_scoring(b'C', b'C');
        assert_eq!(score_1, 5);
        let score_1 = edna_scoring(b'A', b'T');
        assert_eq!(score_1, -4);
        let score_1 = edna_scoring(b'G', b'T');
        assert_eq!(score_1, -4);
        let score_1 = edna_scoring(b'T', b'G');
        assert_eq!(score_1, -4);
    }
}

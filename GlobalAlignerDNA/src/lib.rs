#[allow(dead_code)]
use std::collections::HashSet;
pub mod substitution_matrix;
use substitution_matrix::edna_scoring;
use std::{fmt, ops::Index, ops::IndexMut};




/// needleman wunsch algotithm with affine gap
/// here we use the EDNAFULL Matrix: for the match/mismatch scoring.
/// At the moment there is no option to use others matrix.
/// how to use:
///  ```rust
/// use GlobalAlignerDNA::{AlignmentOperation, needleman_wunsch};
/// let a = "AATTTTTGGG";
/// let b = "AATGGGG";
/// let p = needleman_wunsch(a, b, -10., -0.5);
/// assert_eq!(p.score, 8.);
/// assert_eq!(p.operation, vec![AlignmentOperation::Match(3), AlignmentOperation::Insertion(4), AlignmentOperation::Match(3), AlignmentOperation::Deletion(1)]);
/// ```
///
pub fn needleman_wunsch(seq1: &str, seq2: &str, gap_open: f32, gap_extend: f32) -> AlignerResult {
    let query: &[u8] = seq1.as_bytes();
    let target: &[u8] = seq2.as_bytes();

    let size_query = seq1.len();
    let size_target = seq2.len();

    let aln_mat = compute_matrix(query, target, gap_open, gap_extend, edna_scoring);
    AlignerResult {
        score: aln_mat[(size_query, size_target)],
        operation: backtrack(&aln_mat, query, target),
    }
}

/// Store the alignment operation.
/// Can be used to reconstruct the alignment String
/// Or count the % Match.
/// PartialEq only for test
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignmentOperation {
    Deletion(usize),
    Insertion(usize),
    Match(usize),
    Mismatch(usize),
}




/// Hold the result of the alignment.
/// The score is defined as the global alignment score.
/// the operation indicate how to reconstruct the alignment, it is similar to CIGAR string.
pub struct AlignerResult {
    pub score: f32,
    pub operation: Vec<AlignmentOperation>,
}

/// To be able to pretty print the alignment results.
impl fmt::Display for AlignerResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "alignment score: {}, alignment CIGAR {}",
            self.score,
            self.get_cigar()
        )
    }
}

impl AlignerResult {
    /// Return a compact Cigar String of the Alignment.
    pub fn get_cigar(&self) -> String {
        self.operation
            .iter()
            .map(|operation| match operation {
                AlignmentOperation::Deletion(n) => format!("D{}", n),
                AlignmentOperation::Insertion(n) => format!("I{}", n),
                AlignmentOperation::Match(n) => format!("M{}", n),
                AlignmentOperation::Mismatch(n) => format!("X{}", n),
            })
            .collect::<Vec<String>>()
            .join("")
    }

    /// Return the percentage of identity of the alignment.
    /// It is defined only over the aligned part of the alignment:
    /// number of match / (number of match + number of mismatch)
    pub fn get_percent_identity(&self) -> f32 {
        (self.get_match_n() as f32) / (self.aln_len() as f32)
    }

    /// Return the alignment length, because this is a global alignment.
    /// It is defined as number of match + number of mismatch
    fn aln_len(&self) -> u32 {
        let mut result = 0;
        for op in &self.operation {
            match op {
                AlignmentOperation::Match(n) => {
                    result += n;
                }
                AlignmentOperation::Mismatch(n) => {
                    result += n;
                }
                _ => (),
            }
        }
        result.try_into().unwrap()
    }
    /// Return the number of match.
    fn get_match_n(&self) -> u32 {
        let mut result = 0;
        for op in &self.operation {
            if let AlignmentOperation::Match(n) = op {
                result += n;
            };
        }
        result.try_into().unwrap()
    }
}

/// The algorithm first step is to fill as matrix.
/// based on two sequences and on alignment parameters
/// here we use the EDNAFULL Matrix: for the match/mismatch scoring.
/// At the moment there is no option to use others matrix.
/// TODO implement a Param struct/trait that hold the parameter information.
/// gap open and gap extend expect negative value null or positive value will lead to anything but an alignment.
///
fn compute_matrix(
    query: &[u8],
    target: &[u8],
    gap_open: f32,
    gap_extend: f32,
    scoring_function: fn(u8, u8) -> i32,
) -> Matrix {
    let size_query = query.len();
    let size_target = target.len();

    let mut aln_mat = Matrix::new(size_query + 1, size_target + 1, 0.);
    // We use an hashset to keep track of the position where a gap was open.
    // For now that the easiest way I found, but it may be far from optimal.
    let mut gap_set = HashSet::new();

    // Core of the algorithm see for details:
    // https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    for row_i in 1..size_query + 1 {
        for col_j in 1..size_target + 1 {
            let diag_score = scoring_function(query[row_i - 1], target[col_j - 1]) as f32
                + aln_mat[(row_i - 1, col_j - 1)];

            let top_score = {
                let top_is_gap = if gap_set.contains(&(row_i - 1, col_j)) {
                    gap_extend
                } else {
                    gap_open
                };

                aln_mat[(row_i - 1, col_j)] + top_is_gap
            };

            let left_score = {
                let left_is_gap = if gap_set.contains(&(row_i, col_j - 1)) {
                    gap_extend
                } else {
                    gap_open
                };

                aln_mat[(row_i, col_j - 1)] + left_is_gap
            };

            let best_top_left_score = top_score.max(left_score);

            aln_mat[(row_i, col_j)] = {
                if diag_score >= best_top_left_score {
                    diag_score
                } else {
                    gap_set.insert((row_i, col_j));
                    best_top_left_score
                }
            };
        }
    }
    aln_mat
}

/// Perform the backtracking to recover the alignment
/// return a vec<AlignmentOperation>
///
fn backtrack(aln_mat: &Matrix, query: &[u8], target: &[u8]) -> Vec<AlignmentOperation> {
    let mut i = query.len();
    let mut j = target.len();

    let mut vec_operation: Vec<AlignmentOperation> = Vec::new();
    let mut previous_op: AlignmentOperation = AlignmentOperation::Match(0);
    let mut current_op: AlignmentOperation = AlignmentOperation::Match(0);
    while (i > 0) || (j > 0) {
        match (i > 0, j > 0) {
            (true, true) => {
                if (aln_mat[(i - 1, j - 1)] >= aln_mat[(i, j - 1)])
                    & (aln_mat[(i - 1, j - 1)] >= aln_mat[(i - 1, j)])
                {
                    if query[i - 1] == target[j - 1] {
                        current_op = AlignmentOperation::Match(1);
                    } else {
                        current_op = AlignmentOperation::Mismatch(1);
                    }
                    i -= 1;
                    j -= 1;
                } else if aln_mat[(i - 1, j)] > aln_mat[(i, j - 1)] {
                    current_op = AlignmentOperation::Insertion(1);
                    i -= 1;
                } else {
                    current_op = AlignmentOperation::Deletion(1);
                    j -= 1;
                }
            }
            (false, true) => {
                current_op = AlignmentOperation::Deletion(1);
                j -= 1;
            }
            (true, false) => {
                current_op = AlignmentOperation::Insertion(1);
                i -= 1;
            }
            (false, false) => (),
        }
        previous_op = match (current_op, previous_op) {
            (current, AlignmentOperation::Match(0)) => current,
            (AlignmentOperation::Match(n), AlignmentOperation::Match(m)) => {
                AlignmentOperation::Match(n + m)
            }
            (AlignmentOperation::Mismatch(n), AlignmentOperation::Mismatch(m)) => {
                AlignmentOperation::Mismatch(n + m)
            }
            (AlignmentOperation::Insertion(n), AlignmentOperation::Insertion(m)) => {
                AlignmentOperation::Insertion(n + m)
            }
            (AlignmentOperation::Deletion(n), AlignmentOperation::Deletion(m)) => {
                AlignmentOperation::Deletion(n + m)
            }
            (current, prev) => {
                vec_operation.push(prev);
                current
            }
        };
    }
    vec_operation.push(previous_op);
    vec_operation.iter().copied().rev().collect()
}


/// in my code I use the full ednafull with all IUPAC support,
/// but this is enough to test the code!
/// symbol_query and symbol_target are symbols
/// i32 for compatibility with rust rust HTS-lib. 
/*fn edna_scoring(symbol_query: u8, symbol_target: u8) -> i32 {
    if symbol_query == symbol_target {
        5
    } else {
        -4
    }
}*/

/// Input validation
/// I really appreciated the comment Matthieu M.
/// And understand that I need input validation
/// But I think this is a bit complex, First DNA can and is represented by more than 4 letter this is defined by the IUPAC
/// But shortly when you know: You have either a A or C you will represent it with a M.
/// full list here: https://www.bioinformatics.org/sms/iupac.html
/// Lastly I want this code to be compatible with rust HTS-lib which uses u8...
///
/// Lastly nothing should prevent me to uses this code from other sequences type (like Protein)
/// or event custom set of character in the future
/// So Keeping your word in mind I propose this approach:
pub trait AlphabetTrait {
    fn is_pure(&self, sequence: &str) -> bool;
}

pub struct Alphabet{
    sets: HashSet<char>,
}

impl Alphabet{
    pub fn new(sets: HashSet<char>) -> Self{
        Alphabet{sets}
    }

    pub fn new_dna_alphabet() -> Self{
        let mut sets: HashSet<char> = HashSet::new();
        for e in  ['A',  'C', 'G', 'T', 'U', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', 
        'a',  'c', 'g', 't', 'u', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']{
            sets.insert(e);
        }
        Alphabet::new(sets)
    }
}

impl AlphabetTrait for Alphabet{
    fn is_pure(&self, sequence: &str) -> bool {
        sequence.chars().all(|x| self.sets.contains(&x))
    }
}


/// A row first minimal matrix implementations
struct Matrix {
    data: Vec<f32>,
    rows: usize,
}

impl Matrix {
    /// Create a matrix of size rows * cols filled with default value
    fn new(rows: usize, cols: usize, default_value: f32) -> Self {
        Matrix {
            data: vec![default_value; rows * cols],
            rows,
        }
    }
}

/// Uses a (row, col) for indexing the Matrix, I am unsure if this is a good approach
/// With it we can do Matrix[(row_index, col_index)]
impl Index<(usize, usize)> for Matrix {
    type Output = f32;

    fn index(&self, position: (usize, usize)) -> &Self::Output {
        &self.data[position.0 * self.rows + position.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, position: (usize, usize)) -> &mut Self::Output {
        &mut self.data[position.0 * self.rows + position.1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn needleman_wunsch_() {
        let a = "AATTTTTGGG";
        let b = "AATGGGG";
        let p = needleman_wunsch(a, b, -10., -0.5);
        println!("{}", p);
        println!("{:?}", p.operation);
        assert_eq!(p.score, 8.);
        assert_eq!(
            p.operation,
            vec![
                AlignmentOperation::Match(3),
                AlignmentOperation::Insertion(4),
                AlignmentOperation::Match(3),
                AlignmentOperation::Deletion(1)
            ]
        );
        assert_eq!(true, true);
    }
}

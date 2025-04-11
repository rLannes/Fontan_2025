#![allow(deprecated)]
#[allow(dead_code)]
use pyo3::prelude::*;
use rayon::prelude::*;
use std::sync::Mutex;
use GlobalAlignerDNA::needleman_wunsch;

#[pyfunction]
/// chunck seq1 and seq2 by the windows, the last windows is the remainder
/// optionally you can change alignment parameter,
/// This uses global alignment.
/// be carefull picking the windows size, as this function will return all value.
/// and a too small windows will end up creating a very large matrix.  
/// by default use match and mismatch are derived from EDNAFULL.
///  gap open: -10 gap extend: -0.5
/// number of thread to use
/// match_percent if tru use percent of match instead of score aln.
fn get_matrix(
    seq1: String,
    seq2: String,
    windows_size: usize,
    gap_open: f32,
    gap_extend: f32,
    thread: usize,
    match_percent: bool,
) -> PyResult<Vec<Vec<f32>>> {

    let seq1 = seq1.to_uppercase();
    let seq1_chuncks = &seq1.as_bytes().chunks(windows_size).map(std::str::from_utf8)
    .collect::<Result<Vec<&str>, _>>()
    .expect("unable to chunck sequence 1");

    let seq2 = seq2.to_uppercase();
    let seq2_chuncks = &seq2.as_bytes().chunks(windows_size).map(std::str::from_utf8)
    .collect::<Result<Vec<&str>, _>>()
    .expect("unable to chunck sequence 2");


    let size_s1 = seq1_chuncks.len();
    let size_s2 = seq2_chuncks.len();

    let mut result: Vec<Vec<f32>> =  Vec::with_capacity(size_s1);
    for _i in 0..size_s1{
        result.push(vec![0.; size_s2]);
    }

    //rayon::ThreadPoolBuilder::new().num_threads(thread).build_global().or_else();
    let pool = rayon::ThreadPoolBuilder::new().num_threads(thread).build().unwrap();
    // 
    let mut result = Mutex::new(result);
    pool.install( || {(0..size_s1).into_par_iter().for_each(|row_i| {
    //for row_i in 0..size_s1{
        for col_j in 0..size_s2{
            let value = needleman_wunsch(&seq1_chuncks[row_i], &seq2_chuncks[col_j], gap_open, gap_extend);
            let mut f = value.score;
            if match_percent{
                f = value.get_percent_identity()
            }
            
            let mut result = result.lock().unwrap(); 
            result[row_i][col_j] = f;
        }
    }) });
    /*     
    for row_i in 0..size_s1{
        for col_j in 0..size_s2{
            result[row_i][col_j] = needleman_wunsch(&seq1_chuncks[row_i], &seq2_chuncks[col_j], gap_open, gap_extend).score;
        }
    }; */
    Ok(result.into_inner().expect("unabe to free mutex"))
}


/// A Python module implemented in Rust.
#[pymodule]
fn Rust_alnPairMat(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_matrix, m)?)?;
    Ok(())
}
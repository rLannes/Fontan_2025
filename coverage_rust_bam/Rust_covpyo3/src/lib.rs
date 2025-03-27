use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use rust_htslib::bam::{IndexedReader, Read, Header};
use rust_htslib::bam;
use std::collections::HashMap;
use strand_specifier_lib::{Strand, LibType, check_flag};
use std::str::FromStr;
use std::cmp;
use rust_htslib::bam::record::Record;
use CigarParser::cigar::Cigar;

#[pyfunction]
fn get_header(bam_path: String) -> PyResult<Vec<String>>{

    let bam = bam::Reader::from_path(&bam_path).unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut seq: Vec<String> = Vec::new(); 

    for (key, records) in header.to_hashmap() {
            for record in records {
                if record.contains_key("SN"){
                    seq.push(record["SN"].to_string());
                }
            }
    }
    Ok(seq)
}

#[pyfunction]
fn get_coverage(start:i64, end:i64, chrom: String, strand: String,
     bam_path: String, lib: String, mapq_thr: u8) -> PyResult<Vec<u32>>{
    ///
    /// 
    /// retrieve the Primary aligned reads at between start and end in a given bam
    /// at the moment only rf stranded librairy are supported.
    /// I will update fr ans non starnde dif I need it or people ask for it.
    /// 
    //let mut map: HashMap<i64, u32> = HashMap::with_capacity((end - start) as usize);
    let mut container = vec![0; (end - start) as usize];
    let mut bam = IndexedReader::from_path(&bam_path).unwrap();

    let lib_type = LibType::from(lib.as_str());
    let strand_feature = Strand::from(strand.as_str());

    bam.fetch((&chrom, start, end)).unwrap();
    let mut read_strand: Strand = Strand::Plus; 
    let mut cpt = 0;
    for p in bam.pileup() {
        let pileup = p.unwrap();
        cpt = 0;
        if start <= i64::from(pileup.pos()) && i64::from(pileup.pos()) < end {
            for alignment in pileup.alignments() {
                if !alignment.is_del() && !alignment.is_refskip() && !check_flag(alignment.record().flags(), 256, 0) && !(alignment.record().mapq() < mapq_thr){
                    
                    if strand_feature == Strand::NA{
                        cpt += 1;
                        continue
                    }

                    if let Some(read_strand) = lib_type.get_strand(alignment.record().flags()){
                        if read_strand == Strand::NA {
                            cpt += 1;
                            continue
                        }
                        else if strand_feature == read_strand{
                            cpt += 1;
                            continue
                        }
                    }
                }
            }

            container[(pileup.pos() as i64 - start as i64) as usize ] = cpt;
            
        }

    } 

    return Ok(container)

}

#[pyfunction]
fn get_coverage_algo2(start:i64, end:i64, chrom: String, strand: String,
     bam_path: String, lib: String, mapq_thr: u8, flag_in: u16, flag_exclude: u16) -> PyResult<Vec<u32>>{
    ///
    let mut container = vec![0; (end - start) as usize];
    let mut bam = IndexedReader::from_path(&bam_path).unwrap();

    let lib_type = LibType::from(lib.as_str());
    let strand_feature = Strand::from(strand.as_str());

    bam.fetch((&chrom, start, end)).unwrap();
    let mut read_strand: Strand = Strand::Plus; 
    let mut cpt = 0;

    let mut record: Record;
    let mut pos_s: i64;
    let mut pos_e: i64;
    let mut cig: Cigar;
    let mut flag: u16;

    for p in bam.records() {
        record = p.unwrap();

        pos_s = record.pos();
        cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
        pos_e = cig.get_end_of_aln(&pos_s);
        flag = record.flags();
        // flag_in: u16, flag_exclude
        // !check_flag(flag, 256, 0) && !(record.mapq() < mapq_thr){
        if check_flag(flag, flag_in, flag_exclude) && (mapq_thr == 0 || !(record.mapq() < mapq_thr)){
            
            if let Some(read_strand) = lib_type.get_strand(flag){
                if strand_feature == read_strand{
                    cover_from_intervall(&mut container, start, end, cig.get_reference_cover(pos_s));
                }
            }
        }
    } 
    return Ok(container)
}

pub fn cover_from_intervall(feature_cover : &mut Vec<u32>,
    cover_start: i64,
    cover_end: i64, 
    read_cover: Vec<i64>) -> (){
    
    let mut end: i64 = 0;
    
    for (i, start) in read_cover.iter().enumerate().step_by(2){
         end = read_cover[i + 1];
         if !((*start > cover_end) | (end < cover_start)){
            for ii in (*cmp::max(start, &cover_start) as usize)..(*cmp::min(&end, &cover_end) as usize){
                 //println!("{:?} {:?} {:?}", ii, cmp::max(start, &cover_start), cmp::min(&end, &cover_end));
                 feature_cover[ii - cover_start as usize] += 1;
             } 
             
         }
    }
()
}


/// A Python module implemented in Rust.
#[pymodule]
fn Rust_covpyo3(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(get_header, m)?)?;
    m.add_function(wrap_pyfunction!(get_coverage_algo2, m)?)?;
    Ok(())
}

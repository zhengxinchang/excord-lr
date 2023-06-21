use std::collections::HashMap;


#[derive(Debug)]
pub struct AlignmentPos {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub cigar_map: HashMap<char, u32>,
    pub strand: i32, // true = forward
    pub mapq: u8,
    pub raw_cigar: String,
}

impl AlignmentPos {
    pub fn new(
        chrom: &str,
        start: &i64,
        cigar_map: HashMap<char, u32>,
        strand: &i32,
        mapq: &u8,
        cigar_string: &str,
    ) -> AlignmentPos {
        let end =
            start + (*cigar_map.get(&'D').unwrap()) as i64 + (*cigar_map.get(&'M').unwrap()) as i64
                - 1i64;
        let chrom_clean: String;
        if chrom.starts_with(&"chr".to_string()) {
            chrom_clean = chrom.strip_prefix(&"chr".to_string()).unwrap().to_string();
        } else {
            chrom_clean = chrom.to_string();
        }

        AlignmentPos {
            chrom: String::from(chrom_clean),
            start: *start,
            end: end + 1,
            cigar_map: cigar_map,
            strand: *strand,
            mapq: *mapq,
            raw_cigar: cigar_string.to_string(),
        }
    }
}
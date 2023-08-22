
/// # Type of SVs in alignment event
/// 
/// Have not beed used.
/// 
#[derive(Debug,Clone,PartialEq)]
pub enum AlignEventType {
    Ins,
    Del
}

#[derive(Debug,Clone)]
/// struct for alignment event
/// left and right reference consume will sum the length of D M = X for the left/right part of d
pub struct AlignmentEvent {
    pub lchrom: String,
    pub lstart: u32,
    pub lend: u32,
    pub lstrand: i32,
    pub rchrom: String,
    pub rstart: u32,
    pub rend: u32,
    pub rstrand: i32, // true = forward
    pub events_num: i32,
    pub svtype:AlignEventType
}

impl AlignmentEvent {
    pub fn new(
        chrom: &str,
        left_consume: &u32,
        right_consume: &u32,
        event_len: &u32,
        pos: &i64,
        strand: &i32,
        sv_type: Option<AlignEventType>
    ) -> AlignmentEvent {
        let chrom_clean: String;
        if chrom.starts_with(&"chr".to_string()) {
            chrom_clean = chrom.strip_prefix(&"chr".to_string()).unwrap().to_string();
        } else {
            chrom_clean = chrom.to_string();
        }
        let pos2 = *pos as u32;
        AlignmentEvent {
            lchrom: chrom_clean.clone(),
            lstart: pos2,
            lend: pos2 + left_consume,
            lstrand: *strand,
            rchrom: chrom_clean.clone(),
            rstart: pos2 + left_consume + event_len,
            rend: pos2 + left_consume + event_len + right_consume,
            rstrand: *strand,
            events_num: 1i32,
            svtype:sv_type.unwrap()
        }
    }
}

use crate::split_read_event::SplitReadEvent;
use std::{
    cmp::Ordering,
    collections::HashMap,
    env, io,
    path::{Path, PathBuf},
};

/// # get the position of the first match base.
pub fn find_first_match_pos(cigar_str: &str) -> i64 {
    let mut p = 0i64;
    let mut p_string = "".to_string();
    for c in cigar_str.chars() {
        if c != 'D'
            && c != 'M'
            && c != 'I'
            && c != 'H'
            && c != 'S'
            && c != 'P'
            && c != 'X'
            && c != '='
            && c != 'N'
        {
            p_string += &c.to_string();
        } else {
            if c == 'M' {
                break;
            } else {
                p += p_string.parse::<i64>().unwrap();
                p_string = "".to_string();
            }
        }
    }
    // dbg!(cigar_str, p);
    return p;
}

/// # Compare two Alignment
///
/// From brentp:
///
/// ```
/// output splitters. Splitters are ordered by their offset into the read.
/// given, cigars of:
/// A:20S30M100S
/// B:50S30M50S
/// C:90S30M30S
/// we would order them as they are listed. We would output bedpe intervals
/// for A-B, and B-C
/// ```
pub fn splitter_order_cmp(a: &SplitReadEvent, b: &SplitReadEvent) -> Ordering {
    let a_first_match_pos = find_first_match_pos(&a.raw_cigar);
    let b_first_match_pos = find_first_match_pos(&b.raw_cigar);
    a_first_match_pos.cmp(&b_first_match_pos)
}

pub fn alignment_pos_cmp(a: &SplitReadEvent, b: &SplitReadEvent) -> Ordering {
    if a.chrom.cmp(&b.chrom) != Ordering::Equal {
        a.chrom.as_bytes().cmp(&b.chrom.as_bytes())
    } else {
        a.start.cmp(&b.start)
    }
}

pub fn parse_cigar(cigar_str: &str) -> HashMap<char, u32> {
    let _cigar_str = cigar_str.to_owned();

    let mut n_str = String::new();
    let mut cigar_map: HashMap<char, u32> = HashMap::from([
        ('D', 0),
        ('M', 0),
        ('I', 0),
        ('H', 0),
        ('S', 0),
        ('P', 0),
        ('X', 0),
        ('=', 0),
        ('N', 0),
    ]);

    for x in _cigar_str.chars() {
        if x.is_numeric() {
            n_str += &x.to_string();
        } else {
            // println!("{}",n_str);
            let n: u32 = n_str.parse::<u32>().unwrap();
            let previous_n = cigar_map.entry(x).or_insert(0);
            *previous_n += n;
            n_str = "".to_string();
        }
    }

    cigar_map
}

pub fn parse_supplementary_alignment(s: &str) -> SplitReadEvent {
    let sa_vec: Vec<&str> = s.split(',').collect();
    let chrom = sa_vec[0];
    let pos = sa_vec[1].parse::<i64>().unwrap();
    let strand = match sa_vec[2] {
        "+" => Ok(1),
        "-" => Ok(-1),
        _ => Err("No strand info"),
    };
    let cigar_str = parse_cigar(sa_vec[3]);
    let mapq = sa_vec[4].parse::<u8>().unwrap();
    let _nm = sa_vec[5].parse::<i64>().unwrap();
    SplitReadEvent::new(
        chrom,
        &(pos - 1),
        cigar_str,
        &strand.unwrap(),
        &mapq,
        &sa_vec[3],
    )
}

pub fn absolute_path(path: impl AsRef<Path>) -> io::Result<PathBuf> {
    let path = path.as_ref();
    let absolute_path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir()?.join(path)
    };

    Ok(absolute_path)
}

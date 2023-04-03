use bio_types::{genome::AbstractInterval, strand::ReqStrand};
use clap::Parser;
use rust_htslib::{
    bam,
    bam::{
        ext::BamRecordExtensions,
        record::{Aux, Cigar, CigarStringView},
        Read, Record,
    },
};
use std::{
    cmp::Ordering,
    collections::HashMap,
    env,
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    process::exit,
};
#[derive(Parser, Debug)]
#[command(name = "excord-LR")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.0.1")]
#[command(about = "
Extract Structural Variation Signals from Long-Read BAMs
Contact: Xinchang Zheng <zhengxc93@gmail.com>
", long_about = None)]
struct Cli {
    /// Path to BAM file
    #[arg(short, long)]
    bam: PathBuf,

    /// Minimal MapQ
    #[arg(short = 'F', long, default_value_t = 1)]
    mapq: u8,

    /// Exclude Flags
    #[arg(short = 'Q', long, default_value_t = 1540)]
    excludeflag: u16,

    /// Threads
    #[arg(short, long, default_value_t = 8)]
    thread: usize,

    /// Output file name
    #[arg(short, long)]
    out: PathBuf,
}

#[derive(Debug)]
struct AlignmentPos {
    chrom: String,
    start: i64,
    end: i64,
    cigar_map: HashMap<char, i32>,
    strand: i32, // true = forward
    mapq: u8,
    raw_cigar: String,
}

impl AlignmentPos {
    fn new(
        chrom: &str,
        start: &i64,
        cigar_map: HashMap<char, i32>,
        strand: &i32,
        mapq: &u8,
        cigar_string: &str,
    ) -> AlignmentPos {
        let end =
            start + (*cigar_map.get(&'D').unwrap()) as i64 + (*cigar_map.get(&'M').unwrap()) as i64
                - 1i64;
        let mut chrom_clean: String;
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

fn find_first_match_pos(cigar_str: &str) -> i64 {
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
    dbg!(cigar_str, p);
    return p;
}

// From brentp:
// output splitters. Splitters are ordered by their offset into the read.
// given, cigars of:
// A:20S30M100S
// B:50S30M50S
// C:90S30M30S
// we would order them as they are listed. We would output bedpe intervals
// for A-B, and B-C
fn splitter_order_cmp(a: &AlignmentPos, b: &AlignmentPos) -> Ordering {
    // dbg!(a, b, a.chrom.cmp(&b.chrom));

    // let a_raw_cigar_string = a.raw_cigar.clone();
    // let b_raw_cigar_string = b.raw_cigar.clone();
    let a_first_match_pos = find_first_match_pos(&a.raw_cigar);
    let b_first_match_pos = find_first_match_pos(&b.raw_cigar);
    a_first_match_pos.cmp(&b_first_match_pos)
}

fn alignment_pos_cmp(a: &AlignmentPos, b: &AlignmentPos) -> Ordering {
    // dbg!(a,b,a.chrom.cmp(&b.chrom));

    if a.chrom.cmp(&b.chrom) != Ordering::Equal {
        a.chrom.as_bytes().cmp(&b.chrom.as_bytes())
    } else {
        a.start.cmp(&b.start)
    }
}

fn parse_cigar(cigar_str: &str) -> HashMap<char, i32> {
    let _cigar_str = cigar_str.to_owned();

    let mut n_str = String::new();
    let mut cigar_map: HashMap<char, i32> = HashMap::from([
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
            let n: i32 = n_str.parse::<i32>().unwrap();
            let previous_n = cigar_map.entry(x).or_insert(0);
            *previous_n += n;
            n_str = "".to_string();
        }
    }

    cigar_map
}

fn parse_supplementary_alignment(s: &str) -> AlignmentPos {
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
    AlignmentPos::new(
        chrom,
        &(pos - 1),
        cigar_str,
        &strand.unwrap(),
        &mapq,
        &sa_vec[3],
    )
}

fn absolute_path(path: impl AsRef<Path>) -> io::Result<PathBuf> {
    let path = path.as_ref();
    let absolute_path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir()?.join(path)
    };

    Ok(absolute_path)
}

fn main() {
    let cli = Cli::parse();
    // println!("{:?}", &cli);
    let _bam = cli.bam;
    if !_bam.is_file() {
        println!("Ivalid BAM file path: {} ", _bam.to_str().unwrap());
        exit(1)
    }
    // handel -o option
    let t = cli.out;
    let mut _outprefix_ancestors = t.ancestors().clone();
    _outprefix_ancestors.next();
    let mut _outprefix_parent = _outprefix_ancestors.next().unwrap().to_path_buf();
    let _outprefix_parent_abs = absolute_path(_outprefix_parent).unwrap();
    if !_outprefix_parent_abs.is_dir() {
        println!(
            "Output directory does not exists: {} ",
            _outprefix_parent_abs.to_str().unwrap()
        );
        exit(1);
    }
    let mut f = BufWriter::new(File::create(t).unwrap());
    let mut bam = bam::Reader::from_path(_bam).unwrap();
    bam.set_threads(cli.thread).unwrap();
    let mut record = Record::new();
    // for r in bam.records() {
    while let Some(result) = bam.read(&mut record) {
        match result {
            Err(_) => break,
            Ok(()) => {}
        }
        if record.is_secondary()|| record.is_unmapped() || record.mapq() < cli.mapq || (record.flags() & cli.excludeflag) == 0 {
            continue;
        }
        let st = record.strand().to_owned(); // MUST move out of match, Because of mutable borrow by strand().
        match record.aux("SA".as_bytes()) {
            Ok(_sa) => {
                let mut alignment_vec: Vec<AlignmentPos> = Vec::new();

                /* put the preliminary alignment to  */
                let contig_name = record.contig();
                let pos = record.pos();
                let strand = match st {
                    ReqStrand::Forward => 1,
                    ReqStrand::Reverse => -1,
                };
                let mapq = record.mapq();
                let cigar_stats_nuc = record.cigar_stats_nucleotides();

                let mut cigar_map = HashMap::new();
                cigar_stats_nuc.iter().for_each(|x| {
                    /*
                        Match(u32),    // M
                        Ins(u32),      // I
                        Del(u32),      // D
                        RefSkip(u32),  // N
                        SoftClip(u32), // S
                        HardClip(u32), // H
                        Pad(u32),      // P
                        Equal(u32),    // =
                        Diff(u32),     // X
                        ps:
                        Match(u32)
                                |--------Type of this Enum
                    */

                    match x {
                        (&Cigar::Del(_), &n) => {
                            // println!("Del:is: {}", n);
                            cigar_map.insert('D', n.to_owned());
                        }
                        (&Cigar::Match(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('M', n.to_owned());
                        }
                        (&Cigar::Ins(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('I', n.to_owned());
                        }
                        (&Cigar::RefSkip(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('N', n.to_owned());
                        }
                        (&Cigar::SoftClip(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('S', n.to_owned());
                        }
                        (&Cigar::HardClip(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('H', n.to_owned());
                        }
                        (&Cigar::Pad(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('P', n.to_owned());
                        }
                        (&Cigar::Equal(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('=', n.to_owned());
                        }
                        (&Cigar::Diff(_), &n) => {
                            // println!("Match:is: {}", n);
                            cigar_map.insert('X', n.to_owned());
                        }
                    }
                });

                let mut first_cigar_str = "".to_string();
                // dbg!(&record.cigar());
                record.cigar().iter().for_each(|cigar| {
                    match cigar {
                        &Cigar::Del(n) => {
                            // println!("Del:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"D".to_string();
                        }
                        &Cigar::Match(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"M".to_string();
                        }
                        &Cigar::Ins(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"I".to_string();
                        }
                        &Cigar::RefSkip(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"N".to_string();
                        }
                        &Cigar::SoftClip(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"S".to_string();
                        }
                        &Cigar::HardClip(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"H".to_string();
                        }
                        &Cigar::Pad(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"P".to_string();
                        }
                        &Cigar::Equal(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"=".to_string();
                        }
                        &Cigar::Diff(n) => {
                            // println!("Match:is: {}", n);
                            first_cigar_str += &n.to_owned().to_string();
                            first_cigar_str += &"X".to_string();
                        }
                    }
                });

                alignment_vec.push(AlignmentPos::new(
                    &contig_name,
                    &pos,
                    cigar_map,
                    &strand,
                    &mapq,
                    &first_cigar_str,
                ));

                /* process the supplementary alignments */
                if let Aux::String(sa) = _sa {
                    let sa_list = sa.split(";").collect::<Vec<&str>>();
                    let sa_list_clean = sa_list.iter().filter(|x| x.len() > 0);
                    for single_sa in sa_list_clean {
                        // println!("{}", single_sa);
                        alignment_vec.push(parse_supplementary_alignment(single_sa));
                    }
                }
                dbg!(record.flags());
                alignment_vec.sort_by(|a, b| splitter_order_cmp(a, b));
                // dbg!(alignment_vec);
                // println!("{:?}", alignment_vec);
                // dbg!(&alignment_vec);
                for i in 1..alignment_vec.len() {
                    let j = i - 1;
                    let a: &AlignmentPos = &alignment_vec[j];
                    let b: &AlignmentPos = &alignment_vec[i];
                    let mut bed_line = "".to_string();
                    if alignment_pos_cmp(a, b) == Ordering::Less {
                        bed_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            a.chrom,
                            a.start,
                            a.end,
                            a.strand,
                            b.chrom,
                            b.start,
                            b.end,
                            b.strand,
                            alignment_vec.len() - 1,
                            "excord-lr"
                        );
                    } else {
                        bed_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            b.chrom,
                            b.start,
                            b.end,
                            b.strand,
                            a.chrom,
                            a.start,
                            a.end,
                            a.strand,
                            alignment_vec.len() - 1,
                            "excord-lr"
                        );
                    }

                    f.write(bed_line.as_bytes()).unwrap();

                    // println!(
                    //     "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    //     a.chrom,
                    //     a.start,
                    //     a.end,
                    //     a.strand,
                    //     b.chrom,
                    //     b.start,
                    //     b.end,
                    //     b.strand,
                    //     alignment_vec.len() - 1
                    // );
                }
                alignment_vec.clear();
            }
            Err(_) => {}
        };
    }
}

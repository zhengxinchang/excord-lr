use bio_types::{genome::AbstractInterval, strand::ReqStrand};
use clap::Parser;
use rust_htslib::{
    bam,
    bam::{
        ext::BamRecordExtensions,
        record::{Aux, Cigar},
        Read,
    },
};
use std::{cmp::Ordering, collections::HashMap, path::PathBuf, process::exit};
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
}

impl AlignmentPos {
    fn new(
        chrom: &str,
        start: &i64,
        cigar_map: HashMap<char, i32>,
        strand: &i32,
        mapq: &u8,
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
            end: end,
            cigar_map: cigar_map,
            strand: *strand,
            mapq: *mapq,
        }
    }
}

fn alignment_pos_cmp(a: &AlignmentPos, b: &AlignmentPos) -> Ordering {
    if a.chrom.cmp(&b.chrom) != Ordering::Equal {
        a.chrom.cmp(&b.chrom)
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

    AlignmentPos::new(chrom, &(pos - 1), cigar_str, &strand.unwrap(), &mapq)
}

fn main() {
    let cli = Cli::parse();
    // println!("{:?}", &cli);
    let _bam = cli.bam;
    if !_bam.is_file() {
        println!("Ivalid BAM file path: {} ", _bam.to_str().unwrap());
        exit(1)
    }

    let mut bam = bam::Reader::from_path(_bam).unwrap();

    for r in bam.records() {
        let mut record = r.unwrap();
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
                alignment_vec.push(AlignmentPos::new(
                    &contig_name,
                    &pos,
                    cigar_map,
                    &strand,
                    &mapq,
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

                alignment_vec.sort_by(|a, b| alignment_pos_cmp(a, b));

                // println!("{:?}", alignment_vec);

                for i in 1..alignment_vec.len() {
                    let j = i - 1;
                    let a: &AlignmentPos = &alignment_vec[j];
                    let b: &AlignmentPos = &alignment_vec[i];
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        a.chrom,
                        a.start,
                        a.end,
                        a.strand,
                        b.chrom,
                        b.start,
                        b.end,
                        b.strand,
                        alignment_vec.len() - 1
                    );
                }
            }
            Err(_) => {}
        };

        //
    }
}

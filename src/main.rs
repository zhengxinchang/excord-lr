use bio_types::{genome::AbstractInterval, strand::ReqStrand};
use clap::Parser;
use rust_htslib::{
    bam,
    bam::{
        // ext::BamRecordExtensions,
        record::{Aux, Cigar},
        Read,
        Record,
    },
};
use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::{PathBuf},
    process::exit,
};
mod aligments_pos;
use aligments_pos::AlignmentPos;
mod utils;
use utils::*;
#[derive(Parser, Debug)]
#[command(name = "excord-LR")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.1")]
#[command(about = "
Extract Structural Variation Signals from Long-Read BAMs
Contact: Xinchang Zheng <zhengxc93@gmail.com>
", long_about = None)]
struct Cli {
    /// Path to BAM file
    #[arg(short, long)]
    bam: PathBuf,

    /// Minimal MapQ
    #[arg(short = 'Q', long, default_value_t = 1)]
    mapq: u8,

    /// Exclude Flags
    #[arg(short = 'F', long, default_value_t = 1540)]
    excludeflag: u16,

    /// Threads
    #[arg(short, long, default_value_t = 8)]
    thread: usize,

    /// Minimal len to define an large indels
    #[arg(short, long, default_value_t = 30)]
    indelMin: i64,

    /// Output file name
    #[arg(short, long)]
    out: PathBuf,
}


fn main() {
    let cli = Cli::parse();
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
    let mut alignment_vec: Vec<AlignmentPos> = Vec::new();

    while let Some(result) = bam.read(&mut record) {
        match result {
            Err(_) => break,
            Ok(()) => {}
        }
        if record.is_secondary()
            || record.is_unmapped()
            || record.mapq() < cli.mapq
            || (record.flags() & cli.excludeflag) != 0
        // == 0 means not match with flag.
        {
            continue;
        }
        let st = record.strand().to_owned(); // MUST move out of match, Because of mutable borrow by strand().
                                             // Start to extact SA signals
        match record.aux("SA".as_bytes()) {
            Ok(_sa) => {
                /* put the preliminary alignment to  */
                // dbg!(&record.cigar());
                let contig_name = record.contig();
                let pos = record.pos();
                let strand = match st {
                    ReqStrand::Forward => 1,
                    ReqStrand::Reverse => -1,
                };
                let mapq = record.mapq();
                // let cigar_stats_nuc = record.cigar_stats_nucleotides();
                let mut cigar_map = HashMap::from([
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
                let mut first_cigar_str = "".to_string();
                // dbg!(&record.cigar());
                record.cigar().iter().for_each(|cigar| match cigar {
                    &Cigar::Del(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"D".to_string();
                        cigar_map.insert('D', n.to_owned());
                    }
                    &Cigar::Match(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"M".to_string();
                        cigar_map.insert('M', n.to_owned());
                    }
                    &Cigar::Ins(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"I".to_string();
                        cigar_map.insert('I', n.to_owned());
                    }
                    &Cigar::RefSkip(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"N".to_string();
                        cigar_map.insert('N', n.to_owned());
                    }
                    &Cigar::SoftClip(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"S".to_string();
                        cigar_map.insert('S', n.to_owned());
                    }
                    &Cigar::HardClip(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"H".to_string();
                        cigar_map.insert('H', n.to_owned());
                    }
                    &Cigar::Pad(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"P".to_string();
                        cigar_map.insert('P', n.to_owned());
                    }
                    &Cigar::Equal(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"=".to_string();
                        cigar_map.insert('=', n.to_owned());
                    }
                    &Cigar::Diff(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"X".to_string();
                        cigar_map.insert('X', n.to_owned());
                    }
                });
                // dbg!(&contig_name,pos, std::str::from_utf8(record.qname()).unwrap(), &cigar_map);
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
                    // dbg!(&sa_list_clean);
                    for single_sa in sa_list_clean {
                        alignment_vec.push(parse_supplementary_alignment(single_sa));
                    }
                }
                alignment_vec.sort_by(|a, b| splitter_order_cmp(a, b));
                for i in 1..alignment_vec.len() {
                    let j = i - 1;
                    let a: &AlignmentPos = &alignment_vec[j];
                    let b: &AlignmentPos = &alignment_vec[i];
                    let bed_line: String;
                    // dbg!(&a);
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
                }
            }
            Err(_) => {}
        };
        alignment_vec.clear();

        // parse the CIGAR value for each record
        let cigar = record.cigar();
        // dbg!(&cigar.pos());
        let _ = &cigar.iter().for_each(|x| {
            // dbg!(&x,&x.len(),&x.char());
            // match x.char() {
            //     'D' =>{
            //         if x.len() >= cli.indelMin {
            //             // dbg!("catch them".to_string());
            //             // giggle 其实接受的只是bed文件所以只需要report这个起始和终止就可以了。
            //             let cigar_str = format!("{}",cigar);
            //             let strand = match st {
            //                 ReqStrand::Forward => 1,
            //                 ReqStrand::Reverse => -1,
            //             };

            //             let cigar_map = parse_cigar(&cigar_str);
            //             dbg!(&cigar_map);
            //             let alignment_pos = AlignmentPos::new(
            //                 &record.contig(),
            //                 &record.pos(),
            //                 cigar_map,
            //                 &strand,
            //                 &record.mapq(),
            //                 &cigar_str
            //             );
            //             dbg!(alignment_pos);
            //         }
            //     }
            //     'I' =>{
            //         if x.len() >= cli.indelMin {
            //             let cigar_str = format!("{}",cigar);
            //             let strand = match st {
            //                 ReqStrand::Forward => 1,
            //                 ReqStrand::Reverse => -1,
            //             };

            //             let cigar_map = parse_cigar(&cigar_str);
            //             dbg!(&cigar_map);
            //             let alignment_pos = AlignmentPos::new(
            //                 &record.contig(),
            //                 &record.pos(),
            //                 cigar_map,
            //                 &strand,
            //                 &record.mapq(),
            //                 &cigar_str
            //             );
            //         }
            //     },
            //     _=>{

            //     }
            // }
            // for cigar in record.cigar().iter(){
            //     if cigar == &Cigar::Del(n){

            //     }
            // }
            let mut pos = record.pos() as i64;
            let strand = match record.strand() {
                ReqStrand::Forward => 1i64,
                ReqStrand::Reverse => -1i64,
            };
            let chrom_clean: String;
            if record.contig().starts_with(&"chr".to_string()) {
                chrom_clean = record.contig().strip_prefix(&"chr".to_string()).unwrap().to_string();
            } else {
                chrom_clean = record.contig().to_string();
            }
            record.cigar().iter().for_each(|cigar| match cigar {
                &Cigar::Del(n) => {
                    let m = n as i64;
                    if m >= cli.indelMin {
                        let pt = format!("{}\t{}\t{}\t{}\n", chrom_clean, pos, (pos + m), strand);
                        f.write(pt.as_bytes()).unwrap();
                        pos += strand * m; //finally add the del length;
                    } else {
                        pos += strand * m;
                    }
                }
                &Cigar::Ins(n) => {
                    let m = n as i64;
                    if m >= cli.indelMin {
                        let pt = format!("{}\t{}\t{}\t{}\n", chrom_clean, pos, (pos + m), strand);
                        f.write(pt.as_bytes()).unwrap();
                        pos += strand * m; //finally add the del length;
                    } else {
                        pos += strand * m;
                    }
                }
                _ => {
                    pos += strand * (cigar.len() as i64);
                }
            });
        });

        // dbg!("\n\n\n".to_string());
    }
}

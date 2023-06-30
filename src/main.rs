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
    path::PathBuf,
    process::exit,
};
mod split_read_event;
use split_read_event::SplitReadEvent;
mod utils;
use utils::*;
mod aligments_event;
use aligments_event::*;
#[derive(Parser, Debug)]
#[command(name = "excord-LR")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.2")]
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

    /// Minimal len to define an SV events in CIGAR
    #[arg(short, long, default_value_t = 50)]
    indel_min: u32,

    /// Threshold to merge two adjacent events
    #[arg(short, long, default_value_t = 50)]
    merge_min: u32,

    /// Output file name
    #[arg(short, long)]
    out: PathBuf,

    /// Only report split-read event
    #[arg(short, long, default_value_t = false)]
    split_only: bool,
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
    let mut alignment_vec: Vec<SplitReadEvent> = Vec::new();

    while let Some(result) = bam.read(&mut record) {
        match result {
            Err(_) => break, //exit if the last one was processed.
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
        let contig_name = record.contig();
        let pos = record.pos();
        let strand = match st {
            ReqStrand::Forward => 1,
            ReqStrand::Reverse => -1,
        };
        match record.aux("SA".as_bytes()) {
            Ok(_sa) => {
                /* put the preliminary alignment to  */
                // dbg!(&record.cigar());

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
                // dbg!("xxxxxxxxxxxxxxxx".to_string());
                record.cigar().iter().for_each(|cigar| match cigar {
                    &Cigar::Del(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"D".to_string();
                        cigar_map.entry('D').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::Match(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"M".to_string();
                        cigar_map.entry('M').and_modify(|e| *e += n.to_owned());
                        // dbg!(&n,&cigar_map);
                    }
                    &Cigar::Ins(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"I".to_string();
                        cigar_map.entry('I').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::RefSkip(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"N".to_string();
                        cigar_map.entry('N').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::SoftClip(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"S".to_string();
                        cigar_map.entry('S').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::HardClip(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"H".to_string();
                        cigar_map.entry('H').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::Pad(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"P".to_string();
                        cigar_map.entry('P').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::Equal(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"=".to_string();
                        cigar_map.entry('=').and_modify(|e| *e += n.to_owned());
                    }
                    &Cigar::Diff(n) => {
                        first_cigar_str += &n.to_owned().to_string();
                        first_cigar_str += &"X".to_string();
                        cigar_map.entry('X').and_modify(|e| *e += n.to_owned());
                    }
                });
                // dbg!(&cigar_map,&first_cigar_str);
                // dbg!(&contig_name,pos, std::str::from_utf8(record.qname()).unwrap(), &cigar_map);
                alignment_vec.push(SplitReadEvent::new(
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
                    let a: &SplitReadEvent = &alignment_vec[j];
                    let b: &SplitReadEvent = &alignment_vec[i];
                    let bed_line: String;
                    // dbg!(&a);
                    if alignment_pos_cmp(a, b) == Ordering::Less {
                        bed_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            a.chrom,
                            a.start,
                            a.end,
                            a.strand,
                            b.chrom,
                            b.start,
                            b.end,
                            b.strand,
                            alignment_vec.len() - 1,
                            "excord-lr-split-read",
                            std::str::from_utf8(record.qname()).unwrap(),
                            &strand,
                            record.flags()
                        );
                    } else {
                        bed_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            b.chrom,
                            b.start,
                            b.end,
                            b.strand,
                            a.chrom,
                            a.start,
                            a.end,
                            a.strand,
                            alignment_vec.len() - 1,
                            "excord-lr-split-read",
                            std::str::from_utf8(record.qname()).unwrap(),
                            &strand,
                            record.flags()
                        );
                    }
                    f.write(bed_line.as_bytes()).unwrap();
                }
            }
            Err(_) => {}
        };
        alignment_vec.clear();
        if cli.split_only == false {
            // parse the CIGAR value for each record
            let cigar = record.cigar();
            // dbg!(&cigar.to_string());
            // let cigar_string = cigar.to_string();

            // let alignment_event_vec:Vec<AlignmentPos> = parse_alignment_event(&cigar_string);
            let mut total_consume = 0u32;
            for x in cigar.into_iter() {
                match x {
                    &Cigar::Del(n) => {
                        total_consume += n;
                    }
                    &Cigar::Match(n) => {
                        total_consume += n;
                    }
                    &Cigar::RefSkip(n) => {
                        total_consume += n;
                    }
                    &Cigar::Equal(n) => {
                        total_consume += n;
                    }
                    _ => {}
                }
            }

            let mut aligments_event_vec: Vec<AlignmentEvent> = vec![];
            let mut left_consume = 0u32;
            let mut right_consume = total_consume;
            for x in cigar.into_iter() {
                match x {
                    &Cigar::Del(n) => {
                        right_consume -= n;

                        if n > cli.indel_min {
                            // report one event

                            let aligments_event = AlignmentEvent::new(
                                &contig_name,
                                &left_consume,
                                &right_consume,
                                &n,
                                &pos,
                                &strand,
                            );
                            aligments_event_vec.push(aligments_event);
                            // dbg!(
                            //     "xxx".to_string(),
                            //     &record.pos(),
                            //     &cigar.len(),
                            //     &cigar.to_string(),
                            //     &left_consume,
                            //     &right_consume,
                            //     &n,
                            //     aligments_event.get_bed_record()
                            // );
                        }

                        left_consume += n;
                    }
                    &Cigar::Match(n) => {
                        left_consume += n;
                        right_consume -= n;
                    }
                    &Cigar::RefSkip(n) => {
                        left_consume += n;
                        right_consume -= n;
                    }
                    &Cigar::Equal(n) => {
                        left_consume += n;
                        right_consume -= n;
                    }
                    _ => {}
                }
            }

            // do merge step
            
            // print all record
            aligments_event_vec.into_iter().for_each(|x|{
                f.write(x.get_bed_record().as_bytes()).unwrap();
            });

            dbg!(&total_consume);
        }
    }
}

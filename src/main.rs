use bio_types::{genome::AbstractInterval, strand::ReqStrand};
use clap::Parser;
use rust_htslib::{
    bam,
    bam::{
        // ext::BamRecordExtensions,
        record::{Aux, Cigar},
        Read,
        Record
    },
};
use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    ops::BitAnd,
    path::PathBuf,
    process::exit,
};
mod split_read_event;
use split_read_event::SplitReadEvent;
mod utils;
use utils::*;
mod aligments_event;
use aligments_event::*;
use std::str;

#[derive(Parser, Debug)]
#[command(name = "excord-LR")]
#[command(author = "Xinchang Zheng <zhengxc93@gmail.com>")]
#[command(version = "0.1.17")]
#[command(about = "
Extract Structural Variation Signals from Long-Read BAMs
Contact: Xinchang Zheng <zhengxc93@gmail.com,Xinchang.Zheng@bcm.edu>
", long_about = None)]
struct Cli {
    /// Path to BAM file
    #[arg(short, long)]
    bam: PathBuf,

    /// Path to reference, used for CRAM file
    #[arg(short, long)]
    reference:Option<PathBuf>,

    /// Minimal MapQ
    #[arg(short = 'Q', long, default_value_t = 1)]
    mapq: u8,

    /// Exclude Flags
    #[arg(short = 'F', long, default_value_t = 1796)]
    exclude_flag: u16,

    /// Exclude Secondary Alignment
    #[arg(short = 'S', long, default_value_t = false)]
    exclude_secondary: bool,

    /// Exclude Unmapped Alignment
    #[arg(short = 'U', long, default_value_t = false)]
    exclude_unmapped: bool,

    /// Threads
    #[arg(short, long, default_value_t = 8)]
    thread: usize,

    /// Minimal length to define an SV events in CIGAR
    #[arg(short, long, default_value_t = 50)]
    indel_min: u32,

    /// Threshold to merge two adjacent events
    #[arg(short, long, default_value_t = 5)]
    merge_min: u32,

    /// Minimal length of hard-clip and soft-clip to define a large insertion signal
    #[arg(short, long, default_value_t = 1000)]
    ins_clip_min: u32,

    /// Not merge
    #[arg(short, long, default_value_t = false)]
    not_merge: bool,

    /// Output file name
    #[arg(short, long)]
    out: PathBuf,

    /// Only report split-read event
    #[arg(short, long, default_value_t = false)]
    split_only: bool,

    /// Percent of overlap to discard a potential false positive record[Optional]
    #[arg(short = 'p', long, default_value_t = 0.0)]
    max_pct_overlap: f64,

    /// Maximal number of SA to include a record
    #[arg(short = 'k', long, default_value_t = 4)]
    max_supp_alignm: usize,

    /// Debug
    #[arg(short, long, default_value_t = false)]
    debug: bool,

    /// Verbose output
    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

fn main() {
    let cli = Cli::parse();
    if cli.debug {
        println!("{:?}", &cli);
    }

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

    let mut bam = bam::Reader::from_path(&_bam).unwrap();



    if _bam.to_str().unwrap().to_lowercase().contains(&".cram".to_string()) {
        match cli.reference {
            Some(reference)=>{
                bam.set_reference(reference).unwrap();
            },
            None =>{
                println!("excord-lr is running on CRAM file, reference(-r) is required.");
                exit(1)
            }
        }
    }

    

    bam.set_threads(cli.thread).unwrap();
    let mut record = Record::new();

    while let Some(result) = bam.read(&mut record) {
        // if str::from_utf8(&record.qname()).unwrap() == "4aa7ca6c-d970-4b23-90ec-6e449956685f".to_string() {
        //     dbg!(&record);
        // }

        // dbg!(str::from_utf8(&record.qname()).unwrap());

        match result {
            Err(_) => break, //exit if the last one was processed.
            Ok(()) => {}
        }
        if cli.exclude_secondary & record.is_secondary() {
            // dbg!("found secondary", record.qname());
            continue;
        }

        if cli.exclude_unmapped & record.is_unmapped() {
            // dbg!("found unmapped", record.qname());
            continue;
        }

        if record.mapq() < cli.mapq {
            // dbg!("found low mapq", record.qname());
            continue;
        }

        // dbg!(272u16.bitand(256u16));
        if record.flags().bitand(cli.exclude_flag) != 0u16
        // == 0 means not match with flag.
        {
            dbg!("found filtered flag", record.qname());
            continue;
        }

        // if str::from_utf8(&record.qname()).unwrap() == "4aa7ca6c-d970-4b23-90ec-6e449956685f".to_string() {
        //     dbg!(&record);
        // }

        let st = record.strand().to_owned(); // MUST move out of match, Because of mutable borrow by strand().
                                             // Start to extact SA signals
        let contig_name = record.contig();
        let pos = record.pos();
        let strand = match st {
            ReqStrand::Forward => 1,
            ReqStrand::Reverse => -1,
        };
        let mut alignments_event_vec: Vec<AlignmentEvent> = vec![];

        match record.aux("SA".as_bytes()) {
            Ok(_sa) => {
                /* put the preliminary alignment to  */
                // dbg!(&record.cigar());

                let mut alignment_vec: Vec<SplitReadEvent> = Vec::new();
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
                // let mut left_S:u32 = 0;
                // let mut right_S:u32 = 0;
                // let mut is_have_left_S:bool = false;
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
                        // if is_have_left_S == false{
                        //     is_have_left_S = true;
                        //     left_S = n.to_owned();
                        // }else{
                        //     right_S = n.to_owned();
                        // }
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

                    if sa_list.len() > cli.max_supp_alignm {
                        continue;
                    }

                    let sa_list_clean = sa_list.iter().filter(|x| x.len() > 0);
                    // dbg!(&sa_list_clean);
                    for single_sa in sa_list_clean {
                        alignment_vec.push(parse_supplementary_alignment(single_sa));
                    }
                }
                // dbg!(&std::str::from_utf8(record.qname()),&alignment_vec);
                alignment_vec.sort_by(|a, b| splitter_order_cmp(a, b));

                // iterally remove the potential FP that caused by the secondary alignment which is very close to
                // the primary alignment.
                // e.g.,  primary alignment reported region 1	2367366	2368042; the secondary alignment is 1	2367371	2368042
                // this will introduce a false positive
                //
                // In order to remove it. next step I will compare each adjacent alignment pairs. if the overlap of the region
                // is very long. one of them will be removed.
                //

                // process the large insertion
                //
                // Rules:
                // 1. Only primary and one supplimentary alignment
                // 2. Must have overlap
                // 3. Both alignment should have a soft-clip more than 1kp.
                // 4. Both alignment have same chormosome
                if alignment_vec.len() == 2 {
                    let a = alignment_vec.first().clone().unwrap();
                    let b = alignment_vec.last().clone().unwrap();

                    //2. Must have overlap

                    // 3. Both alignment should have a soft-clip more than 1kp.
                    if *a.cigar_map.get(&'S').unwrap() > cli.ins_clip_min || *a.cigar_map.get(&'H').unwrap() > cli.ins_clip_min {
                        if a.chrom == b.chrom {
                            if a.strand == b.strand {
                                // if overlap > max_pct_overlap. default value of max_pct_overlap is 0 ,means each overlap pass the examination.
                                if overlap(&a.start, &a.end, &b.start, &b.end, cli.max_pct_overlap)
                                {
                                    if *b.cigar_map.get(&'S').unwrap() > cli.ins_clip_min || *b.cigar_map.get(&'H').unwrap() > cli.ins_clip_min {
                                        // will generate a new alignment event
                                        

                                        /*
                                        Linking id (READNAME) = f2720f51-47ee-4dad-91d2-3f2aa26a66d0
                                        Haplotype = 2
                                        # alignments = 2
                                        Total span = 46,603bp
                                        Strands = ++
                                        chr2:9,832,511-9,877,850 (+) = 45,339BP @MAPQ=60 NM=2886 CLIPPING=32S ... 8269S
                                        chr2:9,877,576-9,879,114 (+) = 1,538BP @MAPQ=60 NM=191 CLIPPING=51421H ... 14H
                                         */

                                         // make sure the insertion presents in the middile of two segments
                                        let mut pos_list = vec![a.start.clone(),a.end.clone(),b.start.clone(),b.end.clone()];

                                        pos_list.sort();


                                        let x = AlignmentEvent {
                                            lchrom: a.chrom.clone(),
                                            lstart: pos_list[0].clone() as u32,
                                            lend: pos_list[1].clone() as u32,
                                            lstrand: a.strand,
                                            rchrom: b.chrom.clone(),
                                            rstart: pos_list[1].clone() as u32,
                                            rend: pos_list[1].clone() as u32, // TODO: should estimate the length of insertions and add it at here. [done]
                                            rstrand: b.strand,
                                            events_num: 1,
                                            svtype: AlignEventType::Ins,
                                        };

                                        // alignment_vec.clear();
                                        let rrr = get_alignment_event_record(
                                            &x,
                                            &cli.verbose,
                                            &record,
                                            &strand,
                                            &"excord-lr-alignment-event-large-ins-two-alignments",
                                        );

                                        f.write(rrr.as_bytes()).unwrap();


                                        let x = AlignmentEvent {
                                            lchrom: a.chrom.clone(),
                                            lstart: pos_list[0].clone() as u32,
                                            lend: pos_list[2].clone() as u32,
                                            lstrand: a.strand,
                                            rchrom: b.chrom.clone(),
                                            rstart: pos_list[2].clone() as u32,
                                            rend: pos_list[2].clone() as u32, // TODO: should estimate the length of insertions and add it at here. [done]
                                            rstrand: b.strand,
                                            events_num: 1,
                                            svtype: AlignEventType::Ins,
                                        };

                                        // alignment_vec.clear();
                                        let rrr = get_alignment_event_record(
                                            &x,
                                            &cli.verbose,
                                            &record,
                                            &strand,
                                            &"excord-lr-alignment-event-large-ins-two-alignments",
                                        );

                                        f.write(rrr.as_bytes()).unwrap();


                                    }
                                }
                            }
                        } else {
                            let x = AlignmentEvent {
                                lchrom: a.chrom.clone(),
                                lstart: a.start as u32,
                                lend: a.end as u32,
                                lstrand: a.strand,
                                rchrom: a.chrom.clone(),
                                rstart: a.end as u32,
                                rend: a.end as u32,
                                rstrand: a.strand,
                                events_num: 1,
                                svtype: AlignEventType::Ins,
                            };

                            let rrr = get_alignment_event_record(
                                &x,
                                &cli.verbose,
                                &record,
                                &strand,
                                &"excord-lr-alignment-event-large-ins-one-alignments",
                            );

                            f.write(rrr.as_bytes()).unwrap();
                        }
                    }
                }

                // process the super large insertion. e.g. a insertion larger than average length of long-read
                //
                // Rules:
                // 1. Only primary alignment
                // 3. Alignment should have a soft-clip more than 1kp.

                if alignment_vec.len() == 1 {
                    let a = alignment_vec.first().clone().unwrap();
                    if *a.cigar_map.get(&'S').unwrap() > cli.ins_clip_min || *a.cigar_map.get(&'H').unwrap() > cli.ins_clip_min {
                        // alignments_event_vec.push();
                        let x = AlignmentEvent {
                            lchrom: a.chrom.clone(),
                            lstart: a.start as u32,
                            lend: a.end as u32,
                            lstrand: a.strand,
                            rchrom: a.chrom.clone(),
                            rstart: a.end as u32,
                            rend: a.end as u32,
                            rstrand: a.strand,
                            events_num: 1,
                            svtype: AlignEventType::Ins,
                        };

                        let rrr = get_alignment_event_record(
                            &x,
                            &cli.verbose,
                            &record,
                            &strand,
                            &"excord-lr-alignment-event-large-ins",
                        );

                        f.write(rrr.as_bytes()).unwrap();
                    }
                }

                for i in 1..alignment_vec.len() {
                    let j = i - 1;
                    let a: &SplitReadEvent = &alignment_vec[j];
                    let b: &SplitReadEvent = &alignment_vec[i];
                    let bed_line: String;
                    // dbg!(&a);
                    if alignment_pos_cmp(a, b) == Ordering::Greater {
                        bed_line = get_alignment_split_record(
                            &b,
                            &a,
                            &cli.verbose,
                            &record,
                            &strand,
                            &alignment_vec.len(),
                            &"excord-lr-split-read",
                        );
                    } else {
                        bed_line = get_alignment_split_record(
                            &a,
                            &b,
                            &cli.verbose,
                            &record,
                            &strand,
                            &alignment_vec.len(),
                            &"excord-lr-split-read",
                        );
                    }
                    f.write(bed_line.as_bytes()).unwrap();
                }
            }
            Err(_) => {}
        };

        // Extract Alignment event if split_only option is disabled.

        if cli.split_only == false {
            // parse the CIGAR value for each record
            let cigar = record.cigar();

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

            let mut left_consume = 0u32;
            let mut right_consume = total_consume;
            for x in cigar.into_iter() {
                match x {
                    &Cigar::Del(n) => {
                        right_consume -= n;
                        if n >= cli.indel_min {
                            // report one event
                            let aligments_event = AlignmentEvent::new(
                                &contig_name,
                                &left_consume,
                                &right_consume,
                                &n,
                                &pos,
                                &strand,
                                Some(AlignEventType::Del),
                            );
                            alignments_event_vec.push(aligments_event);
                        }
                        left_consume += n;
                    }
                    &Cigar::Ins(n) => {
                        if n >= cli.indel_min {
                            let aligments_event = AlignmentEvent::new(
                                &contig_name,
                                &left_consume,
                                &n, // for Insetions derived from CIGAR value, encode the length in right region.
                                &(0u32),
                                &pos,
                                &strand,
                                Some(AlignEventType::Ins),
                            );
                            // dbg!("insertion",&aligments_event);
                            alignments_event_vec.push(aligments_event);
                        }
                        // do not recaculate the comsume, beacuse the INS does not take account for the reference.
                        // right_consume -= n;
                        // left_consume += n;
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

            // if str::from_utf8(&record.qname()).unwrap() == "4aa7ca6c-d970-4b23-90ec-6e449956685f".to_string() {
            //     dbg!(&alignments_event_vec);
            //     dbg!(&alignments_event_vec.len());

            // }

            // do merge step
            let mut merged_alignments_event_vec: Vec<AlignmentEvent> = vec![];

            // merge two alignments if needed
            if alignments_event_vec.len() == 2 {
                let a: AlignmentEvent = alignments_event_vec[0].clone();
                let b: AlignmentEvent = alignments_event_vec[1].clone();
                if b.lend.abs_diff(a.rstart) < cli.merge_min
                    && (a.svtype == AlignEventType::Del && b.svtype == AlignEventType::Del)
                {
                    let c: AlignmentEvent = AlignmentEvent {
                        lchrom: a.lchrom,
                        lstart: a.lstart,
                        lend: a.lend,
                        lstrand: a.lstrand,
                        rchrom: b.rchrom,
                        rstart: b.rstart,
                        rend: b.rend,
                        rstrand: b.rstrand,
                        events_num: 1, // hard code 1, STIX will parse it as a split-event.
                        svtype: AlignEventType::Del,
                    };

                    merged_alignments_event_vec.push(c);
                } else {
                    merged_alignments_event_vec.push(a);
                    merged_alignments_event_vec.push(b);
                }
            } else if alignments_event_vec.len() > 2 {
                let mut merge1: Vec<AlignmentEvent> = alignments_event_vec.clone();
                let mut merge2: Vec<AlignmentEvent> = vec![];
                let mut iter_times = 5u32;

                // try to merge several times to handle more than two adjecent dels.
                // 1. merge 5 times, each time will based on the previous one.
                // 2. If the original vec is length as n, we iterally compare previous and next for the 2 to n-1 elements.
                // 3. Only works in Del

                // each round
                //  [ a, b, c, d, e]
                //    |--|--| -> idx =1
                //       |--|--| --> idx =2
                //          |--|--| --> idx =3 

                // 5 round in total
                // [ a, b, c, d, e] --> init
                // [ ab, c,d,e] --> round 1
                // [ abc, d, e] --> round 2
                // ...
                // round 5 ----> merged_alignment_vec and go down.

                loop {
                    iter_times -= 1;
                    let mut idx = 1usize;

                    loop {
                        if idx > (&merge1.len() - 2) {
                            break;
                        }
                        let pidx = idx - 1;
                        let nidx = idx + 1;
                        let prv: AlignmentEvent = merge1[pidx].clone();
                        let target: AlignmentEvent = merge1[idx].clone();
                        let nxt: AlignmentEvent = merge1[nidx].clone();

                        if (prv.lend.abs_diff(target.rstart) < cli.merge_min
                            && (target.svtype == AlignEventType::Del
                                && prv.svtype == AlignEventType::Del))
                            || (target.lend.abs_diff(nxt.rstart) < cli.merge_min
                                && (target.svtype == AlignEventType::Del
                                    && nxt.svtype == AlignEventType::Del))
                        {
                            if prv.lend.abs_diff(target.rstart) < cli.merge_min
                                && (target.svtype == AlignEventType::Del
                                    && prv.svtype == AlignEventType::Del)
                            {
                                let c: AlignmentEvent = AlignmentEvent {
                                    lchrom: prv.lchrom,
                                    lstart: prv.lstart,
                                    lend: prv.lend,
                                    lstrand: prv.lstrand,
                                    rchrom: target.rchrom,
                                    rstart: target.rstart,
                                    rend: target.rend,
                                    rstrand: target.rstrand,
                                    events_num: 1, // hard code 1, STIX will parse it as a split-event.
                                    svtype: AlignEventType::Del,
                                };
                                // merge2.push(c);
                                merge2.push(c);
                            }

                            if target.lend.abs_diff(nxt.rstart) < cli.merge_min
                                && (target.svtype == AlignEventType::Del
                                    && nxt.svtype == AlignEventType::Del)
                            {
                                let c: AlignmentEvent = AlignmentEvent {
                                    lchrom: target.lchrom,
                                    lstart: target.lstart,
                                    lend: target.lend,
                                    lstrand: target.lstrand,
                                    rchrom: nxt.rchrom,
                                    rstart: nxt.rstart,
                                    rend: nxt.rend,
                                    rstrand: nxt.rstrand,
                                    events_num: 1, // hard code 1, STIX will parse it as a split-event.
                                    svtype: AlignEventType::Del,
                                };
                                // merge2.push(c);
                                merge2.push(c);
                            }
                            idx += 1;
                        } else {
                            if idx == 1 {
                                merge2.push(prv);
                                merge2.push(target);
                                merge2.push(nxt);
                            } else {
                                merge2.push(nxt);
                            }

                            idx += 1;
                        }
                    }

                    if merge1.len() == merge2.len() {
                        break;
                    } else if iter_times <= 0 {
                        break;
                    } else {
                        merge1 = merge2.clone();
                        merge2.clear();
                    }
                }
                merged_alignments_event_vec = merge2.clone();

                // if str::from_utf8(&record.qname()).unwrap()
                //     == "4aa7ca6c-d970-4b23-90ec-6e449956685f".to_string()
                // {
                //     dbg!(&merged_alignments_event_vec);
                //     dbg!(&merged_alignments_event_vec.len());
                // }

                // merged_alignments_event_vec.dedup_by(same_bucket)
            } else if alignments_event_vec.len() == 1 {
                let a: AlignmentEvent = alignments_event_vec[0usize].clone();
                merged_alignments_event_vec.push(a);
            }

            merged_alignments_event_vec.into_iter().for_each(|x| {
                let rrr = get_alignment_event_record(
                    &x,
                    &cli.verbose,
                    &record,
                    &strand,
                    &"excord-lr-alignment-event",
                );

                f.write(rrr.as_bytes()).unwrap();
            });
        }
        alignments_event_vec.clear();
    }
}

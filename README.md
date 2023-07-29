# Introduction

Extract Structural Variation signals in Long Reads BAM. This program is a Rust re-implementation of [excord](https://github.com/brentp/excord) which is used to generate indexes for [STIX](https://github.com/ryanlayer/stix). Currently, this program only considers split reads.

# How to use

```
Contact: Xinchang Zheng <zhengxc93@gmail.com,Xinchang.Zheng@bcm.edu>


Usage: excord-lr [OPTIONS] --bam <BAM> --out <OUT>

Options:
  -b, --bam <BAM>                    Path to BAM file
  -Q, --mapq <MAPQ>                  Minimal MapQ [default: 1]
  -F, --exclude-flag <EXCLUDE_FLAG>  Exclude Flags [default: 1540]
  -S, --exclude-secondary            Exclude Secondary Alignment
  -U, --exclude-unmapped             Exclude Unmapped Alignment
  -t, --thread <THREAD>              Threads [default: 8]
  -i, --indel-min <INDEL_MIN>        Minimal len to define an SV events in CIGAR [default: 50]
  -m, --merge-min <MERGE_MIN>        Threshold to merge two adjacent events [default: 5]
  -n, --not-merge                    Not merge
  -o, --out <OUT>                    Output file name
  -s, --split-only                   Only report split-read event
  -d, --debug                        debug
  -h, --help                         Print help
  -V, --version                      Print version
```

### TODO
1. merge more than two adjacent SVs
2. add inverted bed region rather than single one
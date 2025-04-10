# Introduction

Extract Structural Variation signals in Long Reads BAM. This program is a Rust re-implementation of [excord](https://github.com/brentp/excord) which is used to generate indexes for [STIX](https://github.com/ryanlayer/stix). Currently, this program considiers Split-Event and Alignment-Event.

# How to use

```
Contact: Xinchang Zheng <zhengxc93@gmail.com,Xinchang.Zheng@bcm.edu>


Usage: excord-lr [OPTIONS] --bam <BAM> --out <OUT>

Options:
  -b, --bam <BAM>
          Path to BAM file
  -Q, --mapq <MAPQ>
          Minimal MapQ [default: 1]
  -F, --exclude-flag <EXCLUDE_FLAG>
          Exclude Flags [default: 1796]
  -S, --exclude-secondary
          Exclude Secondary Alignment
  -U, --exclude-unmapped
          Exclude Unmapped Alignment
  -t, --thread <THREAD>
          Threads [default: 8]
  -i, --indel-min <INDEL_MIN>
          Minimal len to define an SV events in CIGAR [default: 50]
  -m, --merge-min <MERGE_MIN>
          Threshold to merge two adjacent events [default: 5]
  -n, --not-merge
          Not merge
  -o, --out <OUT>
          Output file name
  -s, --split-only
          Only report split-read event
  -p, --pct-overlap <PCT_OVERLAP>
          percent of overlap to discard a potential false positive record(set 0 to disable) [default: 0.8]
  -k, --max-supp-alignm <MAX_SUPP_ALIGNM>
          maximal number of SA to include a record [default: 4]
  -d, --debug
          debug
  -v, --verbose
          verbose output
  -h, --help
          Print help
  -V, --version
          Print version
```


# Install

Please check the release page to find the latest version.

# Build

We recommend users build excord-lr with [musl libc](https://musl.libc.org/). The binary works for the most of linux distributions.


```
sudo apt install musl-tools
sudo apt install build-essential
sudo apt install clang libclang-dev

rustup target add x86_64-unknown-linux-musl

CC=/usr/bin/musl-gcc  cargo build --release --target=x86_64-unknown-linux-musl
```

You can also build excord-lr in a normal way by this:

```
 cargo build --release
```


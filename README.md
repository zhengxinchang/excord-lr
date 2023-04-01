# Introduction

Extract Structural Variation signals in Long Reads BAM. This program is a Rust re-implementation of [excord](https://github.com/brentp/excord) which is used to generate indexes for [STIX](https://github.com/ryanlayer/stix). Currently, this program only considers split reads.

# How to use

```
Usage: excord-lr [OPTIONS] --bam <BAM> --out <OUT>

Options:
  -b, --bam <BAM>                  Path to BAM file
  -F, --mapq <MAPQ>                Minimal MapQ [default: 1]
  -Q, --excludeflag <EXCLUDEFLAG>  Exclude Flags [default: 1540]
  -o, --out <OUT>                  Output file name
  -h, --help                       Print help
  -V, --version                    Print version
```


# Performance

**Hardware**

* CPU: AMD Ryzen 7 3700X 8-Core Processor
* RAM: 16GB

|Sample|Reads Number|excord(second)|excord-lr(second)|
| --: | --: | --: | --: |
| Sample1 |21526088 |505.24 |127.67|
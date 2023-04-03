# Introduction

Extract Structural Variation signals in Long Reads BAM. This program is a Rust re-implementation of [excord](https://github.com/brentp/excord) which is used to generate indexes for [STIX](https://github.com/ryanlayer/stix). Currently, this program only considers split reads.

# How to use

```
Usage: excord-lr [OPTIONS] --bam <BAM> --out <OUT>

Options:
  -b, --bam <BAM>                  Path to BAM file
  -Q, --mapq <MAPQ>                Minimal MapQ [default: 1]
  -F, --excludeflag <EXCLUDEFLAG>  Exclude Flags [default: 1540]
  -t, --thread <THREAD>            Threads [default: 8]
  -o, --out <OUT>                  Output file name
  -h, --help                       Print help
  -V, --version                    Print version
```


# Performance

**Test1**

* CPU: AMD Ryzen 7 3700X 8-Core Processor
* RAM: 16GB

|Sample|Reads Number|excord(second)|excord-lr(second)|
| --: | --: | --: | --: |
| Sample1 |21526088 |505.24 |127.67|


**Test2**

* CPU: AMD Ryzen 7 PRO 6850HS with Radeon Graphics
* RAM: 32GB

|Sample|Reads Number|excord(second)|excord-lr(second)|
| --: | --: | --: | --: |
| Sample1 |21526088 |548.99 |73.073|
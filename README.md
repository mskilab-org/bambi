# bxBam
indexed file format for barcoded BAMs with API for converting and accessing alignment records

## Table of contents
* [Installation](#installation)
* [Examples to run in command line](#examples-to-run-from-the-command-line)
  * [Index Generation](#index-generation)
  * [Query Bam File and Return Headerless .sam](#query-bam-file-and-return-headerless-sam)
  * [Query Bam File and Return Headered .bam](#query-bam-file-and-return-headered-bam)
* [Examples to run in R](#examples-to-run-in-r)

## Installation

```R
## install packages via devtools::install_github().
library(devtools)
install_github("mskilab/bxBam")
```

## Examples to run from the command line
```bash
git clone git@github.com:mskilab/bxBam.git
make
## LD_LIBRARY_PATH may have to be set to point to lmdb, htslib, and ck.
```

Index Generation
----------------
```bash
bxbam -t "lmdb" -f path_to_bam_file
## for example:
bxbam -t "lmdb" -f HCC1143_BL_phased_possorted.bam
```

Query Bam File and Return Headerless sam
-----------------------------------------
```bash
bxbam bam_file query_string
## for example:
bxbam HCC1143_BL_phased_possorted.bam GTGGTCGCAACGCTTA-1
```

Query Bam File and Return Headered bam
--------------------------------------
```bash
bxbam -h bam_file query_string
## for example:
bxbam -h HCC1143_BL_phased_possorted.bam GTGGTCGCAACGCTTA-1
## the results will be in a file called headered.bam
```

Examples to run in R
--------------------

```R
## via functions available in package.
library(bxBam)
ls.str('package:bxBam')
```

```R
## to generate db file specify bam file generated from 10x Genomics. Use generate_bxi()
generate_bxi(system.file("extdata", "shortened_HCC1143_BL_phased_possorted.bam", package="bxBam")
```

```R
## using that index file, query the bam file for records that match a given BX barcode.
## barcodedReads(bam_file, index_file, barcode)

barcodedReads(system.file("extdata", "shortened_HCC1143_BL_phased_possorted.bam", package="bxBam"),
system.file("extdata", "shortened_HCC1143_BL_phased_possorted.db", package="bxBam"), 'CGGAGCTAGTAAGTAC-1')
```

# bxBam
indexed file format for barcoded BAMs with API for converting and accessing alignment records

## Table of contents
* [Installation](#installation)
* [Examples](#examples)

## Installation

```R
## install packages via devtools::install_github().
library(devtools)
install_github("mskilab/bxBam", auth_token="your_auth_token")
```

## Examples

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

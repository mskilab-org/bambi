# bambi
R package for querying 10x WGS and single-cell BAMs

## Dependencies


```{r}
devtools::install_github('mskilab/gUtils')
```

```{r}
devtools::install_github('mskilab/bamUtils')
```


cbambi
https://github.com/mskilab/cbambi

## bambi commands


Instantiate a bambi object:

Methods:

* grab_bx()
    * `grab_bx(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)`

* grab_cb()
    * `grab_cb(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)`

* grab_ub()
    * `grab_ub(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)`

* fetch_by_tag()
    * `fetch_by_tag(tag, tag_queries, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1)`

* fetch_by_multiple_tag()? 
    * remove
    * the only downside is that `fetch_by_tag()` might take a performance hit....

## Demo

* Instantiate a `bambi` object

```{r}
library(bambi)

> hcc1143_subset = bambi$new(bam_file = "subsetHCC1143_phased_possorted0001.bam", bamdb_path="subsetHCC1143_phased_possorted0001_lmdb")
```

* Call methods

```{r}
> hcc1143_subset$grab_bx('CGACGTGTCCTCTAGC-1')
GRanges object with 2 ranges and 11 metadata columns:
      seqnames                 ranges strand |
         <Rle>              <IRanges>  <Rle> |
  [1]     chr1 [147975454, 147975580]      + |
  [2]     chr1 [147975675, 147975824]      - |
                                         qname      flag      mapq       cigar
                                   <character> <numeric> <numeric> <character>
  [1] ST-K00126:3:H5TL3BBXX:2:2109:25926:37800        99        16        127M
  [2] ST-K00126:3:H5TL3BBXX:2:2109:25926:37800       147        16        150M
            rnext     pnext      tlen
      <character> <numeric> <numeric>
  [1]           = 147975676       371
  [2]           = 147975455      -371
                                                                                                                                                         seq
                                                                                                                                                 <character>
  [1]                        ATGTCTTCTTCCTCATTATCTGGCACTGGTTAGGAAGCACTCATCTCCATGAAGTCATCTTTTGTTAATTCCTCTGGTGTGGTGTGTATTAGCTCTTAAATTCCTCCAAGATCCATATCTTGCAACC
  [2] ATCTGGACACAAATTGTACTTTTGTCCAGCACGAATTTATTGTTTTGAGTTTCATGGTTTTCTATATCAACTGATGACATCTTGAAAGGTGTAAGCCTTCCAGACTTCCATGATGTTCTCTCTATTGGGTTTCTCTTTTGCAATGTTGAC
                                                                                                                                                        qual
                                                                                                                                                 <character>
  [1]                        JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJAJFJJJJJJJJJFJJJJJJJJJJFJJJJFFFJJJFJJJJJJAAJFJJJFAFAFFFJAA<7F<
  [2] A<7FFFJFFFAJJAAAJJF<F<7A-<AA-<<<AFFJJJJJJJJFFJAFFAAFJFJJJAFFJJJJJJJJJJFJFAJJJJJJFJJJJJJ<FFJJJFJJJFJJJJJJJJJJJJJFJJJJFFJ7JJJJF<JJJJJJJJJJJJJJJJJJJFFAA<
                      BX    qwidth
             <character> <integer>
  [1] CGACGTGTCCTCTAGC-1       127
  [2] CGACGTGTCCTCTAGC-1       150
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
  ```





repos:
https://github.com/madler/zlib
https://github.com/concurrencykit/ck 
https://github.com/LMDB/lmdb
https://github.com/samtools/htslib





library(testthat)
library(bambi)
library(bamUtils)
library(gUtils)


example_bam = 'subsetHCC1143_phased_possorted0001.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_bai = 'subsetHCC1143_phased_possorted0001.bam.bai' 
example_lmdb = 'subsetHCC1143_phased_possorted0001_lmdb'


######## utils.R


test_that('countCigar', {
 
    empty_string = 'foobar'
    expect_warning(countCigar(empty_string))   ### should be warning message: 'Warning message: In countCigar(example_bam) : NAs introduced by coercion'
    expect_equal(dim(countCigar(empty_string))[1], 1)
    expect_equal(dim(countCigar(empty_string))[2], 4)
    expect_match(colnames(countCigar(empty_string))[1], "D")
    expect_match(colnames(countCigar(empty_string))[2], "I")
    expect_match(colnames(countCigar(empty_string))[3], "M")
    expect_match(colnames(countCigar(empty_string))[4], "S")
    expect_equal(countCigar(empty_string)[1], 0)
    expect_equal(countCigar(empty_string)[2], 0)
    expect_equal(countCigar(empty_string)[3], 0)
    expect_equal(countCigar(empty_string)[4], 0)

})



test_that('bamflag', {

    gr_example = read.bam(example_bam, all=TRUE)[[42]]  ## random GRanges
    exdt = as.data.table(gr2dt(gr_example))
    ## isPaired  1 1
    expect_equal(as.data.frame(bamflag(exdt))$isPaired[1], 1)
    expect_equal(as.data.frame(bamflag(exdt))$isPaired[2], 1)
    expect_equal(as.data.frame(bamflag(gr_example))$isPaired[1], 1)
    expect_equal(as.data.frame(bamflag(gr_example))$isPaired[2], 1)
    ## isProperPair  1 1
    expect_equal(as.data.frame(bamflag(exdt))$isProperPair[1], 1)
    expect_equal(as.data.frame(bamflag(exdt))$isProperPair[2], 1)  
    ## isUnmappedQuery  0 0 
    expect_equal(as.data.frame(bamflag(exdt))$isUnmappedQuery[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isUnmappedQuery[2], 0) 
    ## hasUnmappedMate  0 0 
    expect_equal(as.data.frame(bamflag(gr_example))$hasUnmappedMate[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$hasUnmappedMate[2], 0) 
    ## isMinusStrand  0 1 
    expect_equal(as.data.frame(bamflag(gr_example))$isMinusStrand[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$isMinusStrand[2], 1) 
    ## isMateMinusStrand  1 0 
    expect_equal(as.data.frame(bamflag(gr_example))$isMateMinusStrand[1], 1)
    expect_equal(as.data.frame(bamflag(gr_example))$isMateMinusStrand[2], 0) 
    ## isFirstMateRead  1 0 
    expect_equal(as.data.frame(bamflag(gr_example))$isFirstMateRead[1], 1)
    expect_equal(as.data.frame(bamflag(gr_example))$isFirstMateRead[2], 0) 
    ## isSecondMateRead  0 1
    expect_equal(as.data.frame(bamflag(gr_example))$isSecondMateRead[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$isSecondMateRead[2], 1) 
    ## isNotPrimaryRead  0 0
    expect_equal(as.data.frame(bamflag(gr_example))$isNotPrimaryRead[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$isNotPrimaryRead[2], 0) 
    ## isNotPassingQualityControls  0 0 
    expect_equal(as.data.frame(bamflag(gr_example))$isNotPassingQualityControls[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$isNotPassingQualityControls[2], 0) 
    ## isDuplicate  0 0 
    expect_equal(as.data.frame(bamflag(gr_example))$isDuplicate[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$isDuplicate[2], 0) 
    ## isSupplementary  0 0
    expect_equal(as.data.frame(bamflag(gr_example))$isSupplementary[1], 0)
    expect_equal(as.data.frame(bamflag(gr_example))$isSupplementary[2], 0) 

})



## parse_outputs

test_that('parse_outputs', {

    ## generate an example output from query_bam_index()
    ## should output a DataFrame with columns "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual"
    foo = read.bam(example_bam, all=TRUE)
    gr2 = foo[[2]]
    dt2 = gr2dt(gr2)
    ## correct how scanBam() names BAM fields
    dt2$rnext = dt2$mrnm
    dt2$pnext = dt2$mpos
    dt2$tlen = dt2$isize
    dt2$rname = dt2$seqnames
    dt2$mrnm = NULL
    dt2$mpos = NULL
    dt2$isize = NULL  
    dt2$MD = NULL
    dt2$MQ = NULL  
    dt2$seqnames = NULL
    dt2$pos = dt2$start
    dt2$start = NULL
    dt2$end = NULL
    dt2$strand = NULL
    dt2$width = NULL
    dt2$qwidth = NULL
    output_list = as.list(dt2)  ## list with 11 mandatory BAM fields
    ## > names(output_list)
    ## [1] "qname"  "flag"   "qwidth" "mapq"   "cigar"  "seq"    "qual"   "rnext" 
    ## [9] "pnext"  "tlen"   "rname"  "pos"  
    expect_equal(width(parse_outputs(output_list))[1], 127)
    expect_equal(width(parse_outputs(output_list))[2], 150)
    expect_match(as.character(seqnames(parse_outputs(output_list))[1]), "chr5")
    expect_match(as.character(seqnames(parse_outputs(output_list))[2]), "chr5")

})



## check_index 
test_that('check_index', {

    expect_true(check_index(example_lmdb, "BX"))
    expect_true(check_index(example_lmdb, "QNAME"))
    expect_false(check_index(example_lmdb, "BZ"))

})


######## R6_class.R



















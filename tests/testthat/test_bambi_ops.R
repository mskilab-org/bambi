
library(testthat)
library(bambi)
library(bamUtils)
library(gUtils)


example_bam = 'subsetHCC1143_phased_possorted0001.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_bai = 'subsetHCC1143_phased_possorted0001.bam.bai' 
example_lmdb = 'subsetHCC1143_phased_possorted0001_lmdb'



#### utils.R


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

    ## isPaired  1 0
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isPaired[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isPaired[2], 0)
    ## isProperPair  1 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isProperPair[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isProperPair[2], 0)
    ## isUnmappedQuery  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isUnmappedQuery[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isUnmappedQuery[2], 0)
    ## hasUnmappedMate  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$hasUnmappedMate[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$hasUnmappedMate[2], 0)
    ## isMinusStrand  1 1 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMinusStrand[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMinusStrand[2], 0)
    ## isMateMinusStrand  1 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMateMinusStrand[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMateMinusStrand[2], 0)
    ## isFirstMateRead  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isFirstMateRead[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isFirstMateRead[2], 0)
    ## isSecondMateRead  1 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSecondMateRead[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSecondMateRead[2], 0)
    ## isNotPrimaryRead  0 0
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPrimaryRead[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPrimaryRead[2], 0)
    ## isNotPassingQualityControls  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPassingQualityControls[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPassingQualityControls[2], 0)
    ## isDuplicate  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isDuplicate[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isDuplicate[2], 0)
    ## isSupplementary  0 0
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSupplementary[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSupplementary[2], 0)

})



## parse_outputs

test_that('parse_outputs', {

    foo = read.bam(example_bam, all=TRUE)
    gr2 = foo[[2]]
    dt2 = gr2dt(gr2)
 

})








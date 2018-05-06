
library(testthat)
library(bambi)
library(bamUtils)
library(gUtils)

example_bam = 'subsetHCC1143_phased_possorted0001.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_bai = 'subsetHCC1143_phased_possorted0001.bam.bai' 
example_lmdb = 'subsetHCC1143_phased_possorted0001_lmdb'




######## R6_class.R



test_that('bambi class constructor', {
 
    ## Error in .subset2(public_bind_env, "initialize")(...) : 
    ##   BAM file not found. A valid BAM for 'bam_file' must be provided.
    expect_error(bambi$new(bam_file = 'foobar'))
    ## Error in .subset2(public_bind_env, "initialize")(...) : 
    ##   Cannot open BAM. A valid BAM for 'bam_file' must be provided.
    expect_error(bambi$new(bam_file = 'fake_bam.txt'))
    ## if (is.null(bamdb_path))
    ## Error in .subset2(public_bind_env, "initialize")(...) : 
    ##   Error: Pathname to bamdb LMDB subdirectory must be provided
    expect_error(bambi$new(bam_file ='bam_with_no_lmdb.bam'))
    ## valid
    foo = bambi$new(bam_file = example_bam)
    expect_match(foo$bam_file, 'subsetHCC1143_phased_possorted0001.bam')
    expect_match(foo$bamdb_path, 'subsetHCC1143_phased_possorted0001_lmdb')

})



test_that('bambi test method grab_bx()', {
 
    ## check works
    foo = bambi$new(bam_file = example_bam)
    expect_equal(foo$grab_bx('TACTCATCACACGCAC-1')$flag[1], 99)
    expect_equal(foo$grab_bx('TACTCATCACACGCAC-1')$flag[2], 147)
    expect_equal(foo$grab_bx('TACTCATCACACGCAC-1')$mapq[1], 60)
    expect_equal(foo$grab_bx('TACTCATCACACGCAC-1')$mapq[2], 60)
    expect_match(foo$grab_bx('TACTCATCACACGCAC-1')$cigar[1], "127M")
    expect_match(foo$grab_bx('TACTCATCACACGCAC-1')$cigar[2], "150M")
    ## if (check_index(self$bamdb_path, 'BX') != TRUE)
    ## wrong = bambi$new(bam_file = example_bam)
    ##
    ## if ((!is.null(barcodes)) & (!is.null(query))){
    expect_equal(foo$grab_bx(), GRanges())
    expect_equal(foo$grab_bx(data.table=TRUE), data.table())
    ## if (!inherits(barcodes, "character"))
    ## foo$grab_bx(barcodes = "2")
    ## one barcode


})





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














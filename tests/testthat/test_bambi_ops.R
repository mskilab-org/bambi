
library(testthat)
library(bambi)
library(bamUtils)
library(gUtils)


example_bam = 'subset_HCC1143_phased_possorted000005.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_bai = 'subset_HCC1143_phased_possorted000005.bam.bai' 
example_lmdb = 'subset_HCC1143_phased_possorted000005_lmdb'

example_CB_bam = 'subset_pbmc8k_possorted_genome_bam000025.CB.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_CB_bai = 'subset_pbmc8k_possorted_genome_bam000025.CB.bam.bai' 
example_CB_lmdb = 'subset_pbmc8k_possorted_genome_bam000025.CB_lmdb'

example_UB_bam = 'subset_pbmc8k_possorted_genome_bam000025.UB.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_UB_bai = 'subset_pbmc8k_possorted_genome_bam000025.UB.bam.bai' 
example_UB_lmdb = 'subset_pbmc8k_possorted_genome_bam000025.UB_lmdb'

example_PS_bam = 'subset_HCC1143_phased_possorted000005.PS.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_PS_bai = 'subset_HCC1143_phased_possorted000005.PS.bam.bai' 
example_PS_lmdb = 'subset_HCC1143_phased_possorted000005.PS_lmdb'

bam_no_lmdb = 'bam_with_no_lmdb.bam'

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
    expect_match(foo$bam_file, 'subset_HCC1143_phased_possorted000005.bam')
    expect_match(foo$bamdb_path, 'subset_HCC1143_phased_possorted000005_lmdb')

})



test_that('bambi test method grab_bx()', {
 
    ## check works
    foo = bambi$new(bam_file = example_bam)
    ## Error in foo$grab_cb() : 
    ## CB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foo$grab_cb())
    ## UB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foo$grab_ub())   
    ###
    expect_equal(foo$grab_bx('CGGCTAGGTAGAATAC-1', verbose=TRUE)$flag[1], 97)
    expect_equal(foo$grab_bx('CGGCTAGGTAGAATAC-1', verbose=TRUE)$flag[2], 145)
    expect_equal(foo$grab_bx('CGGCTAGGTAGAATAC-1', verbose=TRUE)$mapq[1], 60)
    expect_equal(foo$grab_bx('CGGCTAGGTAGAATAC-1', verbose=TRUE)$mapq[2], 60)
    expect_match(foo$grab_bx('CGGCTAGGTAGAATAC-1')$cigar[1], "80M47S")
    expect_match(foo$grab_bx('CGGCTAGGTAGAATAC-1')$cigar[2], "37S113M")
    ##
    ## if ((is.null(barcodes)) & (is.null(query))){
    expect_equal(foo$grab_bx(), GRanges())
    expect_equal(foo$grab_bx(data.table=TRUE), data.table())
    ##
    ## if ((!is.null(barcodes)) & (!is.null(query))){
    ## Error in foo$grab_bx(barcodes = "CGACGTGTCCTCTAGC-1", query = data.table()) : 
    ##   Both 'barcodes' and 'query' parameters cannot be used. Use method 'grab_bx()' by a character vector of BX barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.
    expect_error(foo$grab_bx(barcodes='CGACGTGTCCTCTAGC-1', query=data.table(), verbose=TRUE))
    expect_error(foo$grab_bx(barcodes='CGACGTGTCCTCTAGC-1', query=GRanges("chr5:1053000-1253655", verbose=TRUE)))
    ##
    ## empty query: R standard is NULL for error, NA for no rows returned
    expect_equal(foo$grab_bx(barcodes='2'), NA)
    expect_equal(foo$grab_bx(barcodes='foobar'), NA)
    ## multiple barcodes
    expect_equal(length(foo$grab_bx(barcodes=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1'), mc.cores=2, verbose=TRUE)), 8)
    expect_equal(dim(foo$grab_bx(barcodes=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1'), data.table=TRUE, mc.cores=2, verbose=TRUE))[1], 8)
    expect_equal(dim(foo$grab_bx(barcodes=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1'), data.table=TRUE, mc.cores=2, verbose=TRUE))[2], 22)
    ## 
    ## If there are "empty row" queries, I currently return nothing
    expect_equal(length(foo$grab_bx(barcodes=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1', 'foo', 'foobar', '2'), mc.cores=2)), 8)
    expect_equal(dim(foo$grab_bx(barcodes=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[1], 8)
    expect_equal(dim(foo$grab_bx(barcodes=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[2], 22)
    ##
    ## if (!inherits(barcodes, "character"))
    ## Error in foo$grab_bx(barcodes = 2) : 
    ##  Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.
    expect_error(foo$grab_bx(barcodes=2))
    ##
    ## incorrect query
    ## stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
    expect_error(foo$grab_bx(query=c('foo')))
    ## 
    ## if (length(query)==0){
    expect_equal(foo$grab_bx(query=GRanges('chr5:10-15'), verbose=TRUE), NA)
    expect_equal(foo$grab_bx(query=GRanges("chr5:1053000-1253655")), NA)
    ##
    ## else
    expect_equal(length(foo$grab_bx(query=GRanges('chr19:27842400-27842700'), verbose=TRUE)), 4)
    ##
})




## CB

## 'GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1'
## TCCCGATCACCAGTTA-1
## AGGTCCGAGGTACTCT-1

test_that('bambi test method grab_cb()', {
 
    ## check works
    foocb = bambi$new(bam_file = example_CB_bam)
    ## Error in foo$grab_cb() : 
    ## BX is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foocb$grab_bx())
    ## UB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foocb$grab_ub())   
    ### 
    ## somehow having an issue as it's a vector...
    ## expect_equal(foocb$grab_cb('GTAGTCATCTGGGCCA-1'))$flag, 16)  
    ## expect_equal(foocb$grab_cb('GTAGTCATCTGGGCCA-1'))$mapq, 255)
    ## expect_match(foocb$grab_cb('GTAGTCATCTGGGCCA-1'))$cigar, "98M")
    expect_equal(length(foocb$grab_cb('GTAGTCATCTGGGCCA-1')), 1)
    expect_equal(length(foocb$grab_cb('GTAGTCATCTGGGCCA-1', verbose=TRUE)), 1)
    ##
    ## if ((is.null(barcodes)) & (is.null(query))){
    expect_equal(foocb$grab_cb(), GRanges())
    expect_equal(foocb$grab_cb(data.table=TRUE), data.table())
    ##
    ## if ((!is.null(barcodes)) & (!is.null(query))){
    ## Error in foo$grab_bx(barcodes = "CGACGTGTCCTCTAGC-1", query = data.table()) : 
    ##   Both 'barcodes' and 'query' parameters cannot be used. Use method 'grab_bx()' by a character vector of BX barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.
    expect_error(foocb$grab_cb(barcodes='GTAGTCATCTGGGCCA-1', query=data.table(), verbose=TRUE))
    expect_error(foocb$grab_cb(barcodes='GTAGTCATCTGGGCCA-1', query=GRanges("chr5:1053000-1253655"), verbose=TRUE))
    ##
    ## empty query: R standard is NULL for error, NA for no rows returned
    expect_equal(foocb$grab_cb(barcodes='2'), NA)
    expect_equal(foocb$grab_cb(barcodes='foobar'), NA)
    ##
    ## multiple barcodes
    expect_equal(length(foocb$grab_cb(barcodes=c('GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1'), mc.cores=2, verbose=TRUE)), 10)
    expect_equal(dim(foocb$grab_cb(barcodes=c('GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1'), data.table=TRUE, mc.cores=2, verbose=TRUE))[1], 10)
    expect_equal(dim(foocb$grab_cb(barcodes=c('GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1'), data.table=TRUE, mc.cores=2))[2], 26)
    ## 
    ## If there are "empty row" queries, I currently return nothing
    expect_equal(length(foocb$grab_cb(barcodes=c('GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1', 'foo', 'foobar', '2'), mc.cores=2)), 10)
    expect_equal(dim(foocb$grab_cb(barcodes=c('GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[1], 10)
    expect_equal(dim(foocb$grab_cb(barcodes=c('GTAGTCATCTGGGCCA-1', 'TCCCGATCACCAGTTA-1', 'AGGTCCGAGGTACTCT-1', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[2], 26)
    ##     
    ## if (!inherits(barcodes, "character"))
    ##  Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.
    expect_error(foocb$grab_cb(barcodes=2))
    ##
    ## incorrect query
    ## stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
    expect_error(foocb$grab_cb(query=c('foo')))
    ## 
    ## if (length(query)==0){
    expect_equal(foocb$grab_cb(query=GRanges('5:10-15'), verbose=TRUE), NA)
    expect_equal(foocb$grab_cb(query=GRanges("5:1053000-1253655"), verbose=TRUE), NA)
    ##
    ## else
    expect_equal(length(foocb$grab_cb(query=GRanges("19:1440000-1440500"), verbose=TRUE)), 33)

})




## UB

## ATACAAGCGG
## CGGAGGACGT
## CATAGCGTTT
## 'ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT'

test_that('bambi test method grab_ub()', {
 
    ## check works
    fooub = bambi$new(bam_file = example_UB_bam)
    ## Error in foo$grab_cb() : 
    ## BX is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(fooub$grab_bx())
    ## CB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(fooub$grab_cb())   
    ###
    ## expect_equal(as.data.frame(fooub$grab_ub('ATACAAGCGG', verbose=TRUE))$flag, 0)
    ## expect_equal(as.data.frame(fooub$grab_ub('ATACAAGCGG'))$mapq, 255)
    ## expect_match(as.data.frame(fooub$grab_ub('ATACAAGCGG'))$cigar, "98M")
    expect_equal(length(fooub$grab_ub('ATACAAGCGG')), 1)
    expect_equal(length(fooub$grab_ub('ATACAAGCGG', verbose=TRUE)), 1)
    ##
    ## if ((is.null(barcodes)) & (is.null(query))){
    expect_equal(fooub$grab_ub(), GRanges())
    expect_equal(fooub$grab_ub(data.table=TRUE), data.table())
    ##
    ## empty query: R standard is NULL for error, NA for no rows returned
    expect_equal(fooub$grab_ub(barcodes='2'), NA)
    expect_equal(fooub$grab_ub(barcodes='foobar'), NA)
    ##
    ## if ((!is.null(barcodes)) & (!is.null(query))){
    ## Error in foo$grab_bx(barcodes = "CGACGTGTCCTCTAGC-1", query = data.table()) : 
    ##   Both 'barcodes' and 'query' parameters cannot be used. Use method 'grab_bx()' by a character vector of BX barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.
    expect_error(fooub$grab_ub(barcodes='ATACAAGCGG', query=data.table()))
    expect_error(fooub$grab_ub(barcodes='ATACAAGCGG', query=GRanges("chr5:1053000-1253655")))
    ##
    ## multiple barcodes
    ## expect_equal(length(fooub$grab_ub(barcodes=c('ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT'), mc.cores=2)), 3)
    expect_equal(dim(fooub$grab_ub(barcodes=c('ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT'), data.table=TRUE, mc.cores=2, verbose=TRUE))[1], 0)
    expect_equal(dim(fooub$grab_ub(barcodes=c('ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT'), data.table=TRUE, mc.cores=2, verbose=TRUE))[2], 12)
    ## 
    ## If there are "empty row" queries, I currently return nothing
    ## expect_equal(length(fooub$grab_ub(barcodes=c('ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT', 'foo', 'foobar', '2'), mc.cores=2)), 3)
    expect_equal(dim(fooub$grab_ub(barcodes=c('ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[1], 0)
    expect_equal(dim(fooub$grab_ub(barcodes=c('ATACAAGCGG', 'CGGAGGACGT', 'CATAGCGTTT', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[2], 12)
    ##     
    ## if (!inherits(barcodes, "character")) 
    ##  Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.
    expect_error(fooub$grab_ub(barcodes=2))
    ## if ((is.null(barcodes)) & (!is.null(query))){
    ##
    ## incorrect query
    ## stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
    expect_error(fooub$grab_ub(query=c('foo')))
    ##
    ## if (length(query)==0){
    expect_equal(fooub$grab_ub(query=GRanges('5:10-15'), verbose=TRUE), NA)
    expect_equal(fooub$grab_ub(query=GRanges("5:1053000-1253655")), NA)
    ##
    expect_equal(length(fooub$grab_ub(query=GRanges("19:1440000-1440500"), verbose=TRUE)), 2)
    ## else


})



## another tag, e.g. PS from HCC1143
## fetch_by_tag = function(tag, tag_queries=NULL, query=NULL, data.table=FALSE, verbose=FALSE, mc.cores=1)

test_that('bambi test method fetch_by_tag()', {
 
    ## check works
    foo_fetch_bx = bambi$new(bam_file = example_bam)
    ## Error in foo$grab_cb() : 
    ## CB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foo_fetch_bx$grab_cb())   
    ##
    ## Error in foo_fetch_bx$fetch_by_tag(tag = c("BX", "CB")) : 
    ## Invalid tag input. Multiple tags at once is not currently supported with 'fetch_by_tag()'. Please see documentation for details.
    expect_error(foo_fetch_bx$fetch_by_tag(tag = c('BX', 'CB')))
    ## 
    ## foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1')
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1', verbose=TRUE)$flag[1], 97)
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1', verbose=TRUE)$flag[2], 145)
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1')$mapq[1], 60)
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1')$mapq[2], 60)
    expect_match(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1')$cigar[1], "80M47S")
    expect_match(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1')$cigar[2], "37S113M")
    ##
    ## if (check_index(self$bamdb_path, tag) != TRUE){
    ## Error in foo_fetch_bx$fetch_by_tag(tag = c("CB")) : 
    ## The input 'tag' is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foo_fetch_bx$fetch_by_tag(tag = c('CB')))
    ##
    ## if ((!is.null(tag_queries)) & (!is.null(query))){
    ## Both 'tag_queries' and 'query' parameters cannot be used. Use method 'fetch_by_tag()' by either a character vector of UB barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details
    expect_error(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1', query=data.table(), verbose=TRUE))
    expect_error(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='CGGCTAGGTAGAATAC-1', query=GRanges("chr5:1053000-1253655"), verbose=TRUE))
    ## 
    ## if ((is.null(tag_queries)) & (is.null(query))){
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX'), GRanges())
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', data.table=TRUE), data.table())
    ##
    ## empty query: R standard is NULL for error, NA for no rows returned
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='2'), NA)
    expect_equal(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries='foobar'), NA)
    ##
    ## multiple barcodes
    expect_equal(length(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1'), mc.cores=2, verbose=TRUE)), 8)
    expect_equal(dim(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1'), data.table=TRUE, mc.cores=2, verbose=TRUE))[1], 8)
    expect_equal(dim(foo_fetch_bx$fetch_by_tag(tag = 'BX', tag_queries=c('CGGCTAGGTAGAATAC-1', 'AGCTTCCCACTTCACC-1', 'GGCGTTGAGAGCTGGT-1'), data.table=TRUE, mc.cores=2, verbose=TRUE))[2], 23)
    ## 
    ##
    foo_fetch_ps = bambi$new(bam_file = example_PS_bam)
    ## CB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.
    expect_error(foo_fetch_ps$grab_cb())  
    ##
    expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries='1013620572')$flag[1], 147)
    expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries='1013620572')$flag[2], 1171)
    expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries='1013620572')$mapq[1], 60)
    expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries='1013620572')$mapq[2], 60)
    expect_match(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries='1013620572')$cigar[1], "150M")
    expect_match(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries='1013620572')$cigar[2], "150M")
    ##
    ## multiple barcodes
    expect_equal(length(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries=c('10522462', '1013620572', '870145439217'), mc.cores=2)), 47)
    expect_equal(dim(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries=c('10522462', '1013620572', '870145439217'), data.table=TRUE, mc.cores=2))[1], 47)
    expect_equal(dim(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries=c('10522462', '1013620572', '870145439217'), data.table=TRUE, mc.cores=2))[2], 23)
    ## 
    ## If there are "empty row" queries, I currently return nothing
    expect_equal(length(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries=c('10522462', '1013620572', '870145439217', 'foo', 'foobar', '2'), mc.cores=2)), 47)
    expect_equal(dim(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries=c('10522462', '1013620572', '870145439217', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[1], 47)
    expect_equal(dim(foo_fetch_ps$fetch_by_tag(tag = 'PS', tag_queries=c('10522462', '1013620572', '870145439217', 'foo', 'foobar', '2'), data.table=TRUE, mc.cores=2))[2], 23)
    ##    
    ## if ((is.null(barcodes)) & (!is.null(query))){
    ##
    ## incorrect query
    ## stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
    expect_error(foo_fetch_ps$fetch_by_tag(tag = 'PS', query=c('foo')))
    ##
    ## if (length(query)==0){
    expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', query=GRanges('chr5:10-15'), verbose=TRUE), NA)
    expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', query=GRanges('chr5:1053000-1253655'), verbose=TRUE), NA)
    ##
    ## else
    ## FIX
    ## expect_equal(foo_fetch_ps$fetch_by_tag(tag = 'PS', query=GRanges('chr19:27842400-37842700')), NA)
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
    ## error? 'no method for coercing this S4 class to a vector'
    ##expect_equal(as.data.frame(bamflag(gr_example))$isPaired[1], 1)  ## no method for coercing this S4 class to a vector
    ##expect_equal(as.data.frame(bamflag(gr_example))$isPaired[2], 1)
    ## isProperPair  1 1
    expect_equal(as.data.frame(bamflag(exdt))$isProperPair[1], 1)
    expect_equal(as.data.frame(bamflag(exdt))$isProperPair[2], 1)  
    ## isUnmappedQuery  0 0 
    expect_equal(as.data.frame(bamflag(exdt))$isUnmappedQuery[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isUnmappedQuery[2], 0) 
    ## hasUnmappedMate  0 0 
    expect_equal(as.data.frame(bamflag(exdt))$hasUnmappedMate[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$hasUnmappedMate[2], 0) 
    ## isMinusStrand  0 1 
    expect_equal(as.data.frame(bamflag(exdt))$isMinusStrand[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isMinusStrand[2], 1) 
    ## isMateMinusStrand  1 0 
    expect_equal(as.data.frame(bamflag(exdt))$isMateMinusStrand[1], 1)
    expect_equal(as.data.frame(bamflag(exdt))$isMateMinusStrand[2], 0) 
    ## isFirstMateRead  1 0 
    expect_equal(as.data.frame(bamflag(exdt))$isFirstMateRead[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isFirstMateRead[2], 1) 
    ## isSecondMateRead  0 1
    expect_equal(as.data.frame(bamflag(exdt))$isSecondMateRead[1], 1)
    expect_equal(as.data.frame(bamflag(exdt))$isSecondMateRead[2], 0) 
    ## isNotPrimaryRead  0 0
    expect_equal(as.data.frame(bamflag(exdt))$isNotPrimaryRead[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isNotPrimaryRead[2], 0) 
    ## isNotPassingQualityControls  0 0 
    expect_equal(as.data.frame(bamflag(exdt))$isNotPassingQualityControls[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isNotPassingQualityControls[2], 0) 
    ## isDuplicate  0 0 
    expect_equal(as.data.frame(bamflag(exdt))$isDuplicate[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isDuplicate[2], 0) 
    ## isSupplementary  0 0
    expect_equal(as.data.frame(bamflag(exdt))$isSupplementary[1], 0)
    expect_equal(as.data.frame(bamflag(exdt))$isSupplementary[2], 0) 

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
    expect_equal(width(parse_outputs(output_list))[1], 82)
    expect_equal(width(parse_outputs(output_list))[2], 150)
    expect_match(as.character(seqnames(parse_outputs(output_list))[1]), "chr19")
    expect_match(as.character(seqnames(parse_outputs(output_list))[2]), "chr19")

})



## check_index 
test_that('check_index', {

    expect_true(check_index(example_lmdb, "BX"))
    expect_true(check_index(example_lmdb, "QNAME"))
    expect_false(check_index(example_lmdb, "BZ"))  ## only BX and QNAME should be indexed

})














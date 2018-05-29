#' @import data.table
#' @import GenomicRanges
#' @import IRanges
#' @import gUtils
#' @import bamUtils

#' @name countCigar
#' @rdname internal
#' @title Count bases in cigar string
#' @description
#'
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar() makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#'
#' @param cigar character vector of cigar strings
#' @return a 4-column, length(cigar)-row matrix with the total counts for each type
countCigar = function(cigar) {
    
    cigar = as.character(cigar)

    cigar.vals = unlist(strsplit(cigar, "\\d+"))
    cigar.lens = strsplit(cigar, "[A-Z]")
    lens = nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] = 1
    
    cigar.lens = as.numeric(unlist(cigar.lens))
    cigar.vals = cigar.vals[cigar.vals != ""]
    repr = rep(seq_along(cigar), lens)
    dt  = data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d = dt["D", ][, sum(lens), by=group]
    smr.i = dt["I", ][, sum(lens), by=group]
    smr.m = dt["M", ][, sum(lens), by=group]
    smr.s = dt["S", ][, sum(lens), by=group]
    
    out = matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group, 1] = smr.d$V1
    out[smr.i$group, 2] = smr.i$V1
    out[smr.m$group, 3] = smr.m$V1
    out[smr.s$group, 4] = smr.s$V1
    colnames(out) = c('D','I','M','S')
    
    return(out)
}


#' @name bamflag
#' @title Returns matrix of bits from BAM flags
#' @rdname internal
#' @description
#'
#' Shortcut function in create_granges
#'
#' @param reads data.table holding the reads
#' @return matrix of bits from BAM flags
#' @export
bamflag = function(reads){

    ## should only be type 'data.table'
    if (inherits(reads, 'data.table')){  
        bf = reads$flag
    } else{
        bf = reads
    }

    out = matrix(as.numeric(intToBits(bf)), byrow = TRUE, ncol = 32)[, 1:12, drop = FALSE]
    colnames(out) = c('isPaired', 'isProperPair', 'isUnmappedQuery', 'hasUnmappedMate', 'isMinusStrand', 'isMateMinusStrand', 'isFirstMateRead', 'isSecondMateRead', 'isNotPrimaryRead', 'isNotPassingQualityControls', 'isDuplicate', 'isSupplementary')

    return(out)
}




#' @name parse_outputs
#' @title parses output of BAM queries into GRanges with parsed CIGAR strings
#' @description
#'
#' Internal function to parse CIGAR and output GRanges
#'
#' @param out List output from query_bam_index()
#' @return GRanges with parsed CIGAR 
#' @importFrom IRanges IRanges
#' @import data.table
parse_outputs = function(out){

    out = as.data.table(out)    ### this should already be the case in the source; however, query_bam_index() outputs an R list
    cigs = countCigar(out$cigar)
    out[, pos2 := out$pos + rowSums(cigs[, c('D', 'M'), drop = FALSE], na.rm = TRUE) - 1]
    
    out$qwidth = nchar(as.character(out$seq))
    out$strand = bamflag(out$flag)[, 'isMinusStrand'] == 1
    out$strand = ifelse(out$strand, '-', '+')
    
    unmapped = bamflag(out$flag)[, 'isUnmappedQuery'] == 1

    if (any(unmapped)){
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = '*'
    }
    
    bf = out$flag
    
    ## create data.table of start, end, strand, seqnames
    newdt = data.table(pos = out$pos, pos2 = out$pos2, strand = out$strand, rname = out$rname)   
    
    rr = IRanges(newdt$pos, newdt$pos2)
    sf = factor(newdt$strand, levels=c('+', '-', '*'))
    ff = factor(newdt$rname, levels=unique(newdt$rname))
    gr.fields = c('rname', 'strand', 'pos', 'pos2')
    grobj = GRanges(seqnames=ff, ranges=rr, strand=sf)
    
    vals = out[, setdiff(names(out), gr.fields), with=FALSE]   
    values(grobj) = vals
    return(grobj)
}



#' @name check_index
#' @title checks the LMDB indices to check whether these exist
#' @description
#'
#' Checks the LMDB indices to check whether these exist
#'
#' @param lmdb_path List output from query_bam_index()
#' @param index character string LMDB index to check if exists
#' @return boolean TRUE if LMDB has index 
check_index = function(lmdb_path, index){
    ## the subdirectories in the LMDB path should be the index names
    indices = gsub(".*/", "", list.dirs(path=lmdb_path, recursive=FALSE))  ## re-format list.dirs() outputs
    return(as.logical(is.element(index, indices)))
}







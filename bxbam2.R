#############################################################################
## Evan Biederstedt
## New York Genome Center
## ebiederstedt@nygenome.org
##
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
## New York Genome Center
## mimielinski@nygenome.org
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

#' @import rPython
#' @import data.table
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicAlignments



#' @name countCigar
#' @title countCigar
#' @description
#'
#' Count bases in cigar string
#'
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#' @param cigar character vector of cigar strings
#' @return a 4-column, length(cigar)-row matrix with the total counts for each type
#' @export

countCigar <- function(cigar) {
    
    cigar.vals <- unlist(strsplit(cigar, "\\d+"))
    cigar.lens <- strsplit(cigar, "[A-Z]")
    lens <- nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] <- 1
    
    cigar.lens <- as.numeric(unlist(cigar.lens))
    cigar.vals <- cigar.vals[cigar.vals != ""]
    repr       <- rep(seq_along(cigar), lens)
    dt         <- data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d      <- dt["D",][, sum(lens), by=group]
    smr.i      <- dt["I",][, sum(lens), by=group]
    smr.m      <- dt["M",][, sum(lens), by=group]
    smr.s      <- dt["S",][, sum(lens), by=group]
    
    out <- matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group,1] <- smr.d$V1
    out[smr.i$group,2] <- smr.i$V1
    out[smr.m$group,3] <- smr.m$V1
    out[smr.s$group,4] <- smr.s$V1
    colnames(out) <- c('D','I','M','S')
    
    return(out)
}


#' @name bxbam-class
#' @title bxbam-class
#' @rdname bxbam-class
#' @description
#'
#' class::bxbam
#'
#' Class \code{bxbam} object to load bxbam and query GRanges
#'
#' @import rPython
#' @import data.table
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicAlignments
#' @exportClass bxbam
#' @author Evan Biederstedt

# alls the checks go here????????????????????????



setClass('bxBam', representation(.bxbamfile = 'character', sessionId = 'integer'))

setMethod('initialize', function(.Object,
                                   bxbamfile, ## Input HDF5 file to be read, created by script 'create_bxbam.py')
{
    ## warning: still very low risk of concurrency issue / collision if two parallel processes instantitating objects
    
    .Object@.sessionId <- paste0('session', runif(1))
    python.exec(sprintf("sessions['%s'] = tables.open_file('%s').get_node('/bam_table/bam_fields')", .Object@.sessionId, bxbamfile))
    
    ## onLoad
    #python.exec("import pandas as pd")
    #python.exec("import numpy as np")
    #python.exec("import tables")
    #python.exec("sessions = Set()")
    #python.exec("queries = Set()")
    #python.exec("num_sessions = 0")
    #python.exec("num_queries = 0")
    
    #python.exec("from tables import *")
    #python.assign("bxbamfile", bxbamfile)
    #python.exec("h5f = tables.open_file(bxbamfile)")
    #python.assign("bxbamfile", bxbamfile)
    #python.exec("h5f = tables.open_file(bxbamfile)")
    #python.exec("tb1 = h5f.get_node('/bam_table/bam_fields')")
    # .Object@sessionId <- python.get("num_sessions")
    #python.exec("num_session = num_sessions + 1")
    #python.exec("sessions.add(h5f.get_node('/bam_table/bam_fields'), num_session)")
    #.Object@sessionId <- python.get("num_sessions")
}
# constructor: new(Class="bxbam",bxbamfile = PATHNAME)

setValidity("bxBam", function(object){
    if (!file.exists(bxbamfile)) stop ("'bxbamfile' is missing; please set correct path to bxbam.h5 object")
}


#' @name get_bmates
#' @title get_bmates
#' @description
#' Grab all barcoded reads by bam field BX, return GRanges object
#'
#' @export
#' @author Evan Biederstedt

setMethod("get_bmates", "bxBam", function(.Object, query){


    if (inherits(query, 'GRanges') | inherits(query, 'data.frame'))
    {
        if (is.null(query$BX))
            stop('query must be either BX or GRanges / data.frame with BX field')
        query = query$BX
    }
   
   qstring = paste(paste0('(BX==b"', query, '")'), collapse = "|")
   queryId = paste0('query', runif(1))
   python.exec(sprintf("queries['%s'] = run_query(sessions['%s'], '%s')", queryId, .Object@sessionId, qstring))
   #.Object@.sessionId <- python.get("session_id")
   
    #queryId <- python.get("num_queries")
    #python.exec("num_queries = num_queries+1")
   
   
   query = "Z:GATAACCAGCGCTCAC-1"   # provided by user
   
   #python.exec(sprintf("queries['%s']= pd.DataFrame(tb1.read_where('%s'))", queryId, paste(paste0('(BX==b"', query, '")'), collapse = "|")))
   
   out = data.table(
        bx = python.get(sprintf("queries['%s'].BX.tolist()", queryId)),
        cigar = python.get(sprintf("queries['%s'].CIGAR.tolist()", queryId)),
        flag = python.get(sprintf("queries['%s'].FLAG.tolist()", queryId)),
        mapq = python.get(sprintf("queries['%s'].MAPQ.tolist()", queryId)),
        pnext = python.get(sprintf("queries['%s'].PNEXT.tolist()", queryId)),
        pos = python.get(sprintf("queries['%s'].POS.tolist()", queryId)),
        qname = python.get(sprintf("queries['%s'].QNAME.tolist()", queryId)),
        qual = python.get(sprintf("queries['%s'].QUAL.tolist()", queryId)),
        rname = python.get(sprintf("queries['%s'].RNAME.tolist()", queryId)),
        rnext = python.get(sprintf("queries['%s'].RNEXT.tolist()", queryId)),
        seq = python.get(sprintf("queries['%s'].SEQ.tolist()", queryId)),
        tlen = python.get(sprintf("queries['%s'].TLEN.tolist()", queryId)))


   #    python.exec('BX1 = tb1.read_where(\'BX==b"Z:GATAACCAGCGCTCAC-1"\')')
   # "queries[1]= tb1.read_where('(BX==b\"Z:GATAACCAGCGCTCAC-1\")|(BX==b\"Z:GATAACCAGCGCTCAC-1\")|(BX==b\"Z:GATAACCAGCGCTCAC-1\")')"
   #    python.exec("df1 = pd.DataFrame(BX1)")
   # one approach --- hard code each column---probably quickest---# can only run this once per dataframe!!!
   
//    
//    python.exec('df1.BX = df1.BX.str.decode("utf-8")')
//    python.exec('df1.CIGAR = df1.CIGAR.str.decode("utf-8")')
//    python.exec('df1.QNAME = df1.QNAME.str.decode("utf-8")')
//    python.exec('df1.QUAL = df1.QUAL.str.decode("utf-8")')
//    python.exec('df1.RNAME = df1.RNAME.str.decode("utf-8")')
//    python.exec('df1.RNEXT = df1.RNEXT.str.decode("utf-8")')
//    python.exec('df1.SEQ = df1.SEG.str.decode("utf-8")')      # error here!~!~!!!!!!
//    
//    
//    # convert to Python lists...
//    

    # ...then access as R vectors
//    
//    BX_vec <- python.get("BX_list")
//    CIGAR_vec <- python.get("CIGAR_list")
//    FLAG_vec <- python.get("FLAG_list")
//    MAPQ_vec <- python.get("MAPQ_list")
//    PNEXT_vec <- python.get("PNEXT_list")
//    POS_vec <- python.get("POS_list")
//    QNAME_vec <- python.get("QNAME_list")
//    QUAL_vec <- python.get("QUAL_list")
//    RNAME_vec <- python.get("RNAME_list")
//    RNEXT_vec <- python.get("RNEXT_list")
//    SEQ_vec <- python.get("SEQ_list")
//    TLEN_vec <- python.get("TLEN_list")
//    

//    qname <- python.get("QNAME_list")
//    flag  <- python.get("FLAG_list")
//    rname <- python.get("RNAME_list")
//    pos <- python.get("POS_list")
//    mapq <- python.get("MAPQ_list")
//    cigar <- python.get("CIGAR_list")
//    rnext <- python.get("RNEXT_list")
//    pnext <- python.get("PNEXT_list")
//    tlen <- python.get("TLEN_list")
//    seq <- python.get("SEQ_list")
//    qual <- python.get("QUAL_list")
//    bx  <- python.get("BX_list")
//    
    # create data.table
    
    #out <- data.table(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, bx)
    
    if (any(nnix <<- out$cigar=='*'))
        out$cigar[nnix] = NA
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1  # gives up pos2
    # M" positions in the alignment are encoded in the CIGAR; if add to the BAM pos you get the end position of the interval
    
    
    out$qwidth = nchar(out$seq)
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")

    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = '*'   # doesn't out why?
    }
    
    bf = out$flag
    
    newdt <- data.table(pos = out$pos, pos2 = out$pos2, strand = out$strand, rname = out$rname)   # create data.table of start, end, strand, seqnames
    
    rr <- IRanges(newdt$pos, newdt$pos2)
    sf <- factor(newdt$strand, levels=c('+', '-', '*'))
    ff <- factor(newdt$rname, levels=unique(newdt$rname))
    gr.fields <- c("rname", "strand", "pos", "pos2")
    grobj <- GRanges(seqnames=ff, ranges=rr, strand=sf)
    
    vals = out[, setdiff(names(out), gr.fields), with=FALSE]
    values(grobj) <- vals
    return(grobj)
}


#' @name get_qmates
#' @title get_qmates
#' @description
#' Grab all QNAMEs, return GRanges object
#'
#' @export
#' @author Evan Biederstedt

setMethod("get_qmates",'bxbam', function(reads){                 #
    python.exec('QNAME1 = tb1.read_where(\'QNAME1==b"Z:GATAACCAGCGCTCAC-1"\')')
    
    python.exec("df1 = pd.DataFrame(QNAME1)")
    
    # one approach --- hard code each column---probably quickest---# can only run this once per dataframe!!!
    
    python.exec('df1.BX = df1.BX.str.decode("utf-8")')
    python.exec('df1.CIGAR = df1.CIGAR.str.decode("utf-8")')
    python.exec('df1.QNAME = df1.QNAME.str.decode("utf-8")')
    python.exec('df1.QUAL = df1.QUAL.str.decode("utf-8")')
    python.exec('df1.RNAME = df1.RNAME.str.decode("utf-8")')
    python.exec('df1.RNEXT = df1.RNEXT.str.decode("utf-8")')
    python.exec('df1.SEQ = df1.SEG.str.decode("utf-8")')      # error here!~!~!!!!!!
    
    
    # convert to Python lists...
    
    python.exec(" BX_list = df1.BX.tolist()")
    python.exec(" CIGAR_list = df1.CIGAR.tolist()")
    python.exec(" FLAG_list = df1.FLAG.tolist()")
    python.exec(" MAPQ_list = df1.MAPQ.tolist()")
    python.exec(" PNEXT_list = df1.PNEXT.tolist()")
    python.exec(" POS_list = df1.POS.tolist()")
    python.exec(" QNAME_list = df1.QNAME.tolist()")
    python.exec(" QUAL_list = df1.QUAL.tolist()")
    python.exec(" RNAME_list = df1.RNAME.tolist()")
    python.exec(" RNEXT_list = df1.RNEXT.tolist()")
    python.exec(" SEQ_list = df1.SEQ.tolist()")
    python.exec(" TLEN_list = df1.TLEN.tolist()")
    
    
    # ...then access as R vectors
    
    BX_vec <- python.get("BX_list")
    CIGAR_vec <- python.get("CIGAR_list")
    FLAG_vec <- python.get("FLAG_list")
    MAPQ_vec <- python.get("MAPQ_list")
    PNEXT_vec <- python.get("PNEXT_list")
    POS_vec <- python.get("POS_list")
    QNAME_vec <- python.get("QNAME_list")
    QUAL_vec <- python.get("QUAL_list")table
    RNAME_vec <- python.get("RNAME_list")
    RNEXT_vec <- python.get("RNEXT_list")
    SEQ_vec <- python.get("SEQ_list")
    TLEN_vec <- python.get("TLEN_list")
    
    
    qname <- python.get("QNAME_list")
    flag  <- python.get("FLAG_list")
    rname <- python.get("RNAME_list")
    pos <- python.get("POS_list")
    mapq <- python.get("MAPQ_list")
    cigar <- python.get("CIGAR_list")
    rnext <- python.get("RNEXT_list")
    pnext <- python.get("PNEXT_list")
    tlen <- python.get("TLEN_list")
    seq <- python.get("SEQ_list")
    qual <- python.get("QUAL_list")
    bx  <- python.get("BX_list")
    
    
    # create data.table
    
    out <- data.table(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, bx)
    
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1  # gives up pos2
    # M" positions in the alignment are encoded in the CIGAR; if add to the BAM pos you get the end position of the interval
    
    
    out$qwidth = nchar(out$seq)
    
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")
    
    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = '*'   # doesn't out why?
    }
    
    bf = out$flag
    
    newdt <- data.table(pos = out$pos, pos2 = out$pos2, strand = out$strand, rname = out$rname)   # create data.table of start, end, strand, seqnames
    
    rr <- IRanges(newdt$pos, newdt$pos2)
    sf <- factor(newdt$strand, levels=c('+', '-', '*'))
    ff <- factor(newdt$rname, levels=unique(newdt$rname))
    gr.fields <- c("rname", "strand", "pos", "pos2")
    grobj <- GRanges(seqnames=ff, ranges=rr, strand=sf)
    
    vals = out[, setdiff(names(out), gr.fields), with=FALSE]
    values(grobj) <- vals
    return(grobj)
}





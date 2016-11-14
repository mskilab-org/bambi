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
#' @exportClass bxBam
#' @author Evan Biederstedt
setClass("bxBam", representation(.bxbamfile = 'character', .bamfile = 'BamFile', .sessionId = 'character'))

setMethod('initialize', 'bxBam', function(.Object, bxbamfile, bamfile = NULL)
{
    .Object@.sessionId <- paste0('session', runif(1))
    .Object@.bxbamfile = normalizePath(bxbamfile)
    message(.Object@.bxbamfile)

    if (is.null(bamfile)) # will try and guess the bamfile name from the bxbam
    {
        bamfile = gsub('bxbamfile', 'bxbam', bxbamfile)
        if (!file.exists(bamfile))
            bamfile = paste0(bamfile, '.bam')
        
        .Object@.bamfile = BamFile(bamfile)
    }    
    else ## if bam file is explicitly provided
    {
        if (is(bamfile, 'BamFile'))
            .Object@.bamfile = bamfile
        else
            .Object@.bamfile = BamFile(bamfile)

    }
    python.exec(sprintf("sessions['%s'] = tables.open_file('%s').get_node('/bam_table/bam_fields')", .Object@.sessionId, .Object@.bxbamfile))

    return(.Object)
})

#' @name bxBam
#' @title bxBam
#' @description
#' Initialize bxBam object specifying .bxbam file and optional .bam path
#' @export
#' @author Marcin Imielinski
bxBam = function(bxbamfile, bamfile = NULL)
    new('bxBam', bxbamfile, bamfile)


setValidity("bxBam", function(object){
    if (!file.exists(bxbamfile)) stop ("'bxbamfile' is missing; please set correct path to bxbam.h5 object")
})

#' @name show
#' @title show
#' @description Display a \code{gTrack} object
#' @docType methods
#' @param object \code{gTrack} to display
#' @author Marcin Imielinski
setMethod('show', 'bxBam', function(object)
{
  cat(sprintf('bxBam object stored in HDF5 file\n\t%s\nwith associated BAM file\n\t%s:\n', object@.bxbamfile, Rsamtools::path(object@.bamfile)))
})

#' @name get_bmates
#' @title get_bmates
#' @description
#' Grab all barcoded reads by bam field BX, return GRanges object
#'
#' @export
#' @author Evan Biederstedt
#' @author Marcin Imielinski
setGeneric('get_bmates', function(.Object, query) standardGeneric('get_bmates'))
setMethod("get_bmates", "bxBam", function(.Object, query){
    
    if (inherits(query, 'GRanges') | inherits(query, 'data.frame'))
    {
        if (is.null(query$BX))
            {
                warning("BX field not found, will use read.bam to pull reads under query GRanges from bam file and find their bmates")
                query = read.bam(.Object@.bamfile, gr = query, tag = c('BX', 'MD'), pairs.grl = FALSE)              
            }
        
        query = query$BX            
    }

    if (any(nix <<- is.na(query)))
        query = query[!nix]

    if (length(query)==0)
        stop('Length 0 query, check input')
            
   qstring = paste(paste0('(BX==b"', query, '")'), collapse = "|")
   queryId = paste0('query', runif(1))
   python.exec(sprintf("queries['%s'] = run_query(sessions['%s'], '%s')", queryId, .Object@.sessionId, qstring))
   
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
    
    if (any(nnix <<- out$cigar=='*'))
        out$cigar[nnix] = NA

    if (nrow(out)==0)
        return(GRanges())
    
    browser()
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1
    
    
    out$qwidth = nchar(out$seq)
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")

    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = "*"
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
    return(grobj)
})



#' @name get_qmates
#' @title get_qmates
#' @description
#' Grab all barcoded reads by bam field QNAME, return GRanges object
#'
#' @export
#' @author Evan Biederstedt
#' @author Marcin Imielinski
setGeneric('get_qmates', function(.Object, query) standardGeneric('get_qmates'))
setMethod("get_qmates", "bxBam", function(.Object, query){
    
    if (inherits(query, 'GRanges') | inherits(query, 'data.frame'))
    {
        if (is.null(query$BX))
        {
            warning("BX field not found, will use read.bam to pull reads under query GRanges from bam file and find their bmates")
            query = read.bam(.Object@.bamfile, gr = query, tag = c('BX', 'MD'), pairs.grl = FALSE)    
        }      
        query = query$QNMAE
    }
    
    if (any(nix <<- is.na(query)))
        query = query[!nix]

    if (length(query)==0)
        stop('Length 0 query, check input')

   
    qstring = paste(paste0('(QNAME==b"', query, '")'), collapse = "|")
    queryId = paste0('query', runif(1))
    python.exec(sprintf("queries['%s'] = run_query(sessions['%s'], '%s')", queryId, .Object@.sessionId, qstring))
  
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
    
    
    if (any(nnix <<- out$cigar=='*'))
    out$cigar[nnix] = NA
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1
    
    out$qwidth = nchar(out$seq)
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")
    
    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = "*"
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
    return(grobj)
})




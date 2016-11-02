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
#' @name load_bxbam
#' @title load_bxbam
#' @description
#'
#' Reads in bxbam object
#'
#' Load the bxbam object (HDF5 format) to be queried via pandas HDFStore() object (PyTables)
#' The default key ('HDF5 group') created using the Python script `create_bxbam.py` is hdf_key = "bam_fields"
#' @param bxbamfile Input HDF5 file to be read, created by script 'create_bxbam.py'
#' @param bam Input bam file
#' @param bami Input bam index file
#' @param hdf_key key defined upon creation of bxBam via 'create_bxbam.py'; default set to hdf_key = "bam_fields"
#' @return Noned
#' @examples
#' load_bxbam("pathname/my_bxbam.h5")
#' @export

load_bxbam <- function(bxbamfile, bam, bami = NULL){     # integreates functionality from BamUtils read.bam
    if (missing(bxbamfile)) stop ("'bxbamfile' is missing; please set correct path to bxbam.h5 object")\
    if (missing(bam)) stop ("'bam' is missing; please set correct path to bam file")
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
    python.assign("bxbamfile", bxbamfile)
    python.exec("store = pd.HDFStore(bxbamfile)")
    python.get("print('Now loaded bxBam file')")
    # python.get("store")   # 'store' variable now saved within Python  ; commented out, Nov. 1 2016
    if (!inherits(bam, 'BamFile'))  # check for bam, bami
    {
        if (is.null(bami))
        {
            if (file.exists(bai <- gsub('.bam$', '.bai', bam)))
            bam = BamFile(bam, bai)
            else if (file.exists(bai <- paste(bam, '.bai', sep = '')))
            bam = BamFile(bam, bai)
            else
            bam = BamFile(bam)
        }
        else
        bam = BamFile(bam, index = bami)
    }
}
#' @name get_bmates
#' @title get_bmates
#' @description
#'
#' Reads R vector of BX strings, outputs GRanges or data.table
#'
#' Function used to output all corresponding bam fields and GRanges via data.table
#' @param reads R vector of BX strings, e.g. bmates_vec <- c("BX1", "BX2", "BX3", "BX4")
#' @param hdf_key key defined upon creation of bxBam via 'create_bxbam.py'; use same hdf_key as bxbam in load_bxbam()
#' @param show_bam    Default FALSE
#' @param as.data.table Return reads in the form of a data.table rather than GRanges/GRangesList
#' @param ignore.indels messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
#' @param as.grl Return reads as GRangesList. Controls whether \code{get.pairs.grl} does split. Default TRUE
#' @return Corresponding fields in a bam file in data.table, GRanges, or GRangesList
#' @export
get_bmates <- function(reads, bxbam_key = "bam_fields", pairs.grl = TRUE, ignore.indels= FALSE, as.data.table = FALSE){   # show_bam = FALSE,
    if (is.null(reads))
    {
        stop('Must provide non empty read vector')
    }
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(bxbam_key, where='BX in @reads')")
    df <- python.get("df")    #!!!!!!!!!!!!!! is this crappy R code/style?
    df <- data.table(df)   # outputs data.table of bam rows   #!!!!!!!!!!!!!! is this crappy R code/style?
    
    # we now have data.frame of BAM fields (in the format as out in bamUtils
    #
    
    ## faster CIGAR string parsing with vectorization and data tables
    if (verbose) {
        print(Sys.time() - now)
        print('filling pos2 from cigar')
    }
    if (ignore.indels) {
        cigar <- gsub('[0-9]+D', '', gsub('([0-9]+)I', '\\1M', out$cigar))  ## Remove deletions, turn insertions to matches
        cig <- splitCigar(cigar)
        torun=sapply(cig, function(y) any(duplicated((y[[1]][y[[1]]==M]))))
        M <- charToRaw('M')
        new.cigar <- sapply(cig[torun], function(y) {
            lets <- y[[1]][!duplicated(y[[1]])]
            vals <- y[[2]][!duplicated(y[[1]])]
            vals[lets==M] <- sum(y[[2]][y[[1]]==M])
            lets <- strsplit(rawToChar(lets), '')[[1]]
            paste(as.vector(t(matrix(c(vals, lets), nrow=length(vals), ncol=length(lets)))), collapse='')
        })
        out$cigar[torun] <- new.cigar
    }
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + cigs[, "M"]
    
    if (verbose) {
        print(Sys.time() - now)
        print('fixing seqdata')
    }
    out$qwidth = nchar(out$seq)
    unm = is.na(out$pos)
    if (any(unm))
    {
        out$pos[unm] = 1
        out$pos2[unm] = 0
        out$strand[unm] = '*'
    }
    gr.fields = c('rname', 'strand', 'pos', 'pos2');
    vals = out[, setdiff(names(out), gr.fields)]
    
    if (!as.data.table) {
        out <- GRanges(out$rname, IRanges(out$pos, pmax(0, out$pos2-1)), strand = out$strand, seqlengths = seqlengths(intervals))
        values(out) <- vals;
    } else {
        out <- data.table(seqnames=out$rname, start=out$pos, end= pmax(out$pos2-1, 0), strand=out$strand)
        val <- data.table(vals)
        out <- cbind(out, val)
    }
    #out$uname = paste(out$qname, ifelse(bamflag(out$flag)[, 'isFirstMateRead'], '_r1', '_r2'), sep = '')
    }
    else {
    if (!as.data.table)
    return(GRanges(seqlengths = seqlengths(intervals)))
    else
    return(data.table())
}

if (verbose)
{
    if (as.data.table)
    cat(sprintf('Extracted %s reads\n', nrow(out)))
    else
    cat(sprintf('Extracted %s reads\n', length(out)))
    print(Sys.time() - now)
}

if (pairs.grl)
{
    if (verbose)
    cat('Pairing reads\n')
    out <- get.pairs.grl(out, as.grl=as.grl)
    if (verbose)
    {
        cat('done\n')
        print(Sys.time() - now)
    }
    if (as.grl && !as.data.table) {
        names(out) = NULL;
        values(out)$col = 'gray';
        values(out)$border = 'gray';
    }
}

return(out)
}


}

## Ditto with above----FIX!!!!!

#' @name get_qmates
#' @title get_qmates
#' @description
#'
#' Reads R vector of QNAME strings, outputs GRanges or data.table
#'
#' Function used to output all corresponding bam fields and GRanges via data.table
#' @param reads R vector of QNAME strings, e.g. bmates_vec <- c("QNAME1", "QNAME2", "QNAME3", "QNAME4")
#' @param hdf_key key defined upon creation of bxBam via 'create_bxbam.py'; use same hdf_key as bxbam in load_bxbam()
#' @param show_bam    Default False
#' @return Corresponding fields in a bam file in data.table, GRanges, or GRangesList
#' @export
get_qmates <- function(reads, bxbam_key = "bam_fields", show_bam = FALSE){
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(bxbam_key, where='QNAME in @reads'")
    df <- python.get("df")
    df <- data.table(df)
    return(df)
}


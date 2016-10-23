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
#'
#' An S4 class to represent BxBam
#'
#
#' @slot bxbam_file pathname to name of bxbam file created with create_bxbam.py
#' @slot bam        pathname to associated bam file
#' @slot bai        pathname to name of bxbam file created with create_bxbam.py
#' @slot bxkey      key of bxbam file; default set to  "bxbam_default_key"
setClass(
    Class = "BxBam",
    representation = representation(
        bxbam_file = "character",
        bam = "character",
        bami = "character",
        bxkey = "character"),
    validity= function(object) {
        if (!inherits(object@bam, "BamFile"))  # check for bam, bami
        {
            if (is.null(object@bami))
            {
                if (file.exists(bai <- gsub(".bam$", ".bai", object@bam)))
                bam = BamFile(object@bam, bai)
                else if (file.exists(bai <- paste(object@bam, ".bai", sep = "")))
                bam = BamFile(object@bam, bai)
                else
                bam = BamFile(object@bam)
            }
        else
        bam = BamFile(object@bam, index = bami)
    }
}
# define the constructor to "instantiate/load" the bxbam file
setGeneric(name="load_bxbam",
    def = function(object, bxbam_file, bam, bami = NULL, bxbam_key = "bxbam_default_key")   # here's where we set default?
    {
    standardGeneric("load_bxbam")
    }
    )
# create method for constructor
setMethod(f="load_bxbam",
    signature = "BxBam",
    definition = function(object, bxbam_file, bam, bami = NULL, bxbam_key = "bxbam_default_key")
    {
        python.exec("import pandas as pd")
        python.exec("import numpy as np")
        python.exec("import tables")
        python.exec("from tables import *")
        python.assign("bxbam_file", object@bxbam_file)
        python.exec("store = pd.HDFStore(bxbam_file)")
        python.get("print('Loaded bxBam file')")
        python.get("store")   # 'store' variable now saved within Python
        object@bxbamfile <- python.get("store")
        return(object)
    }
    )
# create method for get_bmates
setGeneric(name="get_bmates",
    def = function(object)
    {
        standardGeneric("get_bmates")
    }
    )
setMethod(f="get_bmates",
    definition = function(object, )
    {
        
    }
    )


# create method for get_qmates
setGeneric(name="get_qmates",
    def = function(object)
    {
        standardGeneric("get_qmates")
    }
    )
setMethod(f="get_bmates",
    definition = function(object, )
    {
    }
    )


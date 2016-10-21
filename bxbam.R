#############################################################################
## Evan Biederstedt
## New York Genome Center
## ebiederstedt@nygenome.org
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
###############################################################################

#' @import rPython
#' @import data.table
#' @name load_bxbam
#' @title load_bxbam
#' @description
#'
#' Reads in bxbam object
#'
#' Load the bxbam object (HDF5 format) to be queried via pandas HDFStore() object (PyTables)
#' The default key ('HDF5 group') created using the Python script `create_bxbam.py` is hdf_key = "bam_fields"
#' @param bxbam_file Input HDF5 file to be read, created by script 'create_bxbam.py'
#' @param hdf_key key defined upon creation of bxBam via 'create_bxbam.py'; default set to hdf_key = "bam_fields"
#' @return None
#' @examples
#' load_bxbam("pathname/my_bxbam.h5")
#' @export
load_bxbam <- function(bxbam_file, bxbam_key = "bam_fields"){     # e.g. bxbam_file is "file1.h5"
    if(missing(bxbam_key)){
        bx
    }
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
    python.assign("bxbam_file", bxbam_file)
    python.exec("store = pd.HDFStore(bxbam_file)")
    python.get("print('Loaded bxBam file')")
    python.get("store")   # 'store' variable now saved within Python
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
#' @param show_bam    Default TRUE
#' @return Corresponding fields in a bam file in data.table, GRanges, or GRangesList
#' @export
get_bmates <- function(reads, bxbam_key = "bam_fields", show_bam = FALSE){
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(bxbam_key, where='BX in @reads'")
    df <- python.get("df")
    df <- data.table(df)
    return(df)
}

#' @name get_qmates
#' @title get_qmates
#' @description
#'
#' Reads R vector of QNAME strings, outputs GRanges or data.table
#'
#' Function used to output all corresponding bam fields and GRanges via data.table
#' @param reads R vector of QNAME strings, e.g. bmates_vec <- c("QNAME1", "QNAME2", "QNAME3", "QNAME4")
#' @param hdf_key key defined upon creation of bxBam via 'create_bxbam.py'; use same hdf_key as bxbam in load_bxbam()
#' @param show_bam    Default TRUE
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


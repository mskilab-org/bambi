#' @import rPython
#' @import data.table
#' @import python.exec("import pandas as pd")
#' @import python.exec("import numpy as np")
#' @import python.exec("import tables")

#' @export hdf_key
#' key defined upon creation of bxBam via 'create_bxbam.py'
#' If users changed this value in 'create_bxbam.py', override here
hdf_key = "bam_fields"
python.assign("hdf_key", hdf_key)

#' @name load_bxbam
#' @title load_bxbam
#' @description
#'
#' Load the bxbam object (HDF5 format) to be queried via pandas HDFStore() object (PyTables)
#' The default key ('HDF5 group') created using the Python script `create_bxbam.py` is hdf_key = "bam_fields"
#' @param bxbam_file Input HDF5 file to be read, created by script 'create_bxbam.py'
#' @return Reads in bxbam objects
#' @export
load_bxbam <- function(bxbam_file){     # e.g. bxbam_file is "file1.h5"
    python.assign("bxbam_file", bxbam_file)
    python.exec("store = pd.HDFStore(bxbam_file)")
    python.get("print('loaded bxBam file')")
    python.get("store")
}


#' @name get_bmates
#' @title get_bmates
#' @description
#'
#' Reads R vector of BX strings, outputs GRanges or data.table
#'
#' Function used to output all corresponding bam fields and GRanges via data.table
#' @param reads R vector of BX strings, e.g. bmates_vec <- c("bx1", "bx2", "bx3", "bx4")
#' @return Corresponding fields in a bam file in data.table, GRanges, or GRangesList
#' @export

get_bmates <- function(reads){
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(hdf_key, where='BX in @reads'")
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
#' @param reads R vector of BX strings, e.g. bmates_vec <- c("bx1", "bx2", "bx3", "bx4")
#' @return Corresponding fields in a bam file in data.table, GRanges, or GRangesList
#' @export

get_qmates <- function(reads){
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(hdf_key, where='QNAME in @reads'")
    df <- python.get("df")
    df <- data.table(df)
    return(df)
}





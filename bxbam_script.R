######
#%install.packages("rPython")
#%install.packages("h5")    necessary???

library(rPython)
library(data.table)


# Load python script
python.load("bxbam_functions.py")

python.exec("import pandas as pd")
python.exec("import numpy as np")
python.exec("import tables")

hdf_key = "bam_fields"

bxbam_file <- readline("Prompt: Enter bxbam file name")

load_bxbam <- function(bxbam_file){     # e.g. bxbam_file is "file1.h5"
    python.assign("bxbam_file", bxbam_file)
    python.exec("store = pd.HDFStore(bxbam_file)")
    python.get("print('loaded bxBam file')")
    python.get("store")
}

#
# Users input R vector of bmates
# bmates_vec <- c("bx1", "bx2", "bx3", "bx4")
#

get_bmates <- function(reads){
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(hdf_key, where='BX in @reads'")
    df <- python.get("df")
    df <- data.table(df)
}

#
# users write/read


get_qmates <- function(reads){
    python.assign("reads", reads)
    pand_func_bmate <- python.exec("df = store.select(hdf_key, where='QNAME in @reads'")
    df <- python.get("df")
    df <- data.table(df)
}





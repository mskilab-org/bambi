import subprocess
import pandas as pd
# import tables    # users must go through downloading this and HDF5 software---put into documentation
import numpy as np
#import numexpr as ne
# also, *must* have Cython installed
import os

bam_path = "/gpfs/commons/home/biederstedte-934/evan_projects/HUMAN_1000Genomes_hs37d5_genomic_RPE1-hTERT_X-33.dupmarked.bam"
bxbam_path = "/gpfs/commons/home/biederstedte-934/evan_projects/HUMAN_1000Genomes_hs37d5_genomic_RPE1-hTERT_X-33.bxbam.h5"

bxbame_name = "HUMAN_1000Genomes_hs37d5_genomic_RPE1-hTERT_X-33.bxbam.h5"

store = pd.HDFStore(bxbame_name)
columns_to_index = [col1, col2]   # list columns to index now, only
col1 = "QNAME"
col2 = "BX"
hdf_key = "default_key"


task = subprocess.Popen("samtools view "+ bam_path, shell=True,  stdout=subprocess.PIPE)
for i, line in enumerate(l.decode(errors='ignore') for l in task.stdout):  # decode binary
    line_split = line.split()           # split line, make list
    mandatory_field = line_split[:11]   # mandatory fields, "values"
    rest_fields = line_split[11:]
    mandatory_field_keys = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
    beginning = '  '.join('{}:{}'.format(k, v) for k, v in zip(mandatory_field_keys, mandatory_field))  # create tsv string
    rest_data = '  '.join(rest_fields)
    total_lines = beginning + rest_data
    # now, parse with awk
    #result = subprocess.check_output(["awk", "-v", "cols=columns.txt", "-f", "bxbam.awk", total_lines])
    result = subprocess.check_output(["awk", "-v", "OFS='\t'", "-f", "bam5.awk", total_lines])
    # result = subprocess.Popen("awk -v OFS='\t' -f bam5.awk "+ total_lines, shell=True,  stdout=subprocess.PIPE)
    store.append(hdf_key, result, data_columns=columns_to_index, index=False)

# index data columns in HDFStore
store.create_table_index(hdf_key, columns=columns_to_index , optlevel=9, kind='full')
store.close()

# awk -v OFS='\t' -f bam5.awk stdout


# It may be best to offer desired order
bxbam_default_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX"]
# then use ` df = df[bxbam_default_columns]
# if MD? Check whether setting parameter in `create_bxbam()` is best
md_default_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX", "MD"]


                        
# functions
             
                        
def get_qmates(param1):
    pd.store(bxbam_path).select("QNAME"==param1)

def get_bmates(param1):
    pd.store(bxbam_path).select("BX"==param1)
                        
def get_field(field, param1):
    field = str(field)
    pd.store(bxbam_path).select(field==param1)
                        
def create_bxbam(bam_path, bxbam_path, column1, columns2):
    columns_to_index = [column1, column2]
    # above code
                        
                        
# finish wrapping into an R package `bxbam` that can use the below functions
                        # library(bxbam)

######
#%install.packages("rPython")
#%install.packages("h5")

library(rPython)
                        

# Load python script
python.load("bxbam.py")   # i.e. this script

bxbam_store <-h5file(bxbam)

bam_path <- readline("Enter bam path")
bxbam_path <- readline("Enter desired location/path for bxbam")
           
column1 <- "QNAME"
column2 <- readline("Enter field to index, e.g. BX")   # default option?
                        
get_bmates <- function(reads, "BX" ) {
    pand_func_bmate <- python.exec(bxbamstore.select("BX"))
    python.call("pand_func_bmate", param1)
}

get_qmates <- function(bxbamstore, param1, param2){
    pand_func_qmate <- python.exec(bxbamstore.select("QNAME")
    python.call("pand_func_qmate", param1)
}

get_field <- function(reads, param1) {
    pand_func <- python.exec(bxbamstore.select(param1))
    python.call("pand_func", param1)
}

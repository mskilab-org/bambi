import subprocess
import pandas as pd
import tables    # users must go through downloading this and HDF5 software---put into documentation
import numpy as np
import numexpr as ne
# also, *must* have Cython installed
import os

bam_path = "pathname"   # users set pathname here, bam_path = input("path name to bam =" )
bxbam_path = "pathname"   # users set pathname here, , bxbam_path = input("path name to bxbam =" ) # use bam_path directory?
columns_to_index = [col1, col2]   # list columns to index now, only
col1 = "QNAME"
col2 = "BX"
hdf_key = "default_key"

task = subprocess.Popen(("samtools view"+str(bam_path), shell=True,  stdout=subprocess.PIPE)
while True:
    line = task.stdout.readline().decode("latin-1")    # annoying---cannot process by byte, only by line
    if len(line) == 0 and task.poll() != None: break
rc = task.wait()
                        
                        
                        
# Set all headers first by running through bam, line by line ---- possibly too expensive; check performance again
headers=set()
with open(bam_path) as file:
    for line in file:
        for record in line.split('\t'):
            head,_,datum=record.partition(":")
            headers.add(head)
            # sort
            bam_columns=sorted(headers, key=lambda e: int(e.partition('_')[2]))

# It may be best to offer desired order
bxbam_default_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX"]
# then use ` df = df[bxbam_default_columns]
# if MD? Check whether setting parameter in `create_bxbam()` is best
md_default_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX", "MD"]


task = subprocess.Popen(("bxbam.awk -F ", shell=True,  stdout=subprocess.PIPE)
while True:   # cols.xt is a flag
    line = task.stdout.readline().decode("latin-1")
    if len(line) == 0 and task.poll() != None: break
    for chunk in pd.read_table(bxbam_path, header=None, sep='\t', chunksize=10**6):
                        # place chunks into a dataframe or HDF
        store.append(hdf_key, chunk, data_columns=columns_to_index, index=False)   # not indexing now
                        
# benchmarking idea: uses internal logic of Pandas to deal with missing columns in a different order
# this is possibly more efficient that awk!
"""
#df=pd.DataFrame()
#with open(bam_path) as file:
    #for i, line in enumerate(file):
         #line_data=pd.DataFrame({k.strip():v.strip() 
             for k,_,v in (e.partition(':') 
                 for e in line.split('\t'))}, index=[i])
#store.append(hdf_key, chunk, data_columns=columns_to_index , index=False)
store.create_table_index(hdf_key, columns=columns_to_index , optlevel=9, kind='full')
store.close()
"""
                        

# index data columns in HDFStore
store.create_table_index(hdf_key, columns=columns_to_index, optlevel=9, kind='full')
store.close()
                        
                        
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

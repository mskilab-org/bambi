import subprocess
import pandas as pd
import tables
import numpy as np

# R script
# users first set directory
# setwd(dir)
# h5_name <- readline("string")

h5_name <- read("")

hdf_key = "bam_fields"    # from create_bxbam.py


python.assign("hdf5name", h5name )            # hdf5name is the python variable
python.exec("store = pd.HDFStore(hdf5name)")  # load into pandas HDFStore
python.get("store")                           # returns variable 'store'


def read_bxbam(filename):
    store = pd.HDFStore(filename)
    return store 

                        
def get_qmates(param1):
    store.select("QNAME"==param1)

def get_bmates(param1):
    store.select("BX"==param1)
                        
def get_field(field, param1):
    field = str(field)
    store.select(field==param1)
                        
def create_bxbam(bam_path, bxbam_path, column1, columns2):
    columns_to_index = [column1, column2]

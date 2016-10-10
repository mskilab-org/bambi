import subprocess
import pandas as pd
import tables    # users must go through downloading this and HDF5 software---put into documentation
from tables import *
import numpy as np
#import numexpr as ne # also, *must* have Cython installed

import os
os.chdir("/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-1-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/")

bam_path = "/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-1-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam"

bxbam_name = "IMI-1-10X-chromium_phased_possorted_bam5.h5"

store = pd.HDFStore(bxbam_name, mode='w')
col1 = "QNAME"
col2 = "BX"
columns_to_index = [col1, col2]   # list columns to index now, only

hdf_key = "bxbam_key"


bxbam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX"]



task = subprocess.Popen("samtools view "+ bam_path, shell=True,  stdout=subprocess.PIPE)
for i, line in enumerate(l.decode(errors='ignore') for l in task.stdout):  # decode binary
    toks = line.split()   # split on all whitespace
    newline = dict(zip(bxbam_columns, toks))  # dictionary of 11 mandatory fields
    toks12 = toks[12:]  # rest of lines
    newline1 = dict(item.split(':', 1) for item in toks12) # 2nd dictionary for rest of field:value's
    merged_dict = newline.copy() # first dictionary...best method for Python pre-3.5
    merged_dict.update(newline1) # ...merged with second dictionary
    dicts_into_df = pd.DataFrame(merged_dict, index=[i])  # place in pandas dictionary
    dicts_into_df = dicts_into_df[bxbam_columns]
    dicts_into_df.QNAME = dicts_into_df.QNAME.astype(str)
    dicts_into_df.FLAG = dicts_into_df.FLAG.astype(int)
    dicts_into_df.RNAME = dicts_into_df.RNAME.astype(str)
    dicts_into_df.POS = dicts_into_df.POS.astype(int)
    dicts_into_df.MAPQ = dicts_into_df.MAPQ.astype(int)
    dicts_into_df.CIGAR = dicts_into_df.CIGAR.astype(str)
    dicts_into_df.RNEXT = dicts_into_df.RNEXT.astype(str)
    dicts_into_df.PNEXT = dicts_into_df.PNEXT.astype(int)
    dicts_into_df.TLEN = dicts_into_df.TLEN.astype(int)
    dicts_into_df.SEQ = dicts_into_df.SEQ.astype(str)
    dicts_into_df.QUAL = dicts_into_df.QUAL.astype(str)
    dicts_into_df.BX = dicts_into_df.BX.astype(str)
    dicts_into_df = dicts_into_df.set_index(["QNAME", "BX"])
    store.append(hdf_key, dicts_into_df,  data_columns=bxbam_columns, index=False, min_itemsize={'BX':175, 'QNAME':275, 'CIGAR':100, 'RNAME':100, 'RNEXT':50})

# index data columns in HDFStore
store.create_table_index(hdf_key, columns=columns_to_index , optlevel=9, kind='full')
store.close()


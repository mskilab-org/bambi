import subprocess
import pandas as pd
import tables
from tables import *
import numpy as np
import argparse

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input", help = "input filename")
parser.add_argument("output", help = "output filename")
parser.add_argument("param1", help = "yada yada")
parser.add_argument("param2", help = "yada yada")
parser.add_argument("final", help = "blah blah yada yada")
args = parser.parse_args()
input = args.input
output = args.output
param1 = args.input
param2 = args.output
final = args.final

parser = argparse.ArgumentParser()
parser.add_argument("input", help = "input bam file path/filename")
parser.add_argument("output", help = "output bxbam file path/filename")
parser.add_argument("index_fields", nargs = "*", help = "bam fields to index")
parser.add_argument("MD", help="store MD tag into bxBam file")
args = parser.parse_args()
bam_path = args.input
bxbam_name = args.output
fields = args.index_fields  # at the moment, users are limited by five fields which they can index
fields = fields.split(',')
field1 = index_fields[0]
field2 = index_fields[1]
if len(fields) == 3:
    field3 = index_fields[2]
elif len(fields) == 4:
    field4 = index_fields[3]
elif len(fields) == 5:
    field5 = index_fields[4]


# chose fields to index
# col1 = "QNAME"
# col2 = "BX"
# columns_to_index = [col1, col2]   # list columns to index now, only

if len(fields) == 2:
    columns_to_index = [field1, field2]
elif len(fields) == 3:
    columns_to_index = [field1, field2, field3]
elif len(fields) == 4:
    columns_to_index = [field1, field2, field3, field4]
elif len(fields) == 5:
    columns_to_index = [field1, field2, field3, field4, field5]



hdf_key = "bam_fields"

bxbam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX"]

bxbam_columns_MD = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "BX", "MD"]

# create bxbam description

class Barcoded_bam_MD(IsDescription):
        QNAME = StringCol(64)
        FLAG  = Int16Col()
        RNAME = StringCol(64)
        POS   = Int32Col()
        MAPQ  = Int8Col()
        CIGAR = StringCol(64)
        RNEXT = StringCol(64)
        PNEXT = Int32Col()
        TLEN  = Int32Col()
        SEG   =  StringCol(256)
        QUAL  = StringCol(256)
        BX    = StringCol(64)
        MD    = StringCol(64)
    
class Barcoded_bam(IsDescription):
        QNAME = StringCol(64)
        FLAG  = Int16Col()
        RNAME = StringCol(64)
        POS   = Int32Col()
        MAPQ  = Int8Col()
        CIGAR = StringCol(64)
        RNEXT = StringCol(64)
        PNEXT = Int32Col()
        TLEN  = Int32Col()
        SEG   =  StringCol(256)
        QUAL  = StringCol(256)
        BX    = StringCol(64)

if args.MD:  # if users set the first class, use that
    bxbam_description = Barcoded_bam_MD()
    bxbam_field_columns = bxbam_columns_MD
else:
    bxbam_description = Barcoded_bam()
    bxbam_field_columns = bxbam_columns


# Open a file in write mode
h5file = open_file(bxbam_name, mode = "w")

# Create group
group = h5file.create_group("/", "bam_table")

table = h5file.create_table(group, hdf_key, bxbam_description, "table of 10x bam field values")

bxbamrow = table.row

task = subprocess.Popen("samtools view "+ bam_path, shell=True,  stdout=subprocess.PIPE)
for i, line in enumerate(l.decode(errors='ignore') for l in task.stdout):  # decode binary
    toks = line.split()   # split on all whitespace
    newline = dict(zip(bxbam_field_columns, toks))  # dictionary of 11 mandatory fields
    toks12 = toks[12:]  # rest of lines
    newline1 = dict(item.split(':', 1) for item in toks12) # 2nd dictionary for rest of field:value's
    merged_dict = newline.copy() # first dictionary...best method for Python pre-3.5
    merged_dict.update(newline1) # ...merged with second dictionary
    bxbam["QNAME"] = merged_dict["QNAME"] if merged_dict["QNAME"] != "" else "NaN"
    bxbam["FLAG"]  = merged_dict["FLAG"] if merged_dict["FLAG"] != "" else "NaN"
    bxbam["RNAME"]  = merged_dict["RNAME"] if merged_dict["RNAME"] != "" else "NaN"
    bxbam["POS"]  = merged_dict["POS"] if merged_dict["POS"] != "" else "NaN"
    bxbam["MAPQ"]  = merged_dict["MAPQ"] if merged_dict["MAPQ"] != "" else "NaN"
    bxbam["CIGAR"]  = merged_dict["CIGAR"] if merged_dict["CIGAR"] != "" else "NaN"
    bxbam["RNEXT"]  = merged_dict["RNEXT"] if merged_dict["RNEXT"] != "" else "NaN"
    bxbam["PNEXT"]  = merged_dict["PNEXT"] if merged_dict["PNEXT"] != "" else "NaN"
    bxbam["TLEN"]  = merged_dict["TLEN"] if merged_dict["TLEN"] != "" else "NaN"
    bxbam["SEG"]  = merged_dict["SEQ"] if merged_dict["SEQ"] != "" else "NaN"
    bxbam["QUAL"]  = merged_dict["QUAL"] if merged_dict["QUAL"] != "" else "NaN"
    bxbam["BX"]  = merged_dict["BX"] if merged_dict["BX"] != "" else "NaN"
    # This injects the Record values
    bxbamrow.append()
    # creates indices
    if len(fields) == 2:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
    elif len(fields) == 3:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
    elif len(fields) == 4:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
    elif len(fields) == 5:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
        table.cols.field5.create_index()
    # Flush the table buffers
    table.flush()

# Finally, close the file (this also will flush all the remaining buffers!)
h5file.close()




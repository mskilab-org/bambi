import subprocess
import pandas as pd
import tables
from tables import *
import numpy as np
import argparse

def h5extension(astr):             # make sure users use 'h5' as bxbam extension
    if not astr.endswith(".h5"):
        astr += ".h5"
    return astr

parser = argparse.ArgumentParser()
parser.add_argument("input", help = "input bam file path/filename")
parser.add_argument("output", help = "output bxbam file path/filename", type=h5extension)
parser.add_argument("index_fields", nargs = "*", help = "bam fields to index")
parser.add_argument("MD", help="store MD tag into bxBam file")
args = parser.parse_args()
bam_path = args.input
bxbam_name = args.output
fields = args.index_fields  # at the moment, users are limited by five fields which they can index
fields = fields[0].split(',')
if  len(fields) == 1:
    field1 = index_fields[0]
elif len(fields) == 2:
    field2 = index_fields[1]
elif len(fields) == 3:
    field3 = index_fields[2]
elif len(fields) == 4:
    field4 = index_fields[3]
elif len(fields) == 5:
    field5 = index_fields[4]
elif len(fields) == 6:
    field6 = index_fields[5]
elif len(fields) == 7:
    field7 = index_fields[6]
elif len(fields) == 8:
    field8 = index_fields[7]
elif len(fields) == 9:
    field9 = index_fields[8]
elif len(fields) == 10:        # we arbitrarily stop at 10 fields to index---you really shouldn't use more---one or two is optimal
    field10 = index_fields[9]



if not args.input:
    paser.error("input bam not given, please provide via -input '<bamname>'")
if not args.output:
    paser.error("output bxbam not named, please provide via -output '<bxbamname.h5>'")


print("input bam entered: "+ bam_path)
print("bxbam file output: "+ bxbam_name)
if len(fields) == 0:
    print("no indexing")
elif len(fields) == 1:
    print("indexing on: "+ field1)
elif len(fields) == 2:
    print("indexing on: "+ field1 + ", " + field2)
elif len(fields) == 3:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3)
elif len(fields) == 4:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4)
elif len(fields) == 5:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5)
elif len(fields) == 6:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field5)
elif len(fields) == 7:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field5 + ", " + field5)
elif len(fields) == 8:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field5 + ", " + field5 + ", " + field5)
elif len(fields) == 9:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field5 + ", " + field5 + ", " + field5 + ", " + field5)
elif len(fields) == 10:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field5 + ", " + field5 + ", " + field5 + ", " + field5 + ", " + field5)


# chose fields to index
# col1 = "QNAME"
# col2 = "BX"
# columns_to_index = [col1, col2]   # list columns to index now, only

if  len(fields) == 1:
    columns_to_index = [field1]
elif len(fields) == 2:
    columns_to_index = [field1, field2]
elif len(fields) == 3:
    columns_to_index = [field1, field2, field3]
elif len(fields) == 4:
    columns_to_index = [field1, field2, field3, field4]
elif len(fields) == 5:
    columns_to_index = [field1, field2, field3, field4, field5]
elif len(fields) == 6:
    columns_to_index = [field1, field2, field3, field4, field5, field6]
elif len(fields) == 7:
    columns_to_index = [field1, field2, field3, field4, field5, field6, field7]
elif len(fields) == 8:
    columns_to_index = [field1, field2, field3, field4, field5, field6, field7, field8]
elif len(fields) == 9:
    columns_to_index = [field1, field2, field3, field4, field5, field6, field7, field8, field9]
elif len(fields) == 10:
    columns_to_index = [field1, field2, field3, field4, field5, field6, field7, field8, field9, field10]





hdf_key = "bxbam_default_key"

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
    bxbam["FLAG"]  = merged_dict["FLAG"]  if merged_dict["FLAG"] != ""  else "NaN"
    bxbam["RNAME"] = merged_dict["RNAME"] if merged_dict["RNAME"] != "" else "NaN"
    bxbam["POS"]   = merged_dict["POS"]   if merged_dict["POS"] != ""   else "NaN"
    bxbam["MAPQ"]  = merged_dict["MAPQ"]  if merged_dict["MAPQ"] != ""  else "NaN"
    bxbam["CIGAR"] = merged_dict["CIGAR"] if merged_dict["CIGAR"] != "" else "NaN"
    bxbam["RNEXT"] = merged_dict["RNEXT"] if merged_dict["RNEXT"] != "" else "NaN"
    bxbam["PNEXT"] = merged_dict["PNEXT"] if merged_dict["PNEXT"] != "" else "NaN"
    bxbam["TLEN"]  = merged_dict["TLEN"]  if merged_dict["TLEN"] != ""  else "NaN"
    bxbam["SEG"]   = merged_dict["SEQ"]   if merged_dict["SEQ"] != ""   else "NaN"
    bxbam["QUAL"]  = merged_dict["QUAL"]  if merged_dict["QUAL"] != ""  else "NaN"
    bxbam["BX"]    = merged_dict["BX"]    if merged_dict["BX"] != ""    else "NaN"
    # This injects the Record values
    bxbamrow.append()
    # creates indices
    if  len(fields) == 1:
        table.cols.field1.create_index()
    elif len(fields) == 2:
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
    elif len(fields) == 6:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
        table.cols.field5.create_index()
        table.cols.field6.create_index()
    elif len(fields) == 7:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
        table.cols.field5.create_index()
        table.cols.field6.create_index()
        table.cols.field7.create_index()
    elif len(fields) == 8:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
        table.cols.field5.create_index()
        table.cols.field6.create_index()
        table.cols.field7.create_index()
        table.cols.field8.create_index()
    elif len(fields) == 9:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
        table.cols.field5.create_index()
        table.cols.field6.create_index()
        table.cols.field7.create_index()
        table.cols.field8.create_index()
        table.cols.field9.create_index()
    elif len(fields) == 10:
        table.cols.field1.create_index()
        table.cols.field2.create_index()
        table.cols.field3.create_index()
        table.cols.field4.create_index()
        table.cols.field5.create_index()
        table.cols.field6.create_index()
        table.cols.field7.create_index()
        table.cols.field8.create_index()
        table.cols.field9.create_index()
        table.cols.field10.create_index()
    # Flush the table buffers
    table.flush()

# Finally, close the file (this also will flush all the remaining buffers!)
h5file.close()




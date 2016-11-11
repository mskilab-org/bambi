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
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field6)
elif len(fields) == 7:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field6 + ", " + field7)
elif len(fields) == 8:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field6 + ", " + field7 + ", " + field8)
elif len(fields) == 9:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field6 + ", " + field7 + ", " + field8 + ", " + field9)
elif len(fields) == 10:
    print("indexing on: "+ field1 + ", " + field2 + ", " + field3 + ", " + field4 + ", " + field5 + ", " + field6 + ", " + field7 + ", " + field8 + ", " + field9 + ", " + field10)


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
        SEQ   =  StringCol(256)
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
        SEQ   =  StringCol(256)
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
    bxbam["SEQ"]   = merged_dict["SEQ"]   if merged_dict["SEQ"] != ""   else "NaN"
    bxbam["QUAL"]  = merged_dict["QUAL"]  if merged_dict["QUAL"] != ""  else "NaN"
    bxbam["BX"]    = merged_dict["BX"]    if merged_dict["BX"] != ""    else "NaN"
    # This injects the Record values
    bxbamrow.append()
    # Flush the table buffers
    table.flush()

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

# Finally, close the file (this also will flush all the remaining buffers!)
h5file.close()



















# must do `$ module remove R/3.3.0
# must do `$  module load R/.3.3.1-milab
# CORRECT R version: /gpfs/commons/groups/imielinski_lab/Software/R-3.3.1/bin/R ; python 3.4.3


/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-1c-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/


import pandas as pd
import numpy as np
import tables
from tables import *




bxbamfile ="/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-3-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/IMI3oct24_bxbam.h5"

bamfile = "/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-3-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bam"

baifile = "/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-3-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/phased_possorted_bam.bai"

import os
os.chdir("/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-3-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/")



h5f = tables.open_file("IMI3oct24_bxbam.h5")
tb1 = h5f.get_node("/bam_table/bam_fields")  # entire key
BX1 = tb1.read_where('BX==b"Z:GATAACCAGCGCTCAC-1"')  # takes like 4 seconds, type=numpy.ndarray


BX1 = tb1.read_where(\\'BX==b"Z:GATAACCAGCGCTCAC-1"\\')
python.exec('BX1 = tb1.read_where(\'BX==b"Z:GATAACCAGCGCTCAC-1"\')')

#len(BX1) = 920
# REPLICABLE RESULTS!!!!
#TypeError: You can't use Python's standard boolean operators in NumExpr expressions. You should use their bitwise counterparts instead: '&' instead of 'and', '|' instead of 'or', and '~' instead of 'not'.
BBB = tb1.read_where('(BX==b"Z:GATAACCAGCGCTCAC-1") |(BX== b"Z:CGATCGGGTACAGATA-1")')
# len(BBB)  # 2491


df1 = pd.DataFrame(BX1)   # shaped (920, 12)

"""
    >>> df1.columns
    Index(['BX', 'CIGAR', 'FLAG', 'MAPQ', 'PNEXT', 'POS', 'QNAME', 'QUAL', 'RNAME',
    'RNEXT', 'SEG', 'TLEN'],
    dtype='object')
    
"""


library(rPython)
library(data.table)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)



setwd("/gpfs/commons/home/biederstedte-934/data/Maciejowski2015/10x/IMI-3-10X-chromium/PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/fork0/files/")

python.exec("import pandas as pd")
python.exec("import numpy as np")
python.exec("import tables")
python.exec("from tables import *")


bxbamfile = "IMI3oct24_bxbam.h5"

python.assign("bxbamfile", bxbamfile)

python.exec("h5f = tables.open_file(bxbamfile)")
python.exec('tb1 = h5f.get_node("/bam_table/bam_fields")')

python.exec('BX1 = tb1.read_where(\'BX==b"Z:GATAACCAGCGCTCAC-1"\')')

python.exec("df1 = pd.DataFrame(BX1)")

# one approach --- hard code each column---probably quickest---# can only run this once per dataframe!!!

python.exec('df1.BX = df1.BX.str.decode("utf-8")')
python.exec('df1.CIGAR = df1.CIGAR.str.decode("utf-8")')
python.exec('df1.QNAME = df1.QNAME.str.decode("utf-8")')
python.exec('df1.QUAL = df1.QUAL.str.decode("utf-8")')
python.exec('df1.RNAME = df1.RNAME.str.decode("utf-8")')
python.exec('df1.RNEXT = df1.RNEXT.str.decode("utf-8")')
python.exec('df1.SEQ = df1.SEG.str.decode("utf-8")')      # error here!~!~!!!!!!


# convert to Python lists...

python.exec(" BX_list = df1.BX.tolist()")
python.exec(" CIGAR_list = df1.CIGAR.tolist()")
python.exec(" FLAG_list = df1.FLAG.tolist()")
python.exec(" MAPQ_list = df1.MAPQ.tolist()")
python.exec(" PNEXT_list = df1.PNEXT.tolist()")
python.exec(" POS_list = df1.POS.tolist()")
python.exec(" QNAME_list = df1.QNAME.tolist()")
python.exec(" QUAL_list = df1.QUAL.tolist()")
python.exec(" RNAME_list = df1.RNAME.tolist()")
python.exec(" RNEXT_list = df1.RNEXT.tolist()")
python.exec(" SEQ_list = df1.SEQ.tolist()")
python.exec(" TLEN_list = df1.TLEN.tolist()")


# ...then access as R vectors

BX_vec <- python.get("BX_list")
CIGAR_vec <- python.get("CIGAR_list")
FLAG_vec <- python.get("FLAG_list")
MAPQ_vec <- python.get("MAPQ_list")
PNEXT_vec <- python.get("PNEXT_list")
POS_vec <- python.get("POS_list")
QNAME_vec <- python.get("QNAME_list")
QUAL_vec <- python.get("QUAL_list")
RNAME_vec <- python.get("RNAME_list")
RNEXT_vec <- python.get("RNEXT_list")
SEQ_vec <- python.get("SEQ_list")
TLEN_vec <- python.get("TLEN_list")


qname <- python.get("QNAME_list")
flag  <- python.get("FLAG_list")
rname <- python.get("RNAME_list")
pos <- python.get("POS_list")
mapq <- python.get("MAPQ_list")
cigar <- python.get("CIGAR_list")
rnext <- python.get("RNEXT_list")
pnext <- python.get("PNEXT_list")
tlen <- python.get("TLEN_list")
seq <- python.get("SEQ_list")
qual <- python.get("QUAL_list")
bx  <- python.get("BX_list")



countCigar <- function(cigar) {
    
    cigar.vals <- unlist(strsplit(cigar, "\\d+"))
    cigar.lens <- strsplit(cigar, "[A-Z]")
    lens <- nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] <- 1
    
    cigar.lens <- as.numeric(unlist(cigar.lens))
    cigar.vals <- cigar.vals[cigar.vals != ""]
    repr       <- rep(seq_along(cigar), lens)
    dt         <- data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d      <- dt["D",][, sum(lens), by=group]
    smr.i      <- dt["I",][, sum(lens), by=group]
    smr.m      <- dt["M",][, sum(lens), by=group]
    smr.s      <- dt["S",][, sum(lens), by=group]
    
    out <- matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group,1] <- smr.d$V1
    out[smr.i$group,2] <- smr.i$V1
    out[smr.m$group,3] <- smr.m$V1
    out[smr.s$group,4] <- smr.s$V1
    colnames(out) <- c('D','I','M','S')
    
    return(out)
}

# create data.table

out <- data.table(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, bx)

cigs <- countCigar(out$cigar)   # columns  D I M S
# Warning message:
# In countCigar(out$cigar) : NAs introduced by coercion
out$pos2 <- out$pos + cigs[, "M"]
# M" positions in the alignment are encoded in the CIGAR; if add to the BAM pos you get the end position of the interval
unmatched = is.na(out$pos)

if (any(unmatched))
{
    out$pos[unmatched] = 1
    out$pos2[unmatched] = 0
    out$strand[unmatched] = '*'
}



gr.fields = c('rname', 'strand', 'pos', 'pos2');
vals = out[, setdiff(names(out), gr.fields), with=FALSE]
# vals ==  this is the metadata of the outgoing GRanges object which is pulled from the data.frame "out"


out <- GRanges(out$rname, IRanges(out$pos, pmax(0, out$pos2-1)), strand = out$strand, seqlengths = seqlengths(intervals))
values(out) <- vals
                             
                             
out <- data.table(seqnames=out$rname, start=out$pos, end= pmax(out$pos2-1, 0), strand=out$strand)
val <- data.table(vals)
out <- cbind(out, val)

# Warning message:
# In data.table::data.table(...) :
# Item 2 is of size 12 but maximum size is 920 (recycled leaving remainder of 8 items)
#

return(out)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intervals = NULL

if (length(intervals)==0)
    intervals = NULL
    
if (is.null(intervals))
    intervals = gr




@@@@@@@@@@@@@@@@@@@@@@@@@


rr <- IRanges(out$pos, out$pos2)
if (!'strand' %in% colnames(out))
    outt$strand <- '*'


sf <- factor(out$strand, levels=c('+', '-', '*'))
ff <- factor(out$seqnames, levels=unique(out$rnames))

# we now have a dataframed


b'Z:GATAACCAGCGCTCAC-1'

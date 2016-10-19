# bxBam
indexed file format for barcoded BAMs with API for converting and accessing alignment records

# dependencies

Python >= 2.7 or Python 3.x

HDF5 >= 1.8.4 (recommend >=1.8.15 , HDF5 v1.10 not supported)

PyTables (https://github.com/PyTables/PyTables)

NumPy >= 1.8.1

Numexpr >= 2.5.2

Cython >= 0.21

Users install Python dependencies with `pip install -r requirements.txt`

# Rough API: 

(1) Given a 10x bam file, create bxBam file with `python create_bxbam.py`. First at the path to your bam in the script, as well as the bxbam name. Users can chose which fields to index---the default is on `BX` and `QUAL` to query bmates against qmates. There exists a flag `--MD` to include the `MD` tag, if this is included (e.g. if user has generated MD tags with `samtools calmd <bam-file>`). Otherwise, the default will include 12 fields: the 11 mandatory bam fields and BX. 

     python indexing_issues_bxbam.py --input="/pathname/path/file.bam" --output "pathname/bxbam.h5" --index_fields="QNAME","BX" 

Question to self: does it matter if they do not use a .h5 extension? If so, let's just write one in.

Optional column is --MD

(2) Run

install.packages("rPython")   This uses the Python version currently loaded---important for `library('rPython')`




# bxBam
indexed file format for barcoded BAMs with API for converting and accessing alignment records

# dependencies

Python >= 2.7 or Python 3.x
HDF5 >= 1.8.4 (recommend >=1.8.15 , HDF5 v1.10 not supported)
PyTables (https://github.com/PyTables/PyTables)
NumPy >= 1.8.1
Numexpr >= 2.5.2
Cython >= 0.21

# Rough API: 
 
Data location for 10x:
Files found here `/gpfs/commons/groups/imielinski_lab/data/10x/bam3/`

Step (1)
Please copy the `.bam` file into your local `~/Projects` directory. If we somehow change the original files, that's not good. 

Gzip the `.bam` file
`samtools view file.bam | gzip > out.gz`
where `out.gz` is the output file. 

Step (2)
load python3.5.1 via the virtual environment:  
    (a)` cd /gpfs/commons/groups/imielinski_lab/lib/python-3/bin/`
    (b) `source activate`
    PS: check you have the correct version with `which python`

Now, execute the following command: 
    `python bxBam.py --filename="out.gz" --HDF_name="file1.h5" --columns_to_index = ["QNAME", "BX"]`

There are two flags: filename and columns (There the default key is `df_key`---feel free to change with the flag `--key==<enter value>`)

This takes time!! I recommend you use the linux command `screen` or `nohup COMMAND &`.


Step (3):
Within python (e.g. a jupyter notebook) import the bxBAM table (the HDF5 file), e.g. 

    import pandas as pd
    store = pd.HDFStore("file1.h5") 
  
Step (4): Perform queries

e.g. query non-local QNAMEs

    df = store.select("df_key", 
       where="QNAME='E00438:46:HMTV7CCXX:3:2202:4330:8956' or QNAME='E00438:46:HMTV7CCXX:8:1216:21917:20893'")

Note the AND/OR keywords

    df2 = store.select("df_key", 
      where=QNAME='E00438:46:HMTV7CCXX:8:1216:21917:20893' and "BX='yyyyyy'")
      
Step (5)
Decide what to do with these queries 

You could save as .csv files

    df.to_csv("query1.csv", index=False)

You could concatenate these together

    df1  # first query
    df2  # second query
    df_together = pd.concat([df1, d2])
    df_togther("all_queries.csv", index=False)



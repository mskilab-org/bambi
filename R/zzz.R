# zzz.R for bxbam

## @import rPython
.onLoad <- function(libname = find.package("bxbam"), pkgname = "bxbam") {
    ## rPython::python.exec("import pandas as pd")
    ## rPython::python.exec("import numpy as np")
    ## rPython::python.exec("import tables")
    ## rPython::python.exec("from tables import *")
    ## rPython::python.exec("sessions = {}")
    ## rPython::python.exec("queries = {}")

    ## # one approach --- hard code each column---probably quickest---# can only run this once per dataframe!!!
    
    ## ## define run_query function
    ## cmd = sprintf('def run_query(tb, qstring):
    ##     tmp = tb.read_where(qstring)
    ##     df1 = pd.DataFrame(tmp)
    ##     df1.BX = df1.BX.str.decode("utf-8")
    ##     df1.CIGAR = df1.CIGAR.str.decode("utf-8")
    ##     df1.QNAME = df1.QNAME.str.decode("utf-8")
    ##     df1.QUAL = df1.QUAL.str.decode("utf-8")
    ##     df1.RNAME = df1.RNAME.str.decode("utf-8")
    ##     df1.RNEXT = df1.RNEXT.str.decode("utf-8")
    ##     df1.SEQ = df1.SEQ.str.decode("utf-8")   # error here!~!~!!!!!!
    ##     return(df1)')
    ## rPython::python.exec(cmd)
}

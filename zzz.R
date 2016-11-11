# zzz.R for bxbam

#' @import rPython

.onLoad <- function(libname = find.package("bxbam"), pkgname = "bxbam") {
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
}

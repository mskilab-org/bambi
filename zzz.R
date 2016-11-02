# zzz.R for bxbam

#' @import rPython

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to my package") # not necessary, delete
    python.exec("import pandas as pd")
    python.exec("import numpy as np")
    python.exec("import tables")
    python.exec("from tables import *")
}


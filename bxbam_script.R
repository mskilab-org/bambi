######
#%install.packages("rPython")
#%install.packages("h5")

library(rPython)
                        

# Load python script
python.load("bxbam_functions.py")  

bxbam_store <-h5file(bxbam)

bam_path <- readline("Enter bam path")
bxbam_path <- readline("Enter desired location/path for bxbam")
           
column1 <- "QNAME"
column2 <- readline("Enter field to index, e.g. BX")   # default option?
                        
get_bmates <- function(reads, "BX" ) {
    pand_func_bmate <- python.exec(bxbam_store.select("BX"))
    python.call("pand_func_bmate", param1)
}

get_qmates <- function(bxbamstore, param1, param2){
    pand_func_qmate <- python.exec(bxbam_store.select("QNAME")
    python.call("pand_func_qmate", param1)
}

#
# get_field <- function(reads, param1) {
#    pand_func <- python.exec(bxbamstore.select(param1))
#    python.call("pand_func", param1)
# }
#

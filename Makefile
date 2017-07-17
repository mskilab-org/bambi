PKG_BASE = $(shell pwd)/bamdb/

all:
	cd $(PKG_BASE) && $(MAKE)

##include code to compile htslib and ck and lmdb here too b/c the user may only want to use the comand line features.

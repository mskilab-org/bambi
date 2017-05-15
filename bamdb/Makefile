CC		= clang
CFLAGS	= -Wall -g -std=gnu99

PKG_BASE = $(shell pwd)/../
SOURCES	:= src/bam_api.c src/bam_lmdb.c
LIB		:= -lhts -lm -llmdb -lpthread
INC		:= -I include -I$(PKG_BASE)/ck/include/ -I$(PKG_BASE)htslib/ -I$(PKG_BASE)lmdb/libraries/liblmdb/

all:
	@mkdir -p bin/
	$(CC) $(CFLAGS) $(SOURCES) bamdb.c -L$(PKG_BASE)lmdb/libraries/liblmdb/ -L$(PKG_BASE)/htslib/ -L$(PKG_BASE)/ck/ -o bin/bxbam $(LIB) $(INC)

clean:
	rm -f bin/bxbam

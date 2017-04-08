CC		= gcc
CFLAGS	= -Wall -g -std=c11 -fPIC
PKG_BASE = $(shell pwd)/../


##SOURCES	:= src/bam_api.c src/bam_sqlite.c src/bam_lmdb.c
##LIB		:= -lhts -lm -lsqlite3 -llmdb
INC		:= -Iinclude -I./ -I$(PKG_BASE)/htslib/ -I$(PKG_BASE)/lmdb/libraries/liblmdb/ -I$(PKG_BASE)/sqlite3/

all:
	@mkdir -p bin/

	$(CC) $(CFLAGS) -c bamdb.c -o bin/bamdb.o $(INC)
	$(CC) $(CFLAGS) -c src/bam_api.c -o bin/bam_api.o $(INC)
	$(CC) $(CFLAGS) -c src/bam_lmdb.c -o bin/bam_lmdb.o $(INC)
	$(CC) $(CFLAGS) -c src/bam_sqlite.c -o bin/bam_sqlite.o $(INC)

	ar rs bin/libbamdb.a bin/bam_sqlite.o bin/bam_api.o bin/bamdb.o

##	$(CC) $(CFLAGS) $(SOURCES) bamdb.c -L /usr/local/lib/ -o bin/bamdb $(LIB) $(INC)

clean:
##	rm -f bin/bamdb
	rm -f bin/bamdb.o
	rm -f bin/bam_api.o
	rm -f bin/bam_lmdb.o
	rm -f bin/bam_sqlite.o
	rm -f bin/libbamdb.a

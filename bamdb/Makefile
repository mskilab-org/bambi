CC		= gcc
CFLAGS	= -Wall -g -std=c11 -fPIC -fpic
PKG_BASE = $(shell pwd)/../


#SOURCES	:= src/bam_api.c src/bam_sqlite.c src/bam_lmdb.c
#LIB		:= -lhts -lsqlite3 -llmdb -lm

##INC		:= -Iinclude -I./ -I/gpfs/commons/home/knagdimov/temp/htslib/ -I/gpfs/commons/home/knagdimov/temp/lmdb/libraries/liblmdb/ -I/gpfs/commons/home/knagdimov/temp/sqlite3/

##PKG_LIBS 	:= -L/gpfs/commons/home/knagdimov/temp/htslib -L/gpfs/commons/home/knagdimov/temp/lmdb/libraries/liblmdb/ -L/gpfs/commons/home/knagdimov/temp/sqlite3/ -lhts -llmdb -lsqlite3

##INC		:= -Iinclude/ -I./ -I/gpfs/commons/home/knagdimov/temp/git/Rbamdb/htslib/ -I/gpfs/commons/home/knagdimov/temp/git/Rbamdb/lmdb/libraries/liblmdb/ -I/gpfs/commons/home/knagdimov/temp/git/Rbamdb/sqlite3/

INC		:= -Iinclude/ -I./ -I$(PKG_BASE)/htslib/ -I$(PKG_BASE)lmdb/libraries/liblmdb/ -I$(PKG_BASE)sqlite3/



#PKG_LIBS	:= -L/gpfs/commons/home/knagdimov/temp/git/Rbamdb//htslib -L/gpfs/commons/home/knagdimov/temp/git/Rbamdb/lmdb/libraries/liblmdb/ -L/gpfs/commons/home/knagdimov/temp/git/Rbamdb/sqlite3 -lhts -llmdb -lsqlite3

all:
	@mkdir -p bin
        ##$(CC) -c $(CFLAGS) $(SOURCES) bamdb.c $(PKG_LIBS) -o bin/bamdb $(LIB) $(INC)

	$(CC) $(CFLAGS) -c bamdb.c -o bin/bamdb.o $(INC)
	$(CC) $(CFLAGS) -c src/bam_api.c -o bin/bam_api.o $(INC)
	$(CC) $(CFLAGS) -c src/bam_lmdb.c -o bin/bam_lmdb.o $(INC)
	$(CC) $(CFLAGS) -c src/bam_sqlite.c -o bin/bam_sqlite.o $(INC)

	ar rs bin/libbamdb.a bin/bam_sqlite.o bin/bam_api.o bin/bamdb.o

#	ar rs bin/libbamdb.a bin/bamdb.o bin/bam_api.o bin/bam_lmdb.o bin/bam_sqlite.o

	#$(CC) $(CFLAGS) $(SOURCES) bamdb.c $(PKG_LIBS) -o bin/bamdb $(LIB) $(INC)
	#gcc -o file.cgi bin/bamdb.o bin/bam_api.o bin/bam_lmdb.o bin/bam_sqlite.o $(PKG_LIBS) -pthread

clean:
	rm -f bin/bamdb.o
	rm -f bin/bam_api.o
	rm -f bin/bam_lmdb.o
	rm -f bin/bam_sqlite.o
	rm -f bin/libbamdb.a

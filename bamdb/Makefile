CC		= gcc
CFLAGS	= -Wall -g -std=c99

SOURCES	:= src/bam_api.c src/bam_sqlite.c src/bam_lmdb.c
LIB		:= -lhts -lm -lsqlite3 -llmdb
INC		:= -I include

all:
	@mkdir -p bin/
	$(CC) $(CFLAGS) $(SOURCES) bamdb.c -L /usr/local/lib/ -o bin/bamdb $(LIB) $(INC)

clean:
	rm -f bin/bamdb

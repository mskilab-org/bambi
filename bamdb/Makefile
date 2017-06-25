.PHONY: clean install uninstall

CC		= clang
CFLAGS	= -Wall -g -std=gnu99

BUILD_DIR = build
INCLUDE_DIR = include
SRC_DIR = src
PREFIX = sdfsdf

TARGET = bxbam

OBJECTS = $(patsubst %.c, %.o, $(wildcard $(SRC_DIR)/*.c))
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)

%.o: %.c $(HEADERS)
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) -c $< -o $(patsubst src/%, $(BUILD_DIR)/%, $@)

BUILT_OBJECTS = $(patsubst src/%, $(BUILD_DIR)/%, $(OBJECTS))

PKG_BASE = $(shell pwd)/../
SOURCES	:= src/bam_api.c src/bam_lmdb.c
LIB		:= -lhts -lm -llmdb -lpthread -lz
INC		:= -I/ -Iinclude -I$(PKG_BASE)/ck/include/ -I$(PKG_BASE)htslib/ -I$(PKG_BASE)lmdb/libraries/liblmdb/

all:
	@mkdir -p bin/
	$(CC) -c bamdb.c $(INC) -o bin/bamdb.o
	$(CC) -c src/bam_api.c $(INC) -o bin/bam_api.o
	$(CC) -c src/bam_lmdb.c $(INC) -o bin/bam_lmdb.o

	ar cr bin/libbamdb.so bin/bamdb.o bin/bam_api.o bin/bam_lmdb.o

	$(CC) $(CFLAGS) $(SOURCES) bamdb.c -L$(PKG_BASE)lmdb/libraries/liblmdb/ -L$(PKG_BASE)/htslib/ -L$(PKG_BASE)/ck/ -o $(PKG_BASE)/bxbam $(LIB) $(INC)

clean:
	rm -f bin/bxbam
	rm -f $(PKG_BASE)/bxbam
	rm -f bin/bamdb.o
	rm -f bin/bam_api.o
	rm -f bin/bam_lmdb.o
	rm -f bin/lbamdb.a

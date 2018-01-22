
temporary repo for bxBam R API
original here: github.com/mskilab/bxbam

## Dependencies

liblmdb-dev

libck-dev

libhts-dev

libz-dev

## MacOS users using brew

brew install lmdb

brew install concurrencykit

brew install htslib  // recommend the latest version via https://github.com/samtools/htslib

brew install zlib    // mac OS already provides this software




For compilers to find this software you may need to set:

LDFLAGS:  -L/usr/local/opt/zlib/lib

CPPFLAGS: -I/usr/local/opt/zlib/include

https://github.com/concurrencykit/ck


## for installing bamdb on mskilab (with cmake)

cd bamdb-master

module load htslib/1.5

module load cmake/3.8.0

cmake .

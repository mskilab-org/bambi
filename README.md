
## bbami

instantiate
check LMDB indicies
grab_bx
grab_cb
grab_ub
grab_tag 

## bamdb

Program: bamdb (Software for indexing and querying 10x BAMs)
Version: 0.9 (using htslib 1.5)

Usage:   bamdb index <required-flag>

Required flags:
WGS            index QNAME and BX tags for 10x WGS BAM
single-cell    index QNAME, CB, and UB tags for single-cell BAM


e.g. create bamdb index for WGS on BX and QNAME field

bamdb index WGS example.bam

e.g. create bamdb index for WGS on BX and QNAME field, and also MD and

bamdb index WGS example.bam  --tag=MD --tag=BZ




### Dependencies

liblmdb-dev

libck-dev

libhts-dev

libz-dev

### MacOS users using brew

brew install lmdb

brew install concurrencykit

brew install htslib  // recommend the latest version via https://github.com/samtools/htslib

brew install zlib    // mac OS already provides this software




For compilers to find this software you may need to set:

LDFLAGS:  -L/usr/local/opt/zlib/lib

CPPFLAGS: -I/usr/local/opt/zlib/include

https://github.com/concurrencykit/ck


### for installing bamdb on mskilab (with cmake)

cd bamdb-master

module load htslib/1.5

module load cmake/3.8.0

cmake .

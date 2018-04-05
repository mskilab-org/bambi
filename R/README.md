# bambi
R package for querying 10x WGS and single-cell BAMs

## Dependencies
Currently, first install `bamdb` package as follows

    git clone https://github.com/D-Lo/bamdb
    cd bamdb
    mkdir build
    cd build
    cmake .. 
    make

If one then uses `sudo make install`, the library/headers with be placed in `/usr/local/include` and gcc will find them for this R package. If not, link to headers in the `Makevars`

On the R side, one needs:

    RCpp
    mskilab/bamUtils
    mskilab/gUtils


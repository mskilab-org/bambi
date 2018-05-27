# Minimum of autoconfg 2.69
AC_PREREQ(2.69)

AC_INIT([bambi], [0.0.1])

# We're using C++11
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX
AC_GNU_SOURCE

AC_SUBST([PKG_CXXFLAGS],["-std=c++11"])

base_dir=$(pwd)

AC_MSG_NOTICE([Building zlib library])
(
    cd "src/zlib"; \
    mkdir -p "build"; \
    cd "build"; \
    cmake .. ; \
    make;
)

AC_MSG_NOTICE([zlib build complete])



AC_MSG_NOTICE([Building ck library])
(
    cd "src/ck"; \
    ./configure; \
    make;
)

AC_MSG_NOTICE([ck build complete])



AC_MSG_NOTICE([Building htslib library])
(
    cd "src/htslib"; \
    make;
)

AC_MSG_NOTICE([htslib build complete])


AC_MSG_NOTICE([Building lmdb library])
(
    cd "src/lmdb/libraries/liblmdb"; \
    make;
)

AC_MSG_NOTICE([lmdb build complete])


AC_MSG_NOTICE([Building bamdb library])
(
    cd "src/bamdb"; \
    mkdir -p "build"; \
    cd "build"; \
    cmake .. -DSTATIC_BAMDB_BUILD=ON; \
    make;
)

AC_MSG_NOTICE([bamdb build complete])



## TODO: Dynamically find htslib
BAMDB_INCL_FLAGS="-I$base_dir/src/zlib/include -I$base_dir/src/ck/include/ -I$base_dir/src/htslib/ -I$base_dir/src/htslib/htslib -I/usr/include/htslib -I/usr/include/htslib/htslib -I$base_dir/src/bamdb/include -I$base_dir/src/lmdb/libraries/liblmdb/"
BAMDB_LIB="$base_dir/src/bamdb/build/libbamdb.a"

bamdb_include="${BAMDB_INCL_FLAGS}"
bamdb_libs="-L$base_dir/src/zlib/build -lz -L$base_dir/src/ck -L$base_dir/src/lmdb/libraries/liblmdb/ -llmdb -L$base_dir/src/htslib/ -L$base_dir/src/htslib/htslib -L/usr/include/htslib  -lhts ${BAMDB_LIB}"

AC_SUBST([PKG_CFLAGS],["${PKG_CFLAGS} $bamdb_include"])
AC_SUBST([PKG_CXXFLAGS],["${PKG_CXXFLAGS} $bamdb_include"])
AC_SUBST([PKG_LIBS],["${PKG_LIBS} $bamdb_libs"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

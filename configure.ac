# Minimum of autoconfg 2.69
AC_PREREQ(2.69)

AC_INIT([bambi], [0.0.1])

# We're using C++11
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX
AC_GNU_SOURCE
AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
AC_SUBST([PKG_CXXFLAGS],["-std=c++11"])

base_dir=$(pwd)

AC_MSG_NOTICE([Updating bamdb submodule])
$(git submodule update --init)

AC_MSG_NOTICE([Building bamdb library])
$(
    cd "src/bamdb"; \
    mkdir -p "build"; \
    cd "build"; \
    cmake .. -DSTATIC_BAMDB_BUILD=ON; \
    make;
)

AC_MSG_NOTICE([bamdb build complete])

# TODO: Dynamically find htslib
BAMDB_INCL_FLAGS="-I$base_dir/src/bamdb/include -I/usr/include/htslib/ -I/usr/local/include/htslib/"
BAMDB_LIB="$base_dir/src/bamdb/build/libbamdb.a"

bamdb_include="${BAMDB_INCL_FLAGS}"
bamdb_libs="${BAMDB_LIB} -llmdb -lhts"

AC_SUBST([PKG_CFLAGS],["${PKG_CFLAGS} $bamdb_include"])
AC_SUBST([PKG_CXXFLAGS],["${PKG_CXXFLAGS} $bamdb_include"])
AC_SUBST([PKG_LIBS],["${PKG_LIBS} $bamdb_libs"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

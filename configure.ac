# Minimum of autoconfg 2.69
AC_PREREQ(2.69)

AC_INIT([bambi], [0.0.1])

# We're using C++11
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX
AC_GNU_SOURCE
AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
AC_SUBST([PKG_CFLAGS],["-fPIC"])
AC_SUBST([PKG_CXXFLAGS],["-fPIC -std=c++11"])

base_dir=$(pwd)

# Create directory to build libraries
$(mkdir -p "$base_dir/src/build";)

AC_MSG_NOTICE([Building lmdb])
$(
    cd "$base_dir/src/lmdb/libraries/liblmdb"; \
    make clean; \
    make XCFLAGS=-fPIC;
)

AC_MSG_NOTICE([lmdb build complete])

LMDB_INCL_PATH="$base_dir/src/lmdb/libraries/liblmdb"
LMDB_INCL_FLAGS="-I${LMDB_INCL_PATH}"
LMDB_LIB="$base_dir/src/lmdb/libraries/liblmdb/liblmdb.a"
LMDB_CMAKE_FLAGS="-DLMDB_INCLUDE_DIR=${LMDB_INCL_PATH} -DLMDB_LIBRARY=${LMDB_LIB}"

AC_MSG_NOTICE([Building bamdb library])
$(
    cd "$base_dir/src/build"; \
	cmake "$base_dir/src/bamdb/" -DSTATIC_BAMDB_BUILD=ON -DBUILD_BAMDB_WRITER=OFF ${LMDB_CMAKE_FLAGS}; \
    make;
)

AC_MSG_NOTICE([bamdb build complete])

BAMDB_INCL_FLAGS="-I$base_dir/src/bamdb/include -I/usr/include/htslib/ -I/usr/local/include/htslib/"
BAMDB_LIB="$base_dir/src/build/libbamdb.a"

bambi_include="${LMDB_INCL_FLAGS} ${BAMDB_INCL_FLAGS}"
bambi_libs="${BAMDB_LIB} ${LMDB_LIB} -lhts"

AC_SUBST([PKG_CFLAGS],["${PKG_CFLAGS} $bambi_include"])
AC_SUBST([PKG_CXXFLAGS],["${PKG_CXXFLAGS} $bambi_include"])
AC_SUBST([PKG_LIBS],["${PKG_LIBS} $bambi_libs"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

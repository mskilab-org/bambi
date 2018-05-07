# - Try to find LMDB
# Once done, this will define
#
#  LMDB_FOUND - system has LMDB
#  LMDB_INCLUDE_DIRS - the LMDB include directories
#  LMDB_LIBRARIES - link these to use LMDB

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(LMDB_PKGCONF LMDB)

# Include dir
find_path(LMDB_INCLUDE_DIR
  NAMES lmdb.h
  PATHS ${LMDB_PKGCONF_INCLUDE_DIRS}
)

# Find the library itself
find_library(LMDB_LIBRARY
  NAMES lmdb
  PATHS ${LMDB_PKGCONF_LIBRARY_DIRS}
)

libfind_process(LMDB)

# - Try to find Concurrency Kit
# Once done, this will define
#
#  CK_FOUND - system has CK
#  CK_INCLUDE_DIRS - the CK include directories
#  CK_LIBRARIES - link these to use CK

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(CK_PKGCONF CK)

# Include dir
find_path(CK_INCLUDE_DIR
  NAMES ck_fifo.h
  PATHS ${CK_PKGCONF_INCLUDE_DIRS}
)

# Find the library itself
find_library(CK_LIBRARY
  NAMES ck
  PATHS ${CK_PKGCONF_LIBRARY_DIRS}
)

libfind_process(CK)

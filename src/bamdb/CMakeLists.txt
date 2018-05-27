cmake_minimum_required(VERSION 2.8)

project(bamdb)

set(CMAKE_C_FLAGS "-Wall -g -std=gnu99 -fPIC")

option(STATIC_BAMDB_BUILD "STATIC_BAMDB_BUILD" OFF)

# External dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

## set variables for include/libraries
## CK_INCLUDE_DIR=<not found>
## CK_LIBRARY=<not found>
## LMDB_INCLUDE_DIR=<not found>
## LMDB_LIBRARY=<not found>
## HTSlib_INCLUDE_DIR=<not found>
## HTSlib_LIBRARY=<not found>
## 
## 

## set(HTSlib_INCLUDE_DIR )

message('here is the package base directory, PROJECT_SOURCE_DIR')
message('$PROJECT_SOURCE_DIR')
message('${PROJECT_SOURCE_DIR}')
message('looking for CMAKE_CURRENT_LIST_DIR')
message("${CMAKE_CURRENT_LIST_DIR}")
message("upper subdirectory")
message("${CMAKE_CURRENT_LIST_DIR}/../")

set(CK_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/../ck")
set(CK_LIBRARY "${CMAKE_CURRENT_LIST_DIR}/../ck/include")


find_package(CK REQUIRED)
find_package(LMDB REQUIRED)
find_package(HTSlib REQUIRED)


set(LIBS
  ${CK_LIBRARIES}
  ${HTSlib_LIBRARIES}
  ${LMDB_LIBRARIES}
  m
  pthread
)

include_directories(include)
file(GLOB SOURCES "src/*.c")

# Include dependency headers
include_directories(
  ${CK_INCLUDE_DIRS}
  ${HTSlib_INCLUDE_DIRS}
  ${LMDB_INCLUDE_DIRS}
)

# Create bamdb executable
add_executable(bamdb ${SOURCES})
target_link_libraries(bamdb ${LIBS})

# Create bamdb library
if(STATIC_BAMDB_BUILD)
  add_library(libbamdb STATIC ${SOURCES})
else(STATIC_BAMDB_BUILD)
  add_library(libbamdb SHARED ${SOURCES})
endif(STATIC_BAMDB_BUILD)

target_link_libraries(libbamdb ${LIBS})
set_target_properties(libbamdb PROPERTIES OUTPUT_NAME "bamdb")

# Default install path. Can be overridden with DESTDIR
install(TARGETS bamdb RUNTIME DESTINATION bin)
install(DIRECTORY include/
  DESTINATION include/bamdb)
install(TARGETS libbamdb
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

cmake_minimum_required(VERSION 2.8.12)

project(bamdb)

set(CMAKE_C_FLAGS "-Wall -g -std=gnu99 -fPIC")

option(STATIC_BAMDB_BUILD "STATIC_BAMDB_BUILD" OFF)
option(BUILD_BAMDB_WRITER "BUILD_BAMDB_WRITER" ON)


# External dependencies
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(HTSlib REQUIRED)
find_package(LMDB REQUIRED)

set(LIBS
  ${HTSlib_LIBRARIES}
  ${LMDB_LIBRARIES}
  m
  pthread
)

if(BUILD_BAMDB_WRITER)
  add_definitions(-DBUILD_BAMDB_WRITER)
  find_package(CK REQUIRED)
  list(APPEND LIBS ${CK_LIBRARIES})
  include_directories(${CK_INCLUDE_DIRS})
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")
file(GLOB SOURCES "src/*.c")

# Include dependency headers
include_directories(
  ${HTSlib_INCLUDE_DIRS}
  ${LMDB_INCLUDE_DIRS}
)

# Create bamdb library
if(STATIC_BAMDB_BUILD)
  add_library(libbamdb STATIC ${SOURCES})
else(STATIC_BAMDB_BUILD)
  add_library(libbamdb SHARED ${SOURCES})
endif(STATIC_BAMDB_BUILD)

target_link_libraries(libbamdb ${LIBS})
set_target_properties(libbamdb PROPERTIES OUTPUT_NAME "bamdb")

# Create bamdb executable
if(BUILD_BAMDB_WRITER)
  add_executable(bamdb "src/bamdb_main.c")
  target_link_libraries(bamdb ${LIBS} libbamdb)
endif()

# Default install path. Can be overridden with DESTDIR
if(BUILD_BAMDB_WRITER)
  install(TARGETS bamdb RUNTIME DESTINATION bin)
endif()
install(DIRECTORY include/
  DESTINATION include/bamdb)
install(TARGETS libbamdb
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

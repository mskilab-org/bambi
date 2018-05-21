cmake_minimum_required(VERSION 2.8)

project(bamdb)

set(CMAKE_C_FLAGS "-Wall -g -std=gnu99 -fPIC")

# External dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CK REQUIRED)
find_package(HTSlib REQUIRED)
find_package(LMDB REQUIRED)

set(LIBS
  ${CK_LIBRARIES}
  ${HTSlib_LIBRARIES}
  ${LMDB_LIBRARIES}
  m
  pthread
)

include_directories(include)
file(GLOB SOURCES "src/*.c")

# Create bamdb executable
add_executable(bamdb ${SOURCES})
target_link_libraries(bamdb ${LIBS})

# Create bamdb library
add_library(libbamdb SHARED ${SOURCES})
target_link_libraries(libbamdb ${LIBS})
set_target_properties(libbamdb PROPERTIES OUTPUT_NAME "bamdb")

# Default install path. Can be overridden with DESTDIR
install(TARGETS bamdb RUNTIME DESTINATION bin)
install(DIRECTORY include/
  DESTINATION include/bamdb) 
install(TARGETS libbamdb
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

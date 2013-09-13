# CMake build script for the Gecode solver for AMPL.

# Disable -pedantic when building Gecode because it causes errors
# when compiling on Mac with older GCC.
string(REPLACE "-pedantic" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(GECODE_DISABLE_WARNINGS True)
add_subdirectory(lib)
include_directories(lib)

add_library(amplgecode gecode.cc)
target_link_libraries(amplgecode util gecodeminimodel)
add_dependencies(amplgecode arith_h)
include_directories(${GECODE_SOURCE_DIR})

add_executable(gecode main.cc)
target_link_libraries(gecode amplgecode)

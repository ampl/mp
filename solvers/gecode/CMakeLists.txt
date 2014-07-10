# CMake build script for the AMPL/Gecode solver.

# Disable -pedantic when building Gecode because it causes errors
# when compiling on Mac with older GCC.
string(REPLACE "-pedantic" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(GECODE_DISABLE_WARNINGS True)
add_subdirectory(lib)
include_directories(lib)

add_library(amplgecode-static gecode.cc)
target_link_libraries(amplgecode-static util gecodeminimodel)
add_dependencies(amplgecode-static arith_h)
include_directories(${GECODE_SOURCE_DIR})

add_library(amplgecode SHARED ../util/solver-c.cc)
target_link_libraries(amplgecode amplgecode-static)

add_executable(gecode main.cc)
target_link_libraries(gecode amplgecode-static)

# CMake build script for the Gecode solver for AMPL.

add_subdirectory(gecode-src)

include_directories(gecode-src)
add_library(amplgecode gecode.cc)
target_link_libraries(amplgecode util gecodeminimodel)
add_dependencies(amplgecode arith_h)
add_executable(gecode main.cc)
target_link_libraries(gecode amplgecode)


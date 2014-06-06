# CMake build script for the cp library.

add_ampl_library(cp cp.cc)
target_link_libraries(cp util)

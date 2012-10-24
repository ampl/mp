# CMake build script for the driver library that contains
# a boilerplate code required to implement an AMPL driver.

add_library(util expr.cc expr.h)
add_dependencies(util arith_h)
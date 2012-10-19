# CMake build script for the driver library that contains
# a boilerplate code required to implement an AMPL driver.

add_library(util util.cc util.h expr.cc expr.h)

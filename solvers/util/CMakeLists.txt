# CMake build script for the driver library that contains
# a boilerplate code required to implement an AMPL driver.

add_library(util expr.cc expr.h format.cc format.h driver.cc driver.h)
add_dependencies(util arith_h)
target_link_libraries(util ${SOLVER_LIBS})

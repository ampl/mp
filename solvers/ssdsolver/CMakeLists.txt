# CMake build script for AMPL/SSDSolver.

add_library(amplssdsolver-static ssdsolver.cc)
target_link_libraries(amplssdsolver-static util)
add_dependencies(amplssdsolver-static arith_h)

add_library(amplssdsolver SHARED ../util/solver-c.cc)
target_link_libraries(amplssdsolver amplssdsolver-static)

add_executable(ssdsolver main.cc)
target_link_libraries(ssdsolver amplssdsolver-static)

add_ampl_library(ssd ssd.cc)

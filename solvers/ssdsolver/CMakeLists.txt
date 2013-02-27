# CMake build script for the SSD solver for AMPL.

add_library(amplssdsolver ssdsolver.cc)
target_link_libraries(amplssdsolver util)
add_dependencies(amplssdsolver arith_h)
add_executable(ssdsolver main.cc)
target_link_libraries(ssdsolver amplssdsolver)
add_ampl_library(ssd ssd.cc)

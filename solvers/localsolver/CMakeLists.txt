# CMake build script for AMPL/LocalSolver.

include_directories(${LOCALSOLVER_INCLUDE_DIRS})

add_library(ampllocalsolver-static localsolver.cc)
target_link_libraries(ampllocalsolver-static util ${LOCALSOLVER_LIBRARIES})
add_dependencies(ampllocalsolver-static arith_h)

add_library(ampllocalsolver SHARED ../util/solver-c.cc)
target_link_libraries(ampllocalsolver ampllocalsolver-static)

add_executable(localsolver main.cc)
target_link_libraries(localsolver ampllocalsolver-static)

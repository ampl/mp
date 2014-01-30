# CMake build script for the driver library that contains
# a boilerplate code required to implement an AMPL driver.

add_definitions(-DUSE_CLOCK_GETTIME)
add_library(util
  clock.cc clock.h error.cc error.h expr.cc expr.h
  expr-factory.cc expr-factory.h format.cc format.h nlreader.cc nlreader.h
  noncopyable.h os.cc os.h problem.cc problem.h rstparser.cc rstparser.h
  solver.cc solver.h)
add_dependencies(util arith_h)
target_link_libraries(util amplsolver)

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(util ${RT_LIBRARY})
endif ()

# CMake build script for the driver library that contains
# a boilerplate code required to implement an AMPL driver.

include(CheckSymbolExists)
check_symbol_exists(mkstemps stdlib.h HAVE_MKSTEMPS)
if (HAVE_MKSTEMPS)
  add_definitions(-DASL_HAVE_MKSTEMPS=1)
else ()
  set(UTIL_EXTRA_SRC mkstemps.c)
endif ()

add_definitions(-DUSE_CLOCK_GETTIME)
add_library(util
  clock.cc clock.h error.h expr.cc expr.h
  expr-factory.cc expr-factory.h format.cc format.h nl.cc nl.h
  noncopyable.h os.cc os.h problem.cc problem.h rstparser.cc rstparser.h
  solver.cc solver.h ${UTIL_EXTRA_SRC})
add_dependencies(util arith_h)
target_link_libraries(util amplsolver)

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(util ${RT_LIBRARY})
endif ()

# CMake build script for the driver library that contains
# a boilerplate code required to implement an AMPL driver.

add_library(util
  clock.h clock.cc error.h expr.cc expr.h filesystem.cc filesystem.h
  format.cc format.h noncopyable.h problem.cc problem.h solver.cc solver.h)
add_dependencies(util arith_h)
target_link_libraries(util amplsolver)

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(util ${RT_LIBRARY})
endif ()

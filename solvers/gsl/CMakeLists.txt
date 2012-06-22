# CMake build script for the GSL wrapper for AMPL.

find_package(GSL)
if (GSL_FOUND)
  if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE Release)
  endif ()

  add_library(amplgsl SHARED amplgsl.c)
  target_link_libraries(amplgsl ${SOLVER_LIBS} ${GSL_LIBRARIES})
endif ()

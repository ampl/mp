# CMake build script for the Gecode driver for AMPL.

find_package(Gecode)
if (GECODE_FOUND)
  include_directories(${GECODE_INCLUDE_DIRS})
  add_library(amplgecode gecode.cc)
  add_dependencies(amplgecode arith_h)
  add_executable(gecode main.cc)
  target_link_libraries(gecode amplgecode ${SOLVER_LIBS} ${GECODE_LIBRARIES})
endif ()

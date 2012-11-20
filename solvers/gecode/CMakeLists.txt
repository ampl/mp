# CMake build script for the Gecode driver for AMPL.

find_package(Gecode)
if (GECODE_FOUND)
  include_directories(${GECODE_INCLUDE_DIRS})
  add_library(gecodedriver gecode.cc)
  target_link_libraries(gecodedriver util ${SOLVER_LIBS} ${GECODE_LIBRARIES})
  add_dependencies(gecodedriver arith_h)
  add_executable(gecode main.cc)
  target_link_libraries(gecode gecodedriver)
endif ()

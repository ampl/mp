# CMake build script for the Gecode solver for AMPL.

find_package(Gecode)
if (GECODE_FOUND)
  include_directories(${GECODE_INCLUDE_DIRS})
  add_library(amplgecode gecode.cc)
  target_link_libraries(amplgecode util ${SOLVER_LIBS} ${GECODE_LIBRARIES})
  add_dependencies(amplgecode arith_h)
  add_executable(gecode main.cc)
  target_link_libraries(gecode amplgecode)
endif ()

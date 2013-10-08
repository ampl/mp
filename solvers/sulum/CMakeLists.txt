# CMake build script for the Sulum solver for AMPL.

find_package(Sulum)
if (SULUM_FOUND)
  add_definitions(-DSLM_NO_BOOST)
  include_directories(amplsulum ${SULUM_INCLUDE_DIRS})
  add_library(amplsulum sulum.cc)
  add_dependencies(amplsulum arith_h)
  target_link_libraries(amplsulum ${SULUM_LIBRARIES} util)

  add_executable(sulum main.cc)
  target_link_libraries(sulum amplsulum)
endif ()

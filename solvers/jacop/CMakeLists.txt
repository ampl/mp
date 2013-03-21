# CMake build script for the JaCoP solver for AMPL.

find_package(JNI)
if (JNI_FOUND)
  include_directories(${JNI_INCLUDE_DIRS})
  add_library(ampljacop jacop.cc)
  target_link_libraries(ampljacop util ${JNI_LIBRARIES})
  add_dependencies(ampljacop arith_h)
  add_executable(jacop main.cc)
  target_link_libraries(jacop ampljacop)
endif ()

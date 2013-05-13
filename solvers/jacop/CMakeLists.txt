# CMake build script for the JaCoP solver for AMPL.

include_directories(${JNI_INCLUDE_DIRS})
add_library(ampljacop jacop.cc jacop.h java.cc java.h)
target_link_libraries(ampljacop util ${JNI_LIBRARIES})
add_dependencies(ampljacop arith_h)
add_executable(jacop main.cc)
target_link_libraries(jacop ampljacop)
if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
  target_link_libraries(jacop delayimp)
endif ()


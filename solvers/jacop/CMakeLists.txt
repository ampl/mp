# CMake build script for the JaCoP solver for AMPL.

include_directories(${JNI_INCLUDE_DIRS})

# Allow "long long" because it it used by JNI.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAVE_WNO_LONG_LONG)
if (HAVE_WNO_LONG_LONG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
endif ()

add_library(ampljacop jacop.cc jacop.h java.cc java.h)
target_link_libraries(ampljacop util ${JNI_LIBRARIES})
add_dependencies(ampljacop arith_h)

add_executable(jacop main.cc)
target_link_libraries(jacop ampljacop)
if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
  target_link_libraries(jacop delayimp)
endif ()

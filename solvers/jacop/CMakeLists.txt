# CMake build script for the JaCoP solver for AMPL.

# Allow "long long" because it it used by JNI.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAVE_WNO_LONG_LONG)
if (HAVE_WNO_LONG_LONG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
endif ()

add_library(ampljacop_lib jacop.cc jacop.h java.cc java.h)
set_target_properties(ampljacop_lib PROPERTIES OUTPUT_NAME ampljacop)
target_link_libraries(ampljacop_lib util ${JNI_LIBRARIES})
add_dependencies(ampljacop_lib arith_h)
include_directories(${JNI_INCLUDE_DIRS})

add_executable(jacop main.cc)
target_link_libraries(jacop ampljacop_lib)
if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
  target_link_libraries(jacop delayimp)
endif ()

include(UseJava)
set(CMAKE_JAVA_INCLUDE_PATH lib/JaCoP-3.2.jar)
add_jar(ampljacop Interrupter.java SolutionListener.java)

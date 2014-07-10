# CMake build script for the AMPL/JaCoP solver.

# Allow "long long" because it it used by JNI.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAVE_WNO_LONG_LONG)
if (HAVE_WNO_LONG_LONG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
endif ()

if (JACOP_JAR MATCHES "jacop-(.*)\\.jar")
  set(JACOP_VERSION ${CMAKE_MATCH_1})
else ()
  message(FATAL_ERROR "Cannot determine JaCoP version.")
endif ()

add_definitions(-DJACOP_VERSION="${JACOP_VERSION}")

add_library(ampljacop-static jacop.cc jacop.h java.cc java.h)
target_link_libraries(ampljacop-static util ${JNI_LIBRARIES})
add_dependencies(ampljacop-static arith_h)
include_directories(${JNI_INCLUDE_DIRS})

add_library(ampljacop-shared SHARED ../util/solver-c.cc)
set_target_properties(ampljacop-shared PROPERTIES OUTPUT_NAME ampljacop)
target_link_libraries(ampljacop-shared ampljacop-static)

add_executable(jacop main.cc)
target_link_libraries(jacop ampljacop-static)
if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
  target_link_libraries(jacop delayimp)
endif ()

include(UseJava)
set(CMAKE_JAVA_INCLUDE_PATH lib/${JACOP_JAR})
add_jar(ampljacop Interrupter.java SolutionListener.java)

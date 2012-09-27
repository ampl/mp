# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.6)
project(AMPL)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE Release)
endif ()

# Set SOLVER_LIBS to point to the AMPL solver library and its dependencies.
include_directories(.)
set(SOLVER_LIBS amplsolver)
if (NOT WIN32)
  set(SOLVER_LIBS ${SOLVER_LIBS} dl)
endif ()

if (NOT MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-missing-field-initializers")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
endif ()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++0x HAS_STD_CPP0X_FLAG)
if (HAS_STD_CPP0X_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif ()

if (EXISTS solvers/amplgsl/gsl/CMakeLists.txt)
  set(HAVE_GSL 1)
endif ()

add_subdirectory(solvers)

enable_testing()
add_subdirectory(tests)

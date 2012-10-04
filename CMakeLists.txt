# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.6)
project(AMPL)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set SOLVER_LIBS to point to the AMPL solver library and its dependencies.
include_directories(. ${CMAKE_CURRENT_BINARY_DIR}
  solvers ${CMAKE_CURRENT_BINARY_DIR}/solvers)
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

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/amplgsl/gsl/CMakeLists.txt)
  set(HAVE_GSL 1)
endif ()

add_subdirectory(solvers)
add_subdirectory(tables)

enable_testing()
add_subdirectory(tests)

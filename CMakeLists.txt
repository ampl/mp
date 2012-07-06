# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.6)
project(AMPL)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find the AMPL solver library.
include_directories(.)
if (NOT WIN32)
  find_path(SOLVER_LIB_DIR amplsolver.a solvers)
  set(SOLVER_LIBS "${SOLVER_LIB_DIR}/amplsolver.a" dl)
elseif (NOT SOLVER_LIBS)
  set(SOLVER_LIBS "${CMAKE_CURRENT_SOURCE_DIR}/solvers/amplsolv.lib")
endif ()

if (NOT MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-missing-field-initializers")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -std=c++0x")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
endif ()

add_subdirectory(solvers/concert)
add_subdirectory(solvers/gsl)

enable_testing()
add_subdirectory(tests)

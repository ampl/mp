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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
endif ()

set(GSL_DIR gsl-1.15)
set(GSL_LIBRARY_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/solvers/amplgsl/${GSL_DIR}/lib")

# Build GSL ourselves with -fPIC to make possible linking of static GSL
# libraries into a dynamic AMPLGSL library.
add_custom_command(OUTPUT solvers/amplgsl/${GSL_DIR}
  COMMAND ./get-gsl WORKING_DIRECTORY solvers/amplgsl)

add_custom_target(gsl-libs ALL DEPENDS solvers/amplgsl/${GSL_DIR})

add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION
  "${GSL_LIBRARY_DIR}/libgsl.a")
add_dependencies(gsl gsl-libs)

add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION
  "${GSL_LIBRARY_DIR}/libgslcblas.a")

add_subdirectory(solvers/amplgsl)
add_subdirectory(solvers/concert)

enable_testing()
add_subdirectory(tests)

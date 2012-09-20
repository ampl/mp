# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.6)
project(AMPL)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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
else ()
  # Allow "long long" even without c++0x because it is used by CPLEX.
  # -Wno-long-long should go after -pedantic.
  check_cxx_compiler_flag(-Wno-long-long HAS_WNO_LONG_LONG_FLAG)
  if (HAS_WNO_LONG_LONG_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
  endif ()
endif ()

if (EXISTS solvers/amplgsl/gsl/configure)
  set(HAVE_GSL 1)

  # IMPORTED GSL libraries are listed here to be accessible both from
  # amplgsl and tests. We build GSL ourselves with -fPIC to make possible
  # linking of static GSL libraries into a dynamic amplgsl library.

  set(GSL_LIBRARY_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/solvers/amplgsl/gsl/lib")

  add_custom_target(build-gsl DEPENDS "${GSL_LIBRARY_DIR}/libgsl.a")
  add_custom_command(OUTPUT "${GSL_LIBRARY_DIR}/libgsl.a"
    COMMAND ./build-gsl WORKING_DIRECTORY solvers/amplgsl)

  add_library(gsl STATIC IMPORTED)
  set_property(TARGET gsl PROPERTY IMPORTED_LOCATION
    "${GSL_LIBRARY_DIR}/libgsl.a")
  add_dependencies(gsl build-gsl)

  add_library(gslcblas STATIC IMPORTED)
  set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION
    "${GSL_LIBRARY_DIR}/libgslcblas.a")
endif ()

add_subdirectory(solvers)

enable_testing()
add_subdirectory(tests)

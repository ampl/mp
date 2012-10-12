# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.6)

# Set the default CMAKE_BUILD_TYPE to Release.
# This should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif ()

# Use static MSVC runtime.
set(CMAKE_USER_MAKE_RULES_OVERRIDE
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cxx_flag_overrides.cmake)

project(AMPL)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

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

include(CheckCCompilerFlag)
check_c_compiler_flag(-fPIC HAVE_FPIC_FLAG)

# Adds a shared AMPL library which by convention doesn't have any prefix
# and has a suffix ".dll".
macro(add_ampl_library name)
  add_library(${name} SHARED ${ARGN})
  set_target_properties(${name} PROPERTIES PREFIX "")
  set_target_properties(${name} PROPERTIES SUFFIX ".dll")
endmacro()

# Change compiler flags to use the DLL runtime - applies to all the binaries
# in the current directory.
macro(use_dll_runtime)
  if (MSVC)
    foreach (var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      if (${var} MATCHES "${match}")
        string(REPLACE "/MT" "/MD" ${var} "${${var}}")
      endif ()
    endforeach ()
  endif ()
endmacro ()

add_subdirectory(solvers)
add_subdirectory(tables)

enable_testing()
add_subdirectory(tests)

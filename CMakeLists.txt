# Use CMake 3.0 as it can be easily obtained and used without superuser
# privileges even if the system has an older version of CMake installed.
cmake_minimum_required (VERSION 3.0)

# Set the default CMAKE_BUILD_TYPE to Release.
# This should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if (NOT CMAKE_BUILD_TYPE)
  set(none "None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used)")
  set(options "${none} Debug Release RelWithDebInfo MinSizeRel")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: ${options}.")
endif ()

# Use static MSVC runtime.
# This should be done before the project command.
set(CMAKE_USER_MAKE_RULES_OVERRIDE
  ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
  ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake/cxx_flag_overrides.cmake)

project(MP)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(MP_CMAKE_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MP_CMAKE_MODULE_DIR})

# Enable C++11.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++11 HAVE_STD_CPP11_FLAG)
if (HAVE_STD_CPP11_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
  check_cxx_compiler_flag(-std=c++0x HAVE_STD_CPP0X_FLAG)
  if (HAVE_STD_CPP0X_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif ()
endif ()

if (NOT MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
else ()
  # Disable useless MSVC warnings suggesting nonportable "secure" alternatives.
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

# Adds a prefix to arguments.
function (add_prefix var prefix)
  set(result ${${var}})
  foreach (arg ${ARGN})
    set(result ${result} "${prefix}${arg}")
  endforeach ()
  set(${var} ${result} PARENT_SCOPE)
endfunction ()

# Adds a shared AMPL library which by convention doesn't have any prefix
# and has a suffix ".dll".
macro(add_ampl_library name)
  add_library(${name} SHARED ${ARGN})
  set_target_properties(${name} PROPERTIES PREFIX "")
  set_target_properties(${name} PROPERTIES SUFFIX ".dll")
  target_include_directories(${name} PUBLIC
    ${PROJECT_SOURCE_DIR}/src/asl/solvers ${PROJECT_BINARY_DIR}/src/asl)
endmacro()

include_directories(include)

include(CheckCCompilerFlag)
set(FPIC_FLAG -fPIC)
check_c_compiler_flag(${FPIC_FLAG} HAVE_FPIC_FLAG)
if (HAVE_FPIC_FLAG)
  add_definitions(${FPIC_FLAG})
endif ()

add_prefix(MP_HEADERS include/mp/
  clock.h error.h format.h nl.h os.h posix.h
  problem-base.h rstparser.h safeint.h)
add_prefix(MP_SOURCES src/ clock.cc format.cc nl.cc os.cc posix.cc rstparser.cc)
add_library(mp ${MP_HEADERS} ${MP_SOURCES})
target_link_libraries(mp asl)

add_subdirectory(src/asl)
add_subdirectory(src/cp)

if (EXISTS ${PROJECT_SOURCE_DIR}/thirdparty/gsl/CMakeLists.txt)
  add_subdirectory(thirdparty/gsl)
  add_subdirectory(src/gsl)
endif ()

enable_testing()
add_subdirectory(test)

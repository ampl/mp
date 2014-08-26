# Use CMake 2.8.12 for target_compile_options.
# It can be easily obtained and used without superuser privileges even
# if the system has an older version of CMake installed.
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})

# Set the path to CMake modules.
set(MP_CMAKE_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MP_CMAKE_MODULE_DIR})

include(init)

project(MP)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Set output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Get optional modules.
set_cache(BUILD "" STRING
  "Comma-separated list of optional modules to build or \"all\""
  "to build all modules.")
if (BUILD)
  if (BUILD STREQUAL all)
    set(MP_MODULES all)
  else ()
    string(REGEX MATCHALL "[^,]+" MP_MODULES "${BUILD}")
  endif ()
endif ()

# Checks if the C compiler supports flag and if it does, adds it to the target.
function (add_c_compiler_flags target)
  foreach (flag ${ARGN})
    check_c_compiler_flag(${flag} supports_flag)
    if (supports_flag)
      target_compile_options(${target} PRIVATE ${flag})
    endif ()
  endforeach ()
endfunction ()

function (add_module name)
  set_property(GLOBAL APPEND PROPERTY
    MP_USED_MODULES ${MP_USED_MODULES} ${name})
endfunction ()

# Sets ${var} to TRUE if module ${name} is enabled.
# Usage:
#   check_module(<name> <var> [EXTERNAL])
# If EXTERNAL is specified, check_module initializes an external dependency
# in thirdparty/${name}.
function(check_module name var)
  get_property(modules GLOBAL PROPERTY MP_USED_MODULES)
  list(FIND modules ${name} module_index)
  if (NOT module_index EQUAL -1)
    set(${var} TRUE PARENT_SCOPE)
    return () # Module already checked.
  endif ()
  add_module(${name})
  list(FIND MP_MODULES ${name} module_index)
  set(${var} FALSE PARENT_SCOPE)
  if (module_index EQUAL -1 AND NOT MP_MODULES STREQUAL all)
    message(STATUS "Module ${name} disabled")
    return ()
  endif ()
  cmake_parse_arguments(check_module EXTERNAL "" "" ${ARGN})
  if (check_module_EXTERNAL)
    # Get a dependency using git.
    find_package(Git)
    if (NOT GIT_FOUND)
      message(SEND_ERROR "Git not found (required by ${name})")
      return ()
    endif ()
    set(dir ${PROJECT_SOURCE_DIR}/thirdparty/${name})
    execute_process(COMMAND
      ${GIT_EXECUTABLE} submodule update --init ${repo} ${dir})
    if (EXISTS ${PROJECT_SOURCE_DIR}/thirdparty/${name}/CMakeLists.txt)
      add_subdirectory(
        ${PROJECT_SOURCE_DIR}/thirdparty/${name}
        ${PROJECT_BINARY_DIR}/thirdparty/build/${name})
    endif ()
  endif ()
  message(STATUS "Module ${name} enabled")
  set(${var} TRUE PARENT_SCOPE)
endfunction()

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

# This is a workaround for ccache bug 8118
# https://bugzilla.samba.org/show_bug.cgi?id=8118
# Should be fixed in ccache version 3.2
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  check_cxx_compiler_flag(-Qunused-arguments HAVE_QUNUSED_ARGUMENTS)
  if (HAVE_QUNUSED_ARGUMENTS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qunused-arguments")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
  endif ()
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
  target_link_libraries(${name} asl)
endmacro()

include_directories(include)

include(CheckCCompilerFlag)
set(FPIC_FLAG -fPIC)
check_c_compiler_flag(${FPIC_FLAG} HAVE_FPIC_FLAG)
if (HAVE_FPIC_FLAG)
  add_definitions(${FPIC_FLAG})
endif ()

set(MP_DATE 20140617)
set(MP_SYSINFO "${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")

add_prefix(MP_HEADERS include/mp/
  clock.h error.h format.h nl.h os.h posix.h
  problem-base.h rstparser.h safeint.h solver.h)
add_prefix(MP_SOURCES src/
  clock.cc format.cc nl.cc os.cc posix.cc problem-base.cc
  rstparser.cc solver.cc)
add_library(mp ${MP_HEADERS} ${MP_SOURCES})
target_compile_definitions(mp PRIVATE
  MP_DATE=${MP_DATE} MP_SYSINFO="${MP_SYSINFO}")
target_include_directories(mp PRIVATE)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "#include <memory>
  int main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
if (HAVE_UNIQUE_PTR)
  target_compile_definitions(mp PUBLIC MP_USE_UNIQUE_PTR)
endif ()

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(mp ${RT_LIBRARY})
endif ()

add_subdirectory(src/asl)
add_subdirectory(src/cp)
add_subdirectory(solvers)

set(GSL_DISABLE_WARNINGS TRUE)
set(GSL_DISABLE_TESTS TRUE) # GSL tests are run separately
check_module(gsl build_gsl EXTERNAL)
if (build_gsl)
  add_subdirectory(src/gsl)
endif ()

# Add a target that generates a file with solver and library versions.
set(AMPL_VERSIONS_FILE versions.txt)
add_custom_target(versions cmake -E remove ${AMPL_VERSIONS_FILE}
  DEPENDS src/asl/tables/amplodbc.c)
foreach (target gsl-info gecode jacop path)
  if (TARGET ${target})
    add_custom_command(TARGET versions POST_BUILD
        COMMAND $<TARGET_FILE:${target}> -v >> ${AMPL_VERSIONS_FILE})
  endif ()
endforeach ()
file(READ src/asl/tables/amplodbc.c amplodbc)
if (amplodbc MATCHES "version ([0-9]+)")
  add_custom_command(TARGET versions POST_BUILD
      COMMAND echo ampltabl ${CMAKE_MATCH_1} >> ${AMPL_VERSIONS_FILE})
endif ()

# Check for invalid module names.
get_property(modules GLOBAL PROPERTY MP_USED_MODULES)
foreach (module ${MP_MODULES})
  if (NOT module STREQUAL all)
    list(FIND modules ${module} module_index)
    if (module_index EQUAL -1)
      message(SEND_ERROR "Unknown module ${module}")
    endif ()
  endif ()
endforeach ()

enable_testing()
add_subdirectory(test)

# Use CMake 2.8.12.2 as it can be easily obtained and used without superuser
# privileges even if the system has an older version of CMake installed.
cmake_minimum_required (VERSION 2.8.12.2)

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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(MP_CMAKE_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MP_CMAKE_MODULE_DIR})

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Get optional modules.
set(BUILD "" CACHE STRING
  "Comma-separated list of optional modules to build or \"all\" to build all modules.")
if (BUILD)
  if (BUILD STREQUAL all)
    set(MP_MODULES all)
  else ()
    string(REGEX MATCHALL "[^,]+" MP_MODULES "${BUILD}")
  endif ()
endif ()

# Sets ${var} to TRUE if module ${name} is enabled.
function (is_module_enabled name var)
  set(MP_USED_MODULES ${MP_USED_MODULES} ${name})
  list(FIND MP_MODULES ${name} module_index)
  if (NOT module_index EQUAL -1 OR MP_MODULES STREQUAL all)
    message(STATUS "Module ${name} enabled")
    set(${var} TRUE PARENT_SCOPE)
  else ()
    message(STATUS "Module ${name} disabled")
    set(${var} FALSE PARENT_SCOPE)
  endif ()
endfunction ()

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

if (CMAKE_GENERATOR MATCHES "Visual Studio")
  # If Microsoft SDK is installed create script run-msbuild.bat that
  # calls SetEnv.cmd to to set up build environment and runs msbuild.
  # It is useful when building Visual Studio projects with the SDK
  # toolchain rather than Visual Studio.
  include(FindSetEnv)
  if (WINSDK_SETENV)
    if (CMAKE_SIZEOF_VOID_P EQUAL 4)
      set(setenv_flag /x86)
    endif()
    # Set FrameworkPathOverride to get rid of MSB3644 warnings.
    file(WRITE run-msbuild.bat "
      call \"${WINSDK_SETENV}\" ${setenv_flag}
      ${CMAKE_MAKE_PROGRAM} -p:FrameworkPathOverride=^\"C:\\Program Files^
\\Reference Assemblies\\Microsoft\\Framework\\.NETFramework\\v4.0^\" %*")
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

add_prefix(MP_HEADERS include/mp/
  clock.h error.h format.h nl.h os.h posix.h
  problem-base.h rstparser.h safeint.h)
add_prefix(MP_SOURCES src/
  clock.cc format.cc nl.cc os.cc posix.cc problem-base.cc rstparser.cc)
add_library(mp ${MP_HEADERS} ${MP_SOURCES})

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(mp ${RT_LIBRARY})
endif ()

add_subdirectory(src/asl)
add_subdirectory(src/cp)
add_subdirectory(solvers)

if (EXISTS ${PROJECT_SOURCE_DIR}/thirdparty/gsl/CMakeLists.txt)
  add_subdirectory(thirdparty/gsl)
  add_subdirectory(src/gsl)
endif ()

# Add a target that generates a file with solver and library versions.
set(AMPL_VERSIONS_FILE versions.txt)
add_custom_target(versions cmake -E remove ${AMPL_VERSIONS_FILE}
  DEPENDS src/asl/tables/amplodbc.c)
foreach (target gen-decls gecode jacop path)
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
foreach (module ${MP_MODULES})
  if (NOT module STREQUAL all)
    list(FIND MP_USED_MODULES ${module} module_index)
    if (module_index EQUAL -1)
      message(SEND_ERROR "Unknown module ${module}")
    endif ()
  endif ()
endforeach ()

enable_testing()
add_subdirectory(test)

# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.8.0)

# Set the default CMAKE_BUILD_TYPE to Release.
# This should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if (NOT CMAKE_BUILD_TYPE)
  set(options
    "None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: ${options}.")
endif ()

# Use static MSVC runtime.
set(CMAKE_USER_MAKE_RULES_OVERRIDE
  ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
  ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake/cxx_flag_overrides.cmake)

project(AMPL)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/support/cmake")

include_directories(. ${CMAKE_CURRENT_BINARY_DIR}
  solvers ${CMAKE_CURRENT_BINARY_DIR}/solvers)

if (NOT MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-missing-field-initializers")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-overlength-strings")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
else ()
  # Disable useless MSVC warnings suggesting nonportable "secure" alternatives.
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

if (APPLE)
  # Ipopt uses thread-local storage which requires OS X Lion or higher target.
  add_definitions(-mmacosx-version-min=10.7)
endif ()

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

include(CheckCCompilerFlag)
check_c_compiler_flag(-fPIC HAVE_FPIC_FLAG)
if (HAVE_FPIC_FLAG)
  add_definitions(-fPIC)
endif ()

# Adds a shared AMPL library which by convention doesn't have any prefix
# and has a suffix ".dll".
macro(add_ampl_library name)
  add_library(${name} SHARED ${ARGN})
  add_dependencies(${name} arith_h)
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

include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
  #include <unordered_map>
  int main() {
    std::unordered_map<int, int> m;
    return 0;
  }" HAVE_UNORDERED_MAP)
if (HAVE_UNORDERED_MAP)
  add_definitions(-DHAVE_UNORDERED_MAP)
endif ()

check_cxx_source_compiles(
  "#include <memory>
  int main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
if (HAVE_UNIQUE_PTR)
  add_definitions(-DHAVE_UNIQUE_PTR)
endif ()

check_cxx_source_compiles("int main() { auto i = 0; }" HAVE_AUTO)

add_subdirectory(doc)
add_subdirectory(solvers)
add_subdirectory(tables)

# Add a target that generates a file with solver and library versions.
set(AMPL_VERSIONS_FILE versions.txt)
add_custom_target(versions cmake -E remove ${AMPL_VERSIONS_FILE}
  DEPENDS tables/amplodbc.c)
foreach (target gen-decls gecode jacop path)
  if (TARGET ${target})
    add_custom_command(TARGET versions POST_BUILD
        COMMAND $<TARGET_FILE:${target}> -v >> ${AMPL_VERSIONS_FILE})
  endif ()
endforeach ()
file(READ tables/amplodbc.c amplodbc)
if (amplodbc MATCHES "version ([0-9]+)")
  add_custom_command(TARGET versions POST_BUILD
      COMMAND echo ampltabl ${CMAKE_MATCH_1} >> ${AMPL_VERSIONS_FILE})
endif ()

enable_testing()
add_subdirectory(tests)

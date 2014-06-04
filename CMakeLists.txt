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

if (APPLE)
  # Unset MACOSX_DEPLOYMENT_TARGET to avoid CMake error
  #   CMAKE_OSX_DEPLOYMENT_TARGET is '...' but CMAKE_OSX_SYSROOT
  #   is not set to a MacOSX SDK with a recognized version.
  set(ENV{MACOSX_DEPLOYMENT_TARGET} )
endif ()

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
check_cxx_compiler_flag(-Qunused-arguments HAVE_QUNUSED_ARGUMENTS)
if (HAVE_QUNUSED_ARGUMENTS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qunused-arguments")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
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
    set(MSBUILD_SETUP "call \"${WINSDK_SETENV}\" ${setenv_flag}")
  endif ()
  # Set FrameworkPathOverride to get rid of MSB3644 warnings.
  file(WRITE run-msbuild.bat "
    ${MSBUILD_SETUP}
    ${CMAKE_MAKE_PROGRAM} -p:FrameworkPathOverride=^\"C:\\Program Files^
\\Reference Assemblies\\Microsoft\\Framework\\.NETFramework\\v4.0^\" %*")
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/amplgsl/gsl/CMakeLists.txt)
  set(HAVE_GSL 1)
  set(get_amplgsl_version $<TARGET_FILE:gen-decls> >> versions)
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/coin/CMakeLists.txt)
  set(HAVE_COIN 1)
  set(get_cbc_version $<TARGET_FILE:cbc-executable> -v -AMPL >> versions)
  if (TARGET ipopt-executable)
    set(get_ipopt_version $<TARGET_FILE:ipopt-executable> -v >> versions)
  endif ()

  # Workaround for f90cache which accept only modified name for GNU fortran
  # compilers in the form gfortran-4.x.
  string(REGEX MATCHALL "[^:]+" search_paths "$ENV{PATH}")
  foreach (search_path ${search_paths})
    file(GLOB paths ${search_path}/gfortran-4.*)
    if (paths)
      list(GET paths 0 CMAKE_Fortran_COMPILER)
      message("Detected Fortran compiler: ${CMAKE_Fortran_COMPILER}")
      break ()
    endif ()
  endforeach ()

  if (WINSDK_SETENV)
    # Create a script that sets up build environment in Fortran subprojects.
    file(WRITE setenv.bat "call \"${WINSDK_SETENV}\"")
    # Extra command line flags to pass to CMake in
    # cmake_add_fortran_subdirectory.
    set(FORTRAN_CMAKE_FLAGS "-DCMAKE_GNUtoMS_VCVARS=../../setenv.bat")
  endif ()
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/gecode/lib/CMakeLists.txt)
  set(HAVE_GECODE 1)
  set(get_gecode_version $<TARGET_FILE:gecode> -v >> versions)
endif ()

file(GLOB JACOP_PATH solvers/jacop/lib/JaCoP-*.jar)
if (JACOP_PATH)
  get_filename_component(JACOP_JAR ${JACOP_PATH} NAME)
  find_package(JNI)
  find_package(Java)
  if (JNI_FOUND AND JAVA_FOUND)
    # Add libmawt.so to JNI_LIBRARIES to avoid linkage errors on Ubuntu 10.04.
    foreach (arch amd64 i386)
      set(LIBMAWT_PATH
        /usr/lib/jvm/java-6-openjdk/jre/lib/${arch}/xawt/libmawt.so)
      if (EXISTS ${LIBMAWT_PATH})
        set(JNI_LIBRARIES ${JNI_LIBRARIES} ${LIBMAWT_PATH})
        break ()
      endif ()
    endforeach ()

    set(HAVE_JACOP 1)
    set(get_jacop_version $<TARGET_FILE:jacop> -v >> versions)
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

add_subdirectory(doc)
add_subdirectory(solvers)
add_subdirectory(tables)

# Add a target that generates versions file listing versions of solvers
# and libraries.
file(READ tables/amplodbc.c amplodbc)
if (amplodbc MATCHES "version ([0-9]+)")
  set(get_ampltabl_version echo ampltabl ${CMAKE_MATCH_1} >> versions)
endif ()
add_custom_target(versions cmake -E remove versions
  COMMAND ${get_amplgsl_version} COMMAND ${get_ampltabl_version}
  COMMAND ${get_cbc_version} COMMAND ${get_ipopt_version}
  COMMAND ${get_gecode_version} COMMAND ${get_jacop_version}
  DEPENDS tables/amplodbc.c)

enable_testing()
add_subdirectory(tests)

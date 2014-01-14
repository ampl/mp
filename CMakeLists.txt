# CMake build script for the AMPL solver library.

cmake_minimum_required(VERSION 2.8.7)

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
  add_definitions(-mmacosx-version-min=10.4)

  # Use libstdc++ instead of libc++ on Mac for compatibility with CPLEX 12.4
  # and earlier.
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")
  set(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} -stdlib=libstdc++")
  set(CMAKE_REQUIRED_FLAGS "-stdlib=libstdc++")

  # The version of libstdc++ available on Mac doesn't support initializer lists.
  add_definitions(-DFMT_USE_INITIALIZER_LIST=0)
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

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/amplgsl/gsl/CMakeLists.txt)
  set(HAVE_GSL 1)
  set(get_amplgsl_version solvers/amplgsl/gen-decls >> versions)
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/cbc/CMakeLists.txt)
  set(HAVE_CBC 1)
  set(get_cbc_version solvers/cbc/bin/cbc -v -AMPL >> versions)
endif ()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/solvers/gecode/lib/CMakeLists.txt)
  set(HAVE_GECODE 1)
  set(get_gecode_version solvers/gecode/gecode -v >> versions)
endif ()

file(GLOB JACOP_PATH solvers/jacop/lib/JaCoP-*.jar)
if (JACOP_PATH)
  get_filename_component(JACOP_JAR ${JACOP_PATH} NAME)
  find_package(JNI)
  find_package(Java)
  if (JNI_FOUND AND JAVA_FOUND)
    set(HAVE_JACOP 1)
    set(get_jacop_version solvers/jacop/jacop -v >> versions)
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

# Adds an AMPL test.
function(add_ampl_test name libs)
  if ("${ARGV2}" STREQUAL "PREFIX")
    list(REMOVE_AT ARGN 0)
    set(prefix "ampl")
  endif ()
  set(target_name "${prefix}${name}")
  add_executable(${target_name} ${ARGN})
  foreach (lib "${libs}")
    target_link_libraries(${target_name} ${lib})
  endforeach()
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  add_test(${name} ${name})
endfunction()

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
  COMMAND ${get_cbc_version} COMMAND ${get_gecode_version}
  COMMAND ${get_jacop_version}
  DEPENDS tables/amplodbc.c)

enable_testing()
add_subdirectory(tests)

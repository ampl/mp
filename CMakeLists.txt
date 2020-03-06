# CMake 3.3 or later is required for COMPILE_LANGUAGE.
# It can be easily obtained and used without superuser privileges
# even if the system has an older version of CMake installed.
cmake_minimum_required(VERSION 3.3)

if (CMAKE_MAJOR_VERSION LESS 3)
  # CMake 3.0 is used for generator expressions in install(FILES).
  message(WARNING
    "CMake version 3.0 or newer is required to install all targets")
endif ()

# Set the path to CMake modules.
set(MP_CMAKE_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MP_CMAKE_MODULE_DIR})

include(init)

project(MP)

# Set the c++11 flag without GNU extensions.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (POLICY CMP0042)
  # Enable @rpath support on OS X.
  cmake_policy(SET CMP0042 NEW)
endif ()

if (POLICY CMP0054)
  # Only interpret `if` arguments as variables or keywords when unquoted.
  cmake_policy(SET CMP0054 NEW)
endif ()

set(MP_VERSION 3.1.0)
if (NOT MP_VERSION MATCHES "^([0-9]+).([0-9]+).([0-9]+)$")
  message(FATAL_ERROR "Invalid version format ${MP_VERSION}.")
endif ()
set(MP_VERSION_MAJOR ${CMAKE_MATCH_1})
set(MP_VERSION_MINOR ${CMAKE_MATCH_2})
set(MP_VERSION_PATCH ${CMAKE_MATCH_3})

message(STATUS "System: ${CMAKE_SYSTEM}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# Set output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set the default Java search path.
if (UNIX)
  set(JAVA_HOME /usr/lib/jvm/java-7-oracle)
endif ()

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
include(CheckCCompilerFlag)
function (add_c_compiler_flags target)
  foreach (flag ${ARGN})
    string(REPLACE "-" "_" var ${flag})
    check_c_compiler_flag(${flag} ${var})
    if (${var})
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

if (NOT MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
else ()
  # Disable useless MSVC warnings suggesting nonportable "secure" alternatives.
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(MP_CLANG TRUE)
endif ()

# Adds a prefix to arguments.
function (add_prefix var prefix)
  set(result ${${var}})
  foreach (arg ${ARGN})
    set(result ${result} "${prefix}${arg}")
  endforeach ()
  set(${var} ${result} PARENT_SCOPE)
endfunction ()

set_cache(AMPL_LIBRARY_DIR bin STRING
  "A directory to install AMPL function libraries relative to "
  "${CMAKE_INSTALL_PREFIX}")

# Adds a shared AMPL library which by convention doesn't have any prefix
# and has a suffix ".dll" on all platforms.
macro(add_ampl_library name)
  cmake_parse_arguments(add_ampl_library PRIVATE "" "" ${ARGN})
  add_library(${name} SHARED ${add_ampl_library_UNPARSED_ARGUMENTS})
  set_target_properties(${name} PROPERTIES PREFIX "")
  set_target_properties(${name} PROPERTIES SUFFIX ".dll")
  target_link_libraries(${name} asl)
  if (NOT add_ampl_library_PRIVATE)
    # Specify RUNTIME DESTINATION and LIBRARY DESTINATION, but not
    # DESTINATION or ARCHIVE_DESTINATION because we don't want import
    # libraries installed.
    install(TARGETS ${name} RUNTIME
      DESTINATION ${AMPL_LIBRARY_DIR} LIBRARY DESTINATION ${AMPL_LIBRARY_DIR})
  endif ()
endmacro()

include_directories(include)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

set(MP_DATE 20141202)
set(MP_SYSINFO "${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")

function (add_mp_library name)
  cmake_parse_arguments(add_mp_library "OBJECT;STATIC" ""
    "COMPILE_DEFINITIONS;INCLUDE_DIRECTORIES;LIBRARIES;OBJECT_LIBRARIES;DEPENDS"
    ${ARGN})
  # Get object files.
  set(objects)
  set(dynrt-objects)
  foreach (lib ${add_mp_library_OBJECT_LIBRARIES})
    set(objects ${objects} $<TARGET_OBJECTS:${lib}>)
    set(dynrt-objects ${dynrt-objects} $<TARGET_OBJECTS:${lib}-dynrt>)
  endforeach ()
  # Add library.
  set(libtype )
  if (add_mp_library_OBJECT)
    set(libtype OBJECT)
  elseif (add_mp_library_STATIC)
    set(libtype STATIC)
  endif ()
  add_library(${name} ${libtype}
    ${add_mp_library_UNPARSED_ARGUMENTS} ${objects})
  target_compile_definitions(${name}
    PUBLIC ${add_mp_library_COMPILE_DEFINITIONS})
  target_include_directories(${name}
    PUBLIC ${add_mp_library_INCLUDE_DIRECTORIES})
  if (add_mp_library_DEPENDS)
    add_dependencies(${name} ${add_mp_library_DEPENDS})
  endif ()
  # Add library linked with dynamic runtime.
  if (MSVC)
    add_library(${name}-dynrt ${libtype} EXCLUDE_FROM_ALL
      ${add_mp_library_UNPARSED_ARGUMENTS} ${dynrt-objects})
    target_compile_options(${name}-dynrt PUBLIC /MD$<$<CONFIG:Debug>:d>)
    target_compile_definitions(${name}-dynrt
      PUBLIC ${add_mp_library_COMPILE_DEFINITIONS})
    target_include_directories(${name}-dynrt
      PUBLIC ${add_mp_library_INCLUDE_DIRECTORIES})
    if (add_mp_library_DEPENDS)
      add_dependencies(${name}-dynrt ${add_mp_library_DEPENDS})
    endif ()
  endif ()
  # Link libraries.
  foreach (lib ${add_mp_library_LIBRARIES})
    target_link_libraries(${name} ${lib})
    if (MSVC)
      target_link_libraries(${name}-dynrt ${lib}-dynrt)
    endif ()
  endforeach ()
endfunction ()

# Link target with libraries built with dynamic runtime.
function (target_link_libraries_dynrt target)
  foreach (lib ${ARGN})
    if (TARGET ${lib}-dynrt)
      target_link_libraries(${target} ${lib}-dynrt)
    else ()
      target_link_libraries(${target} ${lib})
    endif ()
  endforeach ()
endfunction ()

add_mp_library(format OBJECT
  include/mp/format.h include/mp/posix.h src/format.cc src/posix.cc)

set(MP_EXPR_INFO_FILE ${MP_SOURCE_DIR}/src/expr-info.cc)
add_executable(gen-expr-info EXCLUDE_FROM_ALL
  src/gen-expr-info.cc $<TARGET_OBJECTS:format>)
if (MINGW)
  SET_TARGET_PROPERTIES(gen-expr-info PROPERTIES
    LINK_FLAGS "-static-libgcc -static-libstdc++")
endif ()
# use wine if cross compiling from unix to windows
if (CMAKE_HOST_UNIX AND WIN32)
  set(WINE wine)
else ()
  set(WINE "")
endif ()
if (CMAKE_CROSSCOMPILING)
  # Produce a warning because expr-info.cc can be out of date but cannot be
  # re-generated because we are cross compiling.
  add_custom_command(OUTPUT ${MP_EXPR_INFO_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo
      "warning: cannot re-generate ${MP_EXPR_INFO_FILE}")
else ()
  add_custom_command(OUTPUT ${MP_EXPR_INFO_FILE}
    COMMAND ${WINE} $<TARGET_FILE:gen-expr-info> ${MP_EXPR_INFO_FILE}
    DEPENDS gen-expr-info)
endif ()

add_prefix(MP_HEADERS include/mp/
  arrayref.h basic-expr-visitor.h clock.h common.h error.h expr.h
  expr-visitor.h nl.h nl-reader.h option.h os.h problem.h problem-builder.h
  rstparser.h safeint.h sol.h solver.h suffix.h)
set(MP_SOURCES )
add_prefix(MP_SOURCES src/
  clock.cc expr.cc expr-writer.h nl-reader.cc option.cc os.cc
  problem.cc rstparser.cc sol.cc solver.cc solver-c.h sp.h sp.cc)

add_mp_library(mp ${MP_HEADERS} ${MP_SOURCES} ${MP_EXPR_INFO_FILE}
  COMPILE_DEFINITIONS MP_DATE=${MP_DATE} MP_SYSINFO="${MP_SYSINFO}"
  INCLUDE_DIRECTORIES src OBJECT_LIBRARIES format)
set_target_properties(mp PROPERTIES
  VERSION ${MP_VERSION} SOVERSION ${MP_VERSION_MAJOR})

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS ${MP_CXX11_FLAG})

check_cxx_source_compiles(
  "#include <memory>
  int main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
if (HAVE_UNIQUE_PTR)
  target_compile_definitions(mp PUBLIC MP_USE_UNIQUE_PTR)
endif ()

check_cxx_source_compiles(
  "#include <atomic>
  int main() { std::atomic<int> n; }" HAVE_ATOMIC)
if (HAVE_ATOMIC)
  target_compile_definitions(mp PUBLIC MP_USE_ATOMIC)
endif ()

check_cxx_source_compiles(
  "#include <functional>
  int main() { std::hash<int> h; }" HAVE_HASH)
if (HAVE_HASH)
  target_compile_definitions(mp PUBLIC MP_USE_HASH)
endif ()

set(CMAKE_REQUIRED_FLAGS )

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(mp ${RT_LIBRARY})
endif ()

# Check if variadic templates are working and not affected by GCC bug 39653:
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=39653
check_cxx_source_compiles("
  template <class T, class ...Types>
  struct S { typedef typename S<Types...>::type type; };
  int main() {}" MP_VARIADIC_TEMPLATES)

if (MP_VARIADIC_TEMPLATES)
  add_executable(nl-example src/nl-example.cc)
  target_link_libraries(nl-example mp)
endif ()

add_subdirectory(doc)
add_subdirectory(src/amplsig)
add_subdirectory(src/asl)
add_subdirectory(src/cp)
add_subdirectory(solvers)

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

include(CTest)
enable_testing()
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

install(DIRECTORY include/mp DESTINATION include)
install(TARGETS mp DESTINATION lib RUNTIME DESTINATION bin)
install(FILES LICENSE.rst DESTINATION share/mp)

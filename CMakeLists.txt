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

## Use CCache to optimize recompilation
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
  set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
  set(CMAKE_Fortran_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
else()
  message(STATUS "   HINT: Install ccache for faster recompilation")
endif()

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

option(USE_PROFILER "Use Profiler" OFF)
if(USE_PROFILER)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

# Sanitizers
option(USE_SANITIZERS "Use Sanitizers" OFF)
if(USE_SANITIZERS)
  if (WIN32)
    message(STATUS "   MP: Using address sanitizers on Windows")
    set(CMAKE_C_FLAGS
      "${CMAKE_CXX_FLAGS} /fsanitize=address")
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} /fsanitize=address")
  else()
    message(STATUS "   MP: Using sanitizers on *nix")
    set(CMAKE_C_FLAGS
      "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
  endif()
endif()

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
  set(MP_FOLDER_PREFIX "mp/")
else()
  set(MP_FOLDER_PREFIX "")
endif()
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

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

set_property(GLOBAL PROPERTY MP_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
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

# Add files specified as additional arguments to the specified folder.
function (add_to_folder folder)
  foreach (target ${ARGN})
	if(TARGET ${target})
		set_property(TARGET ${target} PROPERTY FOLDER ${folder})
	endif()
  endforeach ()
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
  add_to_folder(${MP_FOLDER_PREFIX}ampl-libraries ${name})
endmacro()

include_directories(include)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(MP_DATE 20231017)

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
  set_property(TARGET ${name}  PROPERTY CXX_STANDARD 11)
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
  if(hasParent)
    add_to_folder(mp ${name} ${name}-dynrt)
  endif()
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

# gen-expr-info
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
add_to_folder(${MP_FOLDER_PREFIX}util gen-expr-info)


# HEADER LISTS
add_prefix(MP_HEADERS include/mp/
  ampls-c-api.h ampls-cpp-api.h arrayref.h basic-expr-visitor.h
  backend-app.h backend-base.h backend-mip.h backend-std.h
  backend-to-model-api.h backend-with-mm.h backend-with-valcvt.h
  clock.h common.h converter-base.h
  easy-modeler.h
  env.h error.h expr.h expr-visitor.h
  model-mgr-base.h model-mgr-with-pb.h
  model-mgr-with-std-pb.h model-mgr-with-std-pb.hpp
  nl.h nl-reader.h option.h os.h problem.h problem-builder.h
  valcvt-base.h valcvt-node.h valcvt-link.h valcvt.h
  rstparser.h safeint.h sol.h
  solver.h solver-opt.h solver-base.h solver-io.h solver-app-base.h solver-app.h
  suffix.h
  utils-file.h utils-hash.h utils-hash-stream.h
  utils-string.h utils-math.h)

add_prefix(MP_FLAT_HEADERS include/mp/flat/
  backend_flat.h
  constr_base.h constr_eval.h constr_keeper.h constr_hash.h
  constr_std.h constr_static.h constr_functional.h
  constr_algebraic.h constr_general.h
  constr_prepro.h constr_prop_down.h
  context.h
  converter.h converter_model.h converter_model_base.h
  convert_functional.h
  eexpr.h expr_algebraic.h expr_affine.h expr_quadratic.h
  expr_bounds.h
  model_api_base.h model_api_connect.h
  model_info.h model_info.hpp
  obj_std.h preprocess.h problem_flattener.h)

add_prefix(MP_FLAT_REDEF_HEADERS include/mp/flat/redef/
  encodings.h redef_base.h)

add_prefix(MP_FLAT_REDEF_STD_HEADERS include/mp/flat/redef/std/
  range_con.h)

add_prefix(MP_FLAT_REDEF_STD_HEADERS include/mp/flat/redef/conic/
  cones.h qcones2qc.h)

add_prefix(MP_FLAT_REDEF_MIP_HEADERS include/mp/flat/redef/MIP/
  converter_mip.h redefs_mip_std.h abs.h alldiff.h
  cond_eq.h cond_ineq.h complement.h count.h
  div.h
  ifthenelse.h impl.h
  indicator_eq.h indicator_le.h indicator_ge.h
  indicator_quad.h
  lin_approx.h core/lin_approx_core.h
  logical_and.h logical_or.h logical_not.h
  min_max.h mul.h numberof_const.h numberof_var.h
  piecewise_linear.h
  power_const.h sos2.h)


set(MP_ALL_HEADERS ${MP_HEADERS} ${MP_FLAT_HEADERS}
  ${MP_FLAT_REDEF_HEADERS}
  ${MP_FLAT_REDEF_STD_HEADERS} ${MP_FLAT_REDEF_MIP_HEADERS})


set(MP_SOURCES )
add_prefix(MP_SOURCES src/
  clock.cc expr.cc expr-writer.h nl-reader.cc option.cc os.cc
  problem.cc rstparser.cc sol.cc solver.cc solver-c.h sp.h sp.cc
  std_constr.cc
  utils_file.cc utils_string.cc)

set(MP_FLAT_SOURCES )
add_prefix(MP_FLAT_SOURCES src/mp/flat/
  encodings.cpp piecewise_linear.cpp)

add_mp_library(mp ${MP_ALL_HEADERS}
    ${MP_SOURCES} ${MP_FLAT_SOURCES} ${MP_EXPR_INFO_FILE}
  COMPILE_DEFINITIONS MP_DATE=${MP_DATE} MP_SYSINFO="${MP_SYSINFO}"
  INCLUDE_DIRECTORIES src include OBJECT_LIBRARIES format)
set_target_properties(mp PROPERTIES
  VERSION ${MP_VERSION} SOVERSION ${MP_VERSION_MAJOR})

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS ${MP_CXX11_FLAG})

check_cxx_source_compiles(
  "#include <memory>
  int main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
if (HAVE_UNIQUE_PTR)
  target_compile_definitions(mp PUBLIC MP_USE_UNIQUE_PTR)
  if(TARGET mp-dynrt)
    target_compile_definitions(mp-dynrt PUBLIC MP_USE_UNIQUE_PTR)
  endif()
endif ()

check_cxx_source_compiles(
  "#include <atomic>
  int main() { std::atomic<int> n; }" HAVE_ATOMIC)
if (HAVE_ATOMIC)
  target_compile_definitions(mp PUBLIC MP_USE_ATOMIC)
  if(TARGET mp-dynrt)
    target_compile_definitions(mp-dynrt PUBLIC MP_USE_ATOMIC)
  endif()
endif ()

check_cxx_source_compiles(
  "#include <functional>
  int main() { std::hash<int> h; }" HAVE_HASH)
if (HAVE_HASH)
  target_compile_definitions(mp PUBLIC MP_USE_HASH)
   if(TARGET mp-dynrt)
    target_compile_definitions(mp-dynrt PUBLIC MP_USE_HASH)
  endif()
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

# Option BUILD_EXAMPLES
option(BUILD_EXAMPLES "Build examples: nl reader, nl writer" ON)

# Build NLReader example
if (MP_VARIADIC_TEMPLATES AND BUILD_EXAMPLES)
  message(STATUS
    "Building nl-reader-example. Set -DBUILD_EXAMPLES=off to switch off.")
  add_executable(nl-reader-example examples/nl-reader-example.cc)
  add_to_folder(${MP_FOLDER_PREFIX}util nl-reader-example)
  target_link_libraries(nl-reader-example mp)
endif ()

# Build NLWriter2 example
if (BUILD_EXAMPLES)
  message(STATUS
    "Building NLWriter2 example. Set -DBUILD_EXAMPLES=off to switch off.")
  set(NLW2_DIR nl-writer2)
  add_subdirectory(${NLW2_DIR})
  ## Setting source files for docu. TODO keep up2date
  set(NLW2_LIB_FILES
          ${NLW2_DIR}/src/nl-writer2.cc
          ${NLW2_DIR}/src/nl-utils2.cc
          ${NLW2_DIR}/src/dtoa.c)
  set(NLW2_INC_FILES
                ${NLW2_DIR}/include/mp/nlsol.h
                ${NLW2_DIR}/include/mp/nl-header.h
                ${NLW2_DIR}/include/mp/nl-feeder2.h
                ${NLW2_DIR}/include/mp/nl-writer2.h
                ${NLW2_DIR}/include/mp/nl-writer2.hpp
                ${NLW2_DIR}/include/mp/nl-writer2-misc.h
                ${NLW2_DIR}/include/mp/nl-utils2.h
                ${NLW2_DIR}/include/mp/sol-reader2.h
                ${NLW2_DIR}/include/mp/sol-reader2.hpp
                ${NLW2_DIR}/include/mp/sol-handler2.h
                )
  set(NLSOL_EX_FILES
    ${NLW2_DIR}/examples/nlsol_ex.cc
    ${NLW2_DIR}/examples/nlsol_ex_mdl.h
    ${NLW2_DIR}/examples/nlsol_ex_nl.h
    ${NLW2_DIR}/examples/nlsol_ex_sol.h)
  # Export source files
  set(NLW2_DOC_SRC
    ${NLW2_LIB_FILES} ${NLW2_INC_FILES} ${NLSOL_EX_FILES})
endif()

option(BUILD_DOC "Build documentation" ON)
if(BUILD_DOC)
  add_subdirectory(doc)
endif()
add_subdirectory(src/asl)
add_subdirectory(src/amplsig)
add_subdirectory(src/cp)
add_subdirectory(solvers)

# Add a target that generates a file with solver and library versions.
set(AMPL_VERSIONS_FILE versions.txt)
add_custom_target(versions cmake -E remove ${AMPL_VERSIONS_FILE})
add_to_folder(${MP_FOLDER_PREFIX}util versions)
foreach (target gsl-info gecode jacop)
  if (TARGET ${target})
    add_custom_command(TARGET versions POST_BUILD
        COMMAND $<TARGET_FILE:${target}> -v >> ${AMPL_VERSIONS_FILE})
  endif ()
endforeach ()

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

# Option BUILD_TESTS
option(BUILD_TESTS "Build unit tests" ON)
include(CTest)
enable_testing()
if(BUILD_TESTS)
  message(STATUS "Building unit tests. Set -DBUILD_TESTS=off to switch off.")
  add_subdirectory(test)
endif()

install(DIRECTORY include/mp DESTINATION include)
install(TARGETS mp DESTINATION lib RUNTIME DESTINATION bin)
install(FILES LICENSE.rst DESTINATION share/mp)

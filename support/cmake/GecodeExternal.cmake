# Gecode, built from the sources in <mp_dir>/thirdparty/gecode (submodule)
# as an external project and installed into <binary_dir>/thirdparty/gecode:
# its CMake project defines a target named "gecode", which would clash with
# the driver target under add_subdirectory().
#
# add_gecode_external(<mp_dir> <binary_dir>) defines
#   gecode-libs     the external project (configure, build, install);
#   mp-gecode-libs  INTERFACE library with the installed headers, the static
#                   libraries and Threads, depending on gecode-libs.
# Static only, PIC, Release; FlatZinc, examples, Gist, MPFR, CPProfiler and
# tests off.
function(add_gecode_external mp_dir binary_dir)
  if (TARGET mp-gecode-libs)
    return()
  endif ()
  find_package(Git)
  if (GIT_FOUND AND NOT EXISTS ${mp_dir}/thirdparty/gecode/CMakeLists.txt)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init thirdparty/gecode
      WORKING_DIRECTORY ${mp_dir})
  endif ()
  include(ExternalProject)
  set(GECODE_PREFIX ${binary_dir}/thirdparty/gecode)
  set(GECODE_LIBS )
  foreach (lib minimodel driver float set int search kernel support)
    list(APPEND GECODE_LIBS ${GECODE_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gecode${lib}${CMAKE_STATIC_LIBRARY_SUFFIX})
  endforeach ()
  set(GECODE_CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DGECODE_BUILD_SHARED=OFF -DGECODE_BUILD_STATIC=ON
    -DGECODE_ENABLE_FLATZINC=OFF -DGECODE_ENABLE_EXAMPLES=OFF
    -DGECODE_ENABLE_GIST=OFF -DGECODE_ENABLE_MPFR=OFF
    -DGECODE_ENABLE_CPPROFILER=OFF -DBUILD_TESTING=OFF)
  get_property(multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if (NOT multi_config)
    list(APPEND GECODE_CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release)
  endif ()
  foreach (var CMAKE_C_COMPILER CMAKE_CXX_COMPILER
      CMAKE_C_COMPILER_LAUNCHER CMAKE_CXX_COMPILER_LAUNCHER
      CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET
      CMAKE_MSVC_RUNTIME_LIBRARY)
    if (${var})
      string(REPLACE ";" "|" value "${${var}}")
      list(APPEND GECODE_CMAKE_ARGS "-D${var}=${value}")
    endif ()
  endforeach ()
  find_package(Threads REQUIRED)
  ExternalProject_Add(gecode-libs
    SOURCE_DIR ${mp_dir}/thirdparty/gecode
    BINARY_DIR ${binary_dir}/thirdparty/build/gecode
    INSTALL_DIR ${GECODE_PREFIX}
    LIST_SEPARATOR |
    CMAKE_ARGS ${GECODE_CMAKE_ARGS}
    BUILD_BYPRODUCTS ${GECODE_LIBS})
  file(MAKE_DIRECTORY ${GECODE_PREFIX}/include)
  add_library(mp-gecode-libs INTERFACE)
  target_include_directories(mp-gecode-libs INTERFACE ${GECODE_PREFIX}/include)
  target_link_libraries(mp-gecode-libs INTERFACE ${GECODE_LIBS} Threads::Threads)
  add_dependencies(mp-gecode-libs gecode-libs)
endfunction()

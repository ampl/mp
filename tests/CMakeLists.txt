# CMake build script for the AMPL solver tests.

# GTest doesn't detect <tuple> with clang.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
endif ()

# Enable std::this_thread::yield.
add_definitions(-D_GLIBCXX_USE_SCHED_YIELD=1)

include_directories(.)

# Add gtest library with some extra features.
find_package(Threads)
add_library(gtest gtest/gtest-all.cc gtest/gtest.h util.h util.cc test-main.cc)
target_link_libraries(gtest util ${CMAKE_THREAD_LIBS_INIT})

# Adds an AMPL test.
# Usage:
#   add_ampl_test(name [PREFIX] sources... [LIBS libraries])
function(add_ampl_test name)
  set(srcs )
  set(libs )
  set(library FALSE)
  foreach (arg ${ARGN})
    if (library)
      set (libs ${libs} ${arg})
    elseif (arg STREQUAL LIBS)
      set (library TRUE)
    elseif (arg STREQUAL PREFIX)
      set(prefix ampl)
    else ()
      set (srcs ${srcs} ${arg})
    endif ()
  endforeach ()
  set(target_name "${prefix}${name}")
  add_executable(${target_name} ${srcs})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} gtest util ${libs})
  add_test(${name} ${name})
endfunction()

# Check C++ compiler capabilities.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "int main() { int *p = nullptr; }" HAVE_NULLPTR)
check_cxx_source_compiles(
  "int main() { [](int x) { return x; }; }" HAVE_LAMBDAS)
check_cxx_source_compiles(
  "#include <thread>
  void run() {}
  int main() { std::thread t(run); }" HAVE_THREADS)

find_program(HAVE_LSOF lsof)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

add_ampl_library(testlib testlib.cc)
add_library(function function.cc function.h)
add_dependencies(function arith_h)
target_link_libraries(function amplsolver)

add_ampl_test(function-test function-test.cc LIBS function)

if (HAVE_CBC)
  add_ampl_test(cbc-test cbc-test.cc)
endif ()

if (HAVE_GSL)
  include_directories(${GSL_BINARY_DIR})
  add_ampl_test(gsl-test gsl-test.cc LIBS function gsl gslcblas)
endif ()

add_subdirectory(ilogcp)
add_subdirectory(solvers)
add_subdirectory(tables)
add_subdirectory(util)

add_ampl_test(util-test util-test.cc)

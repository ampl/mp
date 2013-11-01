# CMake build script for the AMPL solver tests.

# GTest doesn't detect <tuple> with clang.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
endif ()

# Setup testing.
add_subdirectory(gtest)
include_directories(gtest/include)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/gtest)

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
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

add_ampl_library(testlib testlib.cc)
add_library(function function.cc function.h)
add_dependencies(function arith_h)
target_link_libraries(function amplsolver)

cxx_test_with_flags(function_test "" "gtest;function"
  test_main.cc function_test.cc)

if (HAVE_CBC)
  cxx_test_with_flags(cbc_test "" "gtest;util" test_main.cc cbc_test.cc util.cc)
endif ()

if (HAVE_GSL)
  include_directories(${GSL_BINARY_DIR})
  cxx_test_with_flags(gsl_test ""
    "gtest;function;gsl;gslcblas" test_main.cc gsl_test.cc)
endif ()

add_subdirectory(ilogcp)
add_subdirectory(solvers)
add_subdirectory(tables)
add_subdirectory(util)

cxx_test_with_flags(util_test "" "gtest;util" test_main.cc util.cc util_test.cc)

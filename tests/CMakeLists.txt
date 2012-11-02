# CMake build script for the AMPL solver tests.

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAS_WNO_LONG_LONG_FLAG)
if (HAS_WNO_LONG_LONG_FLAG)
  # Required if -pedantic is used. The long long type is used by GTest.
  add_definitions(-Wno-long-long)
endif ()

# Setup testing.
add_subdirectory(gtest)
include_directories(gtest/include)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/gtest)

# Check C++ compiler capabilities.
include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
  "int main() { int *p = nullptr; }" HAVE_NULLPTR)
CHECK_CXX_SOURCE_COMPILES(
  "#include <memory>\nint main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
CHECK_CXX_SOURCE_COMPILES(
  "int main() { [](int x) { return x; }; }" HAVE_LAMBDAS)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

cxx_test_with_flags(expr_test "" "gtest;util"
  test_main.cc expr_builder.cc expr_test.cc)

add_ampl_library(testlib testlib.cc)
add_library(function function.cc function.h)
add_dependencies(function arith_h)
target_link_libraries(function ${SOLVER_LIBS})

cxx_test_with_flags(function_test "" "gtest;function"
  test_main.cc function_test.cc)

if (HAVE_GSL)
  include_directories(${GSL_BINARY_DIR})
  cxx_test_with_flags(gsl_test ""
    "gtest;function;gsl;gslcblas" test_main.cc gsl_test.cc)
endif ()

find_package(ODBC)
if (ODBC_FOUND)
  cxx_test_with_flags(odbc_test ""
    "gtest;function;${ODBC_LIBRARIES}" test_main.cc odbc_test.cc)
endif ()

add_subdirectory(ilogcp)

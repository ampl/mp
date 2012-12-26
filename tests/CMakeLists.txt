# CMake build script for the AMPL solver tests.

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAVE_WNO_LONG_LONG_FLAG)
if (HAVE_WNO_LONG_LONG_FLAG)
  # Required if -pedantic is used. The long long type is used by GTest.
  add_definitions(-Wno-long-long)
endif ()

check_cxx_compiler_flag(-Wno-variadic-macros HAVE_WNO_VARIADIC_MACROS_FLAG)
if (HAVE_WNO_VARIADIC_MACROS_FLAG)
  # Variadic macros are used by GTest.
  add_definitions(-Wno-variadic-macros)
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
  "#include <memory>
  int main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

cxx_test_with_flags(expr_test "" "gtest;util"
  test_main.cc expr_builder.cc expr_test.cc no_unordered_map_test.cc)

cxx_test_with_flags(format_test "" "gtest;util" format_test.cc)
cxx_test_with_flags(driver_test "" "gtest;util"
  driver_test.cc args.cc test_main.cc)

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
  include(CheckCXXSourceRuns)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/..)
  set(CMAKE_REQUIRED_LIBRARIES ${ODBC_LIBRARIES})
  check_cxx_source_runs(
    "#include \"tests/odbc.cc\"
     int main(int argc, char **argv) {
       return odbc::Env().FindDriver(\"mysql\").empty() ? 1 : 0;
     }" HAVE_MYSQL)

  cxx_test_with_flags(odbc_test ""
    "gtest;function;${ODBC_LIBRARIES}" test_main.cc odbc_test.cc odbc.cc)

  if (HAVE_MYSQL)
    if (WIN32)
      set(MYSQL_TEST_LIBRARIES Ws2_32.lib)
    endif ()
    cxx_test_with_flags(mysql_test ""
      "gtest;function;${ODBC_LIBRARIES};${MYSQL_TEST_LIBRARIES}"
      test_main.cc mysql_test.cc odbc.cc)
  endif ()
endif ()

find_package(Gecode)
if (GECODE_FOUND)
  cxx_test_with_flags(gecode_test "" "gtest;gecodedriver"
    test_main.cc args.cc expr_builder.cc gecode_test.cc)
endif ()

add_subdirectory(ilogcp)

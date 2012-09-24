# CMake build script for the AMPL solver tests.

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif ()

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

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS})
  include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIRS} ${CPLEX_CP_INCLUDE_DIRS})

  cxx_test_with_flags(ilogcp_test ""
    "amplilogcp;gtest;gtest_main;${SOLVER_LIBS}" ilogcp_test.cc)

  target_link_libraries(ilogcp_test
    ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})
endif ()

if (HAVE_GSL)
  include_directories(../solvers/amplgsl/gsl/include)
  cxx_test_with_flags(gsl_test ""
    "gtest;gtest_main;${SOLVER_LIBS}"
    gsl_test.cc function.cc function.h function_test.cc)
  target_link_libraries(gsl_test gsl gslcblas)
endif ()

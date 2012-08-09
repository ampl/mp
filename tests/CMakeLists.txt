# CMake build script for the AMPL solver tests.

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif ()

# Setup testing.
add_subdirectory(gtest)
include_directories(gtest/include)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/gtest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Check C++ compiler capabilities.
include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
  "int main() { int *p = nullptr; }" HAVE_NULLPTR)
CHECK_CXX_SOURCE_COMPILES(
  "#include <memory>\nint main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS})
  include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIRS} ${CPLEX_CP_INCLUDE_DIRS})

  cxx_test_with_flags(concert_test ""
    "amplconcert;gtest;gtest_main;${SOLVER_LIBS}" concert_test.cc)

  target_link_libraries(concert_test
    ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})
endif ()

find_package(GSL)
if (GSL_FOUND)
  cxx_test_with_flags(gsl_test ""
    "gtest;gtest_main;${SOLVER_LIBS};${GSL_LIBRARIES}"
    gsl_test.cc functional.h)
endif ()
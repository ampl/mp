# CMake build script for the AMPL solver tests.

cmake_minimum_required(VERSION 2.6)
project(TESTS)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Setup testing.
enable_testing()
add_subdirectory(gtest)
include_directories(gtest/include)
link_directories(${TESTS_BINARY_DIR}gtest)

# Find the AMPL solver library.
include_directories(..)
if (NOT WIN32)
  find_path(SOLVER_LIB_DIR amplsolver.a ../solvers)
  set(SOLVER_LIBS "${SOLVER_LIB_DIR}/amplsolver.a" dl)
else ()
  set(SOLVER_LIBS ../solvers/amplsolv)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
if (NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -std=c++0x")
endif ()

include(CheckCXXSourceCompiles)
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
    "gtest;gtest_main;${SOLVER_LIBS}"
    ../solvers/concert/build_constr.cpp
    ../solvers/concert/build_expr.cpp
    ../solvers/concert/build_numberof.cpp
    ../solvers/concert/concert.cpp
    ../solvers/concert/concert.h
    ../solvers/concert/objconst0.c
    concert_test.cc)

  target_link_libraries(concert_test
    ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})
endif ()


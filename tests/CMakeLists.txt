# CMake build script for the AMPL solver tests.

cmake_minimum_required(VERSION 2.6)
project(TESTS)

# Setup testing.
enable_testing()
add_subdirectory(gtest)
include_directories(gtest/include)
link_directories(${TESTS_BINARY_DIR}gtest)

include_directories(../solvers)
find_path(SOLVER_LIB_DIR amplsolver.a ../solvers)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# TODO: move to FindCPLEX.cmake
set(SYSTEM "x86-64_sles10_4.1")
set(LIBFORMAT static_pic)

set(CPLEX_STUDIO_DIR /opt/ibm/ILOG/CPLEX_Studio124)
set(CPLEX_DIR ${CPLEX_STUDIO_DIR}/cplex)
set(CONCERT_DIR ${CPLEX_STUDIO_DIR}/concert)
set(CPOPTIMIZER_DIR ${CPLEX_STUDIO_DIR}/cpoptimizer)

add_definitions(-DIL_STD)

include_directories(
  ${CPLEX_DIR}/include
  ${CONCERT_DIR}/include
  ${CPOPTIMIZER_DIR}/include)

link_directories(
  ${CPLEX_DIR}/lib/${SYSTEM}/${LIBFORMAT}
  ${CONCERT_DIR}/lib/${SYSTEM}/${LIBFORMAT}
  ${CPOPTIMIZER_DIR}/lib/${SYSTEM}/${LIBFORMAT})

cxx_test_with_flags(cplexcp_test ""
  "gtest;gtest_main;cp;concert;ilocplex;cplex;${SOLVER_LIB_DIR}/amplsolver.a;dl"
  ../solvers/cplexcp/build_constr.cpp
  ../solvers/cplexcp/build_expr.cpp
  ../solvers/cplexcp/build_numberof.cpp
  ../solvers/cplexcp/concert.cpp
  ../solvers/cplexcp/objconst0.c
  cplexcp_test.cc)

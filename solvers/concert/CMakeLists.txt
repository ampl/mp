# CMake build script for the IBM ILOG Concert/CP driver for AMPL.

cmake_minimum_required(VERSION 2.6)
project(CONCERT)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake")

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE Release)
endif ()

# Find the AMPL solver library.
include_directories(../..)
if (NOT WIN32)
  find_path(SOLVER_LIB_DIR amplsolver.a ..)
  set(SOLVER_LIBS "${SOLVER_LIB_DIR}/amplsolver.a" dl)
elseif (NOT SOLVER_LIBS)
  set(SOLVER_LIBS ../amplsolv)
endif ()

if (NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
endif ()

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS})
  include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIRS} ${CPLEX_CP_INCLUDE_DIRS})

  add_executable(concert
    build_expr.cpp
    concert.cpp
    concert.h
    main.cpp
    util.cpp
    util.h)

  target_link_libraries(concert ${SOLVER_LIBS}
    ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})
endif ()



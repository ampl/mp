# CMake build script for the IBM ILOG Concert/CP driver for AMPL.

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE Release)
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

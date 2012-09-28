# CMake build script for the IBM ILOG Concert/CP driver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})
  include_directories(.. ${CPLEX_ILOCPLEX_INCLUDE_DIRS} ${CPLEX_CP_INCLUDE_DIRS})

  add_library(amplilogcp
    build_expr.cpp
    ilogcp.cpp
    ilogcp.h
    ilogcp_date.h
    util.cpp
    util.h)

  add_executable(ilogcp main.cpp)

  target_link_libraries(ilogcp amplilogcp ${SOLVER_LIBS}
    debug ${CPLEX_ILOCPLEX_LIBRARIES_DEBUG} ${CPLEX_CP_LIBRARIES_DEBUG}
    ${CPLEX_LIBRARIES_DEBUG}
    optimized ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES}
    ${CPLEX_LIBRARIES})
endif ()

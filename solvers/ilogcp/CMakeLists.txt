# CMake build script for the IBM ILOG Concert/CP driver for AMPL.

add_definitions(${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})
include_directories(.. ${CPLEX_ILOCPLEX_INCLUDE_DIRS} ${CPLEX_CP_INCLUDE_DIRS})

if (MSVC)
  # CP Optimizer only provides libraries built with DLL MSVC runtime
  # so use the same runtime for ilogcp.
  replace_flags("/MT" "/MD")
endif ()

add_library(amplilogcp
  build_expr.cpp
  ilogcp.cpp
  ilogcp.h
  ilogcp_date.h
  util.cpp
  util.h)
add_dependencies(amplilogcp arith_h)

add_executable(ilogcp main.cpp)

target_link_libraries(ilogcp amplilogcp ${SOLVER_LIBS}
  ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})

# CMake build script for the IBM ILOG CPLEX solver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND)
  include_directories(.. ${CPLEX_INCLUDE_DIRS})
  add_executable(cplex cplex.c ../funcadd0.c)
  target_link_libraries(cplex ${CPLEX_LIBRARIES} ${SOLVER_LIBS})
  add_dependencies(cplex arith_h)
endif ()

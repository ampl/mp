# CMake build script for the IBM ILOG CPLEX solver for AMPL.

# Adds a cplex solver.
function(add_cplex_solver name)
  add_executable(${name} cplex.c
    ${PROJECT_SOURCE_DIR}/src/asl/solvers/funcadd0.c)
  target_include_directories(${name} PRIVATE ${CPLEX_INCLUDE_DIR}/ilcplex)
  target_link_libraries(${name} cplex-library asl)

  # Suppress warnings.
  if (CMAKE_COMPILER_IS_GNUCC)
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-Wno-unused-result HAVE_WNO_UNUSED_RESULT)
    if (HAVE_WNO_UNUSED_RESULT)
      target_compile_options(${name} PRIVATE -Wno-unused-result)
    endif ()
    target_compile_options(${name} PRIVATE
      -Wno-deprecated-declarations -Wno-long-long -Wno-sign-compare
      -Wno-uninitialized -Wno-unused-parameter -Wno-missing-field-initializers)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(${name} PRIVATE
      -Wno-self-assign -Wno-sign-compare -Wno-sometimes-uninitialized
      -Wno-string-plus-int -Wno-unused-parameter -Wno-missing-field-initializers)
  elseif (MSVC)
    target_compile_options(${name} PRIVATE /wd4244 /wd4267 /wd4996)
  endif ()
endfunction ()

find_package(CPLEX)
if (CPLEX_FOUND)
  add_cplex_solver(cplex)

  # Version of CPLEX with SOCP transformations.
  add_library(socp STATIC socp.cc)
  target_link_libraries(socp asl)
  add_cplex_solver(cplex-socp)
  target_compile_definitions(cplex-socp PRIVATE
    jac0dim_ASL=socp_jac0dim write_sol_ASL=socp_write_sol
    qp_read_ASL=socp_qp_read)
  target_link_libraries(cplex-socp socp)
endif ()

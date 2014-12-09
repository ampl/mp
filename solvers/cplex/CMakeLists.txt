# CMake build script for the IBM ILOG CPLEX solver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND)
  add_library(amplcplex OBJECT cplex.c
    ${PROJECT_SOURCE_DIR}/src/asl/solvers/funcadd0.c)
  target_include_directories(amplcplex PRIVATE
    $<TARGET_PROPERTY:asl,INCLUDE_DIRECTORIES> ${CPLEX_INCLUDE_DIR}/ilcplex)

  if (CMAKE_COMPILER_IS_GNUCC)
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-Wno-unused-result HAVE_WNO_UNUSED_RESULT)
    if (HAVE_WNO_UNUSED_RESULT)
      target_compile_options(amplcplex PRIVATE -Wno-unused-result)
    endif ()
    target_compile_options(amplcplex PRIVATE
      -Wno-deprecated-declarations -Wno-long-long -Wno-sign-compare
      -Wno-uninitialized -Wno-unused-parameter -Wno-missing-field-initializers)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(amplcplex PRIVATE
      -Wno-self-assign -Wno-sign-compare -Wno-sometimes-uninitialized
      -Wno-string-plus-int -Wno-unused-parameter -Wno-missing-field-initializers)
  elseif (MSVC)
    target_compile_options(amplcplex PRIVATE /wd4244 /wd4267 /wd4996)
  endif ()

  add_executable(cplex $<TARGET_OBJECTS:amplcplex>)
  target_link_libraries(cplex cplex-library asl)
endif ()

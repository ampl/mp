# CMake build script for the IBM ILOG CPLEX solver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND)
  include_directories(.. ${CPLEX_INCLUDE_DIRS}/ilcplex)
  if (CMAKE_COMPILER_IS_GNUCC)
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-Wno-unused-result HAVE_WNO_UNUSED_RESULT)
    if (HAVE_WNO_UNUSED_RESULT)
      add_definitions(-Wno-unused-result)
    endif ()
    add_definitions(-Wno-deprecated-declarations -Wno-long-long
                    -Wno-sign-compare -Wno-uninitialized -Wno-unused-parameter
                    -Wno-missing-field-initializers)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_definitions(-Wno-self-assign -Wno-sign-compare
                    -Wno-sometimes-uninitialized -Wno-string-plus-int
                    -Wno-unused-parameter)
  elseif (MSVC)
    add_definitions(/wd4244 /wd4267 /wd4996)
  endif ()
  add_executable(cplex cplex.c ${PROJECT_SOURCE_DIR}/src/asl/solvers/funcadd0.c)
  target_link_libraries(cplex ${CPLEX_LIBRARIES} asl)
endif ()

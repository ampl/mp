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
                    -Wno-sign-compare -Wno-uninitialized -Wno-unused-parameter)
  elseif (MSVC)
    foreach (flag_var
             CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
             CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
      string(REPLACE "/W3" "/w" ${flag_var} "${${flag_var}}")
    endforeach()
  endif ()
  add_executable(cplex cplex.c ../funcadd0.c)
  target_link_libraries(cplex ${CPLEX_LIBRARIES} ${SOLVER_LIBS})
  add_dependencies(cplex arith_h)
endif ()

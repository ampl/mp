# CMake build script for the GSL wrapper for AMPL.

include_directories(${CMAKE_CURRENT_BINARY_DIR}/gsl)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(GSL_DISABLE_WARNINGS 1)
add_subdirectory(gsl)

add_ampl_library(amplgsl amplgsl.c)
add_dependencies(amplgsl gsl gslcblas)
target_link_libraries(amplgsl ${SOLVER_LIBS} gsl gslcblas)

include(CheckLibraryExists)
check_library_exists(m cos "" HAVE_LIBM)
if (HAVE_LIBM)
  target_link_libraries(amplgsl m)
endif ()

add_executable(gen-decls gen-decls.c)
target_link_libraries(gen-decls amplgsl)

add_custom_command(OUTPUT gsl.ampl COMMAND gen-decls DEPENDS amplgsl gen-decls)
add_custom_target(gsl-ampl ALL DEPENDS gsl.ampl)

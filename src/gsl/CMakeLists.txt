# CMake build script for the GSL wrapper for AMPL.

set(GSL_DISABLE_WARNINGS 1)

add_ampl_library(amplgsl amplgsl.c)
target_link_libraries(amplgsl asl gsl gslcblas)
target_include_directories(amplgsl
  PRIVATE ${PROJECT_BINARY_DIR}/thirdparty/build/gsl)

add_executable(gen-decls gen-decls.c)
target_link_libraries(gen-decls amplgsl)

add_custom_command(OUTPUT gsl.ampl COMMAND gen-decls DEPENDS amplgsl gen-decls)
add_custom_target(gsl-ampl ALL DEPENDS gsl.ampl)

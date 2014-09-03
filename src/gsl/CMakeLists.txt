# CMake build script for the GSL wrapper for AMPL.

add_ampl_library(amplgsl amplgsl.c)
target_link_libraries(amplgsl asl gsl gslcblas)
target_include_directories(amplgsl
  PRIVATE ${PROJECT_BINARY_DIR}/thirdparty/build/gsl)

add_executable(gsl-info gsl-info.cc)
target_link_libraries(gsl-info amplgsl)

add_custom_command(OUTPUT gsl.ampl COMMAND gsl-info DEPENDS amplgsl gsl-info)
add_custom_target(gsl-ampl ALL DEPENDS gsl.ampl)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gsl.ampl DESTINATION share/mp)

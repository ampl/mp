# CMake build script for the GSL wrapper for AMPL.

include_directories(${CMAKE_CURRENT_BINARY_DIR}/gsl)

include(CheckCCompilerFlag)
check_c_compiler_flag(-fPIC HAS_FPIC_FLAG)
set(SAVED_FLAGS ${CMAKE_C_FLAGS})
if (HAS_FPIC_FLAG)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif ()
add_subdirectory(gsl)
set(CMAKE_C_FLAGS ${SAVED_FLAGS})

add_library(amplgsl SHARED amplgsl.c)
add_dependencies(amplgsl gsl gslcblas)
target_link_libraries(amplgsl ${SOLVER_LIBS} gsl gslcblas -lm)
set_target_properties(amplgsl PROPERTIES PREFIX "")
set_target_properties(amplgsl PROPERTIES SUFFIX ".dll")

add_executable(gen-decls gen-decls.c)
target_link_libraries(gen-decls amplgsl)

add_custom_command(OUTPUT gsl.ampl
  COMMAND ./gen-decls > gsl.ampl DEPENDS amplgsl gen-decls)
add_custom_target(gsl-ampl ALL DEPENDS gsl.ampl)

add_custom_command(OUTPUT doc/index.html
  COMMAND python extract-docs.py
  COMMAND sphinx-build -b html . doc DEPENDS amplgsl.c conf.py extract-docs.py)
add_custom_target(doc DEPENDS doc/index.html)

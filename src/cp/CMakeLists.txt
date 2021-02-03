# CMake build script for the cp library.
if(build_asl)
  add_ampl_library(cp cp.cc)
  target_link_libraries(cp mp)
endif()
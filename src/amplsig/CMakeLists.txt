# CMake build script for the amplsig library.

if (WIN32)
  add_ampl_library(amplsig amplsig.cc)
  target_link_libraries(amplsig mp Ws2_32.lib)
endif ()

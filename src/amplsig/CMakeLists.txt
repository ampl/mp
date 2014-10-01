# CMake build script for the amplsig library.

if (NOT WIN32)
  # amplsig is only needed on Windows due to broken implementation
  # of signals on this platform.
  return ()
endif ()

find_library(WS2_32_LIBRARY Ws2_32)
if (WS2_32_LIBRARY)
  add_ampl_library(amplsig amplsig.cc)
  target_link_libraries(amplsig mp ${WS2_32_LIBRARY})
endif ()

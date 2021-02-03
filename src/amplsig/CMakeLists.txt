# CMake build script for the amplsig library.
if(NOT build_asl)
  message(STATUS "Not building ASL, skipping amplsig")
  return()
endif()
if (NOT WIN32)
  # amplsig is only needed on Windows due to broken implementation
  # of signals on this platform.
  return ()
endif ()

find_package(WindowsSDK)
GetUMWindowsSKDLibraryDir(WIN_LIBRARY_DIR)

unset(WS2_32_LIBRARY CACHE)
find_library(WS2_32_LIBRARY Ws2_32
PATHS ${WIN_LIBRARY_DIR} NO_DEFAULT_PATH)
if (WS2_32_LIBRARY)
  add_ampl_library(amplsig amplsig.cc)
  target_link_libraries(amplsig mp ${WS2_32_LIBRARY})
endif ()

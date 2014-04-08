if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  # Only 64-bit PATH library for Linux is available.
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    find_library(PATH_LIB path47 PATHS pathlib/lib/linux/dynamic)
    set(PATH_EXTRA_LIBS m)
  endif ()
elseif (WIN32)
  # Only 32-bit PATH library for Windows is available.
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    find_library(PATH_LIB path46 PATHS pathlib/lib/windows/dynamic)
  endif ()
endif ()

if (PATH_LIB)
  include_directories(pathlib/include)
  add_executable(path path.c)
  target_link_libraries(path amplsolver ${PATH_LIB} ${PATH_EXTRA_LIBS})
endif ()

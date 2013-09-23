if (NOT CMAKE_SIZEOF_VOID_P EQUAL 4)
  # PATH only supports 32-bit systems.
elseif (CMAKE_SYSTEM_NAME MATCHES "Linux")
  find_library(PATH_LIB path46 PATHS pathlib/lib/linux/static)
  set(PATH_EXTRA_LIBS m)
elseif (WIN32)
  find_library(PATH_LIB path46 PATHS pathlib/lib/windows/static)
endif ()

if (PATH_LIB)
  include_directories(pathlib/include)
  add_executable(path path.c)
  target_link_libraries(path ${SOLVER_LIBS} ${PATH_LIB} ${PATH_EXTRA_LIBS})
endif ()

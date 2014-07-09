set(PATH_VERSION 47)

if (APPLE)
  set(PATH_LIB osx/libpath${PATH_VERSION}.a)
else ()
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(suffix 64)
  else ()
    set(suffix 32)
  endif ()
  if (WIN32)
    set(PATH_LIB win${suffix}/path${PATH_VERSION}.lib)
  else ()
    set(PATH_LIB linux${suffix}/libpath${PATH_VERSION}.a)
    set(PATH_EXTRA_LIBS m)
  endif ()
endif ()

set(PATH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pathlib)
if (EXISTS ${PATH_DIR}/include)
  include_directories(${PATH_DIR}/include)
  add_executable(path path.c)
  target_link_libraries(path
    amplsolver ${PATH_DIR}/lib/${PATH_LIB} ${PATH_EXTRA_LIBS})
  if (WIN32)
    add_custom_command(TARGET path
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${PATH_DIR}/lib/win${suffix}/path${PATH_VERSION}.dll
              $<TARGET_FILE_DIR:path>/path${PATH_VERSION}.dll)
  endif ()
endif ()

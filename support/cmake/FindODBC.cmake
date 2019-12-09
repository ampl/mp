# Try to find the ODBC libraries.
#
# Once done this will define
#
#  ODBC_FOUND - System has ODBC
#  ODBC_INCLUDE_DIRS - The ODBC include directories
#  ODBC_LIBRARIES - The libraries needed to use ODBC

find_library(ODBC_LIBRARY NAMES odbc odbc32
  PATHS "[${winsdk_key};CurrentInstallFolder]/lib/${libdir_suffix}")

  if(WIN32)
    unset(ODBC_LIBRARY CACHE)
    find_package(WindowsSDK)
    GetUMWindowsSKDLibraryDir(WIN_LIBRARY_DIR)
    GetUMWindowsSKDIncludeDir(WIN_INCLUDE_DIR)
    find_path(ODBC_INCLUDE_DIR sql.h
      PATHS "${WIN_INCLUDE_DIR}" NO_DEFAULT_PATH)
    find_library(ODBC_LIBRARY NAMES odbc odbc32
      PATHS "${WIN_LIBRARY_DIR}" NO_DEFAULT_PATH) 
else()
  # keeping old method on non-windows, can probably be removed altogether
  set(winsdk_key
    "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows")
  find_path(ODBC_INCLUDE_DIR sql.h
    PATHS "[${winsdk_key};CurrentInstallFolder]/include")
  set(libdir_suffix )
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(libdir_suffix x64)
  endif ()
  find_library(ODBC_LIBRARY NAMES odbc odbc32
    PATHS "[${winsdk_key};CurrentInstallFolder]/lib/${libdir_suffix}")
endif()




if (APPLE)
  find_library(COREFOUNDATION_LIBRARY CoreFoundation)
  set(ODBC_EXTRA_LIBS_VAR COREFOUNDATION_LIBRARY)
endif ()

set(ODBC_INCLUDE_DIRS ${ODBC_INCLUDE_DIR})
set(ODBC_LIBRARIES ${ODBC_LIBRARY} ${${ODBC_EXTRA_LIBS_VAR}})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set ODBC_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  ODBC DEFAULT_MSG ODBC_LIBRARY ${ODBC_EXTRA_LIBS_VAR} ODBC_INCLUDE_DIR)

mark_as_advanced(ODBC_LIBRARY ${ODBC_EXTRA_LIBS_VAR} ODBC_INCLUDE_DIR)

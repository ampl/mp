# Try to find the ODBC libraries.
#
# Once done this will define
#
#  ODBC_FOUND - System has ODBC
#  ODBC_INCLUDE_DIRS - The ODBC include directories
#  ODBC_LIBRARIES - The libraries needed to use ODBC

set(winsdk_key
 "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows")
find_path(ODBC_INCLUDE_DIR sql.h
  PATHS "[${winsdk_key};CurrentInstallFolder]/include")
if(ODBC_INCLUDE_DIR)
set(libdir_suffix )
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(libdir_suffix x64)
endif ()
find_library(ODBC_LIBRARY NAMES odbc odbc32
  PATHS "[${winsdk_key};CurrentInstallFolder]/lib/${libdir_suffix}")
else ()  
# If not found (e.g. using windows 10 machines with visual studio)
# TODO REMOVE
	unset(BASE_WINKIT CACHE)
	set(windowskit_key 	"HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows Kits\\Installed Roots")
	find_path(BASE_WINKIT . NO_DEFAULT_PATH PATHS "[${windowskit_key};KitsRoot10]")
	if(NOT BASE_WINKIT)
        message(SEND_ERROR "[ODBC] Could not find Windows SDK")
	endif()
	file(GLOB WINKIT_VERSIONS "${BASE_WINKIT}/include/*")
	# Use newest installed Windows kit
	list(SORT WINKIT_VERSIONS)
    list(REVERSE WINKIT_VERSIONS)
	if(WINKIT_VERSIONS)
	list(GET WINKIT_VERSIONS 0 WINKIT_INCLUDE_DIR)
	message("[ODBC] Using Windows SDK include directory: ${WINKIT_INCLUDE_DIR}")
	get_filename_component(WINKIT_VERSION ${WINKIT_INCLUDE_DIR} NAME)
	set(WINKIT_INCLUDE_DIR ${WINKIT_INCLUDE_DIR}/um)
	set(WINKIT_LIB_DIR ${BASE_WINKIT}/lib/${WINKIT_VERSION}/um)
	
	# Use correct architecture
	set(libdir_suffix )
	if (CMAKE_SIZEOF_VOID_P EQUAL 8)
		set(libdir_suffix x64)
	else ()
		set(libdir_suffix x86)
	endif ()
	
	# Find library and set correct include dirs to correct variable
	find_library(ODBC_LIBRARY NAMES odbc odbc32
	PATHS "${WINKIT_LIB_DIR}/${libdir_suffix}")
	set(ODBC_INCLUDE_DIR ${WINKIT_INCLUDE_DIR})
	endif()
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

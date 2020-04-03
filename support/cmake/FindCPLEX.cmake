# Try to find the CPLEX, Concert, IloCplex and CP Optimizer libraries.
#
# Once done this will add the following imported targets:
#
#  cplex-library - the CPLEX library
#  cplex-concert - the Concert library
#  ilocplex - the IloCplex library
#  cplex-cp - the CP Optimizer library

include(FindPackageHandleStandardArgs)

# Find the path to CPLEX Studio.
# CPLEX Studio 12.4 can be installed in the following default locations:
#   /opt/ibm/ILOG/CPLEX_Studio<edition>124 - Linux
#   /opt/IBM/ILOG/CPLEX_Studio<edition>124 - UNIX
#   ~/Applications/IBM/ILOG/CPLEX_Studio<edition>124 - Mac OS X
#   C:\Program Files\IBM\ILOG\CPLEX_Studio<edition>124 - Windows
if (UNIX)
  set(CPLEX_ILOG_DIRS /opt/ibm/ILOG /opt/IBM/ILOG)
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CPLEX_ARCH x86-64)
  else ()
    set(CPLEX_ARCH x86)
  endif ()
  if (APPLE)
    set(CPLEX_ILOG_DIRS $ENV{HOME}/Applications/IBM/ILOG ${CPLEX_ILOG_DIRS})
    foreach (suffix "osx" "darwin9_gcc4.0")
      set(CPLEX_LIB_PATH_SUFFIXES
          ${CPLEX_LIB_PATH_SUFFIXES} lib/${CPLEX_ARCH}_${suffix}/static_pic)
    endforeach ()
  else ()
    set(CPLEX_LIB_PATH_SUFFIXES
      lib/${CPLEX_ARCH}_sles10_4.1/static_pic lib/${CPLEX_ARCH}_linux/static_pic)
  endif ()
else ()
  set(CPLEX_ILOG_DIRS "C:/Program Files/IBM/ILOG")
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CPLEX_ARCH x64)
  else ()
    set(CPLEX_ARCH x86)
    set(CPLEX_ILOG_DIRS "C:/Program Files (x86)/IBM/ILOG" ${CPLEX_ILOG_DIRS})
  endif ()
   # Amended for VS and its various toolsets
  # https://cmake.org/cmake/help/v3.11/variable/MSVC_VERSION.html
  # Can use GREATER_EQUAL instead of the mess below if cmake version >= 3.7
  if(NOT (MSVC_VERSION LESS 1910)) 
	set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_vs2017/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_vs2017/stat_mdd)
  elseif (NOT (MSVC_VERSION LESS 1900)) # to support VS2015 with 2013 libraries change below
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_vs2015/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_vs2015/stat_mdd)
  elseif (NOT (MSVC_VERSION LESS 1800))
  set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_vs2013/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_vs2013/stat_mdd)
  elseif (NOT (MSVC_VERSION LESS 1700))
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_vs2012/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_vs2012/stat_mdd)
  elseif (NOT (MSVC_VERSION LESS 1600))
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_vs2010/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_vs2010/stat_mdd)
   elseif (NOT (MSVC_VERSION LESS 1500))
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_vs2008/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_vs2008/stat_mdd)
  endif ()
endif ()
if (NOT CPLEX_STUDIO_DIR)
  foreach (dir ${CPLEX_ILOG_DIRS})
    file(GLOB CPLEX_STUDIO_DIRS "${dir}/CPLEX_Studio*")
    list(SORT CPLEX_STUDIO_DIRS)
    list(REVERSE CPLEX_STUDIO_DIRS)
    message("Found studio dirs: ${CPLEX_STUDIO_DIRS}")
    if (CPLEX_STUDIO_DIRS)
      list(GET CPLEX_STUDIO_DIRS 0 CPLEX_STUDIO_DIR_)
      string(REGEX MATCH "[0-9][0-9][0-9][0-9]" CPXVERSION ${CPLEX_STUDIO_DIR_})
      message(STATUS "Found CPLEX Studio: ${CPLEX_STUDIO_DIR_}")
      message(STATUS "Detected CPLEX version ${CPXVERSION}")
      break ()
    endif ()
  endforeach ()
  if (NOT CPLEX_STUDIO_DIR_)
    set(CPLEX_STUDIO_DIR_ CPLEX_STUDIO_DIR-NOTFOUND)
  endif ()
  set(CPLEX_STUDIO_DIR ${CPLEX_STUDIO_DIR_} CACHE PATH
    "Path to the CPLEX Studio directory")
endif ()

# On windows, CPLEX 12.10 brought a big semplification in terms of libraries:
# only one version is used for VS2015, 2017 and 2019 due to the maintained 
# ABI compatibility. Therefore, override the directory
if(MSVC AND (NOT "${CPXVERSION}" STREQUAL ""))
  if(NOT(${CPXVERSION} LESS 1210))
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_windows_msvc14/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG lib/${CPLEX_ARCH}_windows_msvc14/stat_mdd)
  endif()
endif()
find_package(Threads)

# ----------------------------------------------------------------------------
# CPLEX

set(CPLEX_DIR ${CPLEX_STUDIO_DIR_}/cplex)

# Find the CPLEX include directory.
find_path(CPLEX_INCLUDE_DIR ilcplex/cplex.h PATHS ${CPLEX_DIR}/include)

macro(find_win_cplex_library var path_suffixes)
  foreach (s ${path_suffixes})
	file(GLOB CPLEX_LIBRARY_CANDIDATES "${CPLEX_DIR}/${s}/cplex*.lib")
    if (CPLEX_LIBRARY_CANDIDATES)
      list(GET CPLEX_LIBRARY_CANDIDATES 0 ${var})
      break ()
    endif ()
  endforeach ()
  if (NOT ${var})
    set(${var} NOTFOUND)
  endif ()
endmacro()

# Find the CPLEX library.
if (UNIX)
  find_library(CPLEX_LIBRARY NAMES cplex
    PATHS ${CPLEX_DIR} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})
  set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIBRARY})
elseif (NOT CPLEX_LIBRARY)
  # On Windows the version is appended to the library name which cannot be
  # handled by find_library, so search manually.
  find_win_cplex_library(CPLEX_LIB "${CPLEX_LIB_PATH_SUFFIXES}")
  message("CPLEX_LIB " ${CPLEX_LIB})
  message("CPLEX_LIBRARY " ${CPLEX_LIBRARY})
  set(CPLEX_LIBRARY ${CPLEX_LIB} CACHE FILEPATH "Path to the CPLEX library")
  find_win_cplex_library(CPLEX_LIB "${CPLEX_LIB_PATH_SUFFIXES_DEBUG}")
  set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIB} CACHE
    FILEPATH "Path to the debug CPLEX library")
	message("CPLEX_LIBRARY " ${CPLEX_LIBRARY})
  if (CPLEX_LIBRARY MATCHES ".*/(cplex.*)\\.lib")
    file(GLOB CPLEX_DLL_ "${CPLEX_DIR}/bin/*/${CMAKE_MATCH_1}.dll")
    set(CPLEX_DLL ${CPLEX_DLL_} CACHE PATH "Path to the CPLEX DLL.")
	message("CPLEX_DLL " ${CPLEX_DLL})
  endif ()
endif ()

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR)

mark_as_advanced(CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR)

if (CPLEX_FOUND AND NOT TARGET cplex-library)
  set(CPLEX_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})
  check_library_exists(m floor "" HAVE_LIBM)
  if (HAVE_LIBM)
    set(CPLEX_LINK_LIBRARIES ${CPLEX_LINK_LIBRARIES} m)
  endif ()
  add_library(cplex-library STATIC IMPORTED GLOBAL)
  set_target_properties(cplex-library PROPERTIES
    IMPORTED_LOCATION "${CPLEX_LIBRARY}"
    IMPORTED_LOCATION_DEBUG "${CPLEX_LIBRARY_DEBUG}"
    INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${CPLEX_LINK_LIBRARIES}")
endif ()

# ----------------------------------------------------------------------------
# Concert

macro(find_cplex_library var name paths)
  find_library(${var} NAMES ${name}
    PATHS ${paths} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})
  if (UNIX)
    set(${var}_DEBUG ${${var}})
  else ()
    find_library(${var}_DEBUG NAMES ${name}
      PATHS ${paths} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES_DEBUG})
  endif ()
endmacro()

set(CPLEX_CONCERT_DIR ${CPLEX_STUDIO_DIR}/concert)

# Find the Concert include directory.
find_path(CPLEX_CONCERT_INCLUDE_DIR ilconcert/ilosys.h
  PATHS ${CPLEX_CONCERT_DIR}/include)

# Find the Concert library.
find_cplex_library(CPLEX_CONCERT_LIBRARY concert ${CPLEX_CONCERT_DIR})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CONCERT_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CONCERT DEFAULT_MSG CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_DEBUG
  CPLEX_CONCERT_INCLUDE_DIR)

mark_as_advanced(CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_DEBUG
  CPLEX_CONCERT_INCLUDE_DIR)

if (CPLEX_CONCERT_FOUND AND NOT TARGET cplex-concert)
  add_library(cplex-concert STATIC IMPORTED GLOBAL)
  set_target_properties(cplex-concert PROPERTIES
    IMPORTED_LOCATION "${CPLEX_CONCERT_LIBRARY}"
    IMPORTED_LOCATION_DEBUG "${CPLEX_CONCERT_LIBRARY_DEBUG}"
    INTERFACE_COMPILE_DEFINITIONS IL_STD # Require standard compliance.
    INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_CONCERT_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
endif ()

# ----------------------------------------------------------------------------
# IloCplex - depends on CPLEX and Concert

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAVE_WNO_LONG_LONG_FLAG)
if (HAVE_WNO_LONG_LONG_FLAG)
  # Required if -pedantic is used.
  set(CPLEX_ILOCPLEX_DEFINITIONS -Wno-long-long)
endif ()

# Find the IloCplex include directory - normally the same as the one for CPLEX
# but check if ilocplex.h is there anyway.
find_path(CPLEX_ILOCPLEX_INCLUDE_DIR ilcplex/ilocplex.h
  PATHS ${CPLEX_INCLUDE_DIR})

# Find the IloCplex library.
find_cplex_library(CPLEX_ILOCPLEX_LIBRARY ilocplex ${CPLEX_DIR})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_ILOCPLEX_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_ILOCPLEX DEFAULT_MSG
  CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_DEBUG
  CPLEX_ILOCPLEX_INCLUDE_DIR CPLEX_FOUND CPLEX_CONCERT_FOUND)

mark_as_advanced(CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_DEBUG
  CPLEX_ILOCPLEX_INCLUDE_DIR)

if (CPLEX_ILOCPLEX_FOUND AND NOT TARGET ilocplex)
  add_library(ilocplex STATIC IMPORTED GLOBAL)
  set_target_properties(ilocplex PROPERTIES
    IMPORTED_LOCATION "${CPLEX_ILOCPLEX_LIBRARY}"
    IMPORTED_LOCATION_DEBUG "${CPLEX_ILOCPLEX_LIBRARY_DEBUG}"
    INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_ILOCPLEX_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "cplex-concert;cplex-library")
endif ()

# ----------------------------------------------------------------------------
# CP Optimizer - depends on Concert

set(CPLEX_CP_DIR ${CPLEX_STUDIO_DIR}/cpoptimizer)

# Find the CP Optimizer include directory.
find_path(CPLEX_CP_INCLUDE_DIR ilcp/cp.h PATHS ${CPLEX_CP_DIR}/include)

# Find the CP Optimizer library.
find_cplex_library(CPLEX_CP_LIBRARY cp ${CPLEX_CP_DIR})

if (WIN32)
  set(CPLEX_CP_EXTRA_LIBRARIES Ws2_32.lib)
endif ()

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CP_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CP DEFAULT_MSG CPLEX_CP_LIBRARY CPLEX_CP_LIBRARY_DEBUG
  CPLEX_CP_INCLUDE_DIR CPLEX_CONCERT_FOUND)

mark_as_advanced(CPLEX_CP_LIBRARY CPLEX_CP_LIBRARY_DEBUG CPLEX_CP_INCLUDE_DIR)

if (CPLEX_CP_FOUND AND NOT TARGET cplex-cp)
  add_library(cplex-cp STATIC IMPORTED GLOBAL)
  set_target_properties(cplex-cp PROPERTIES
    IMPORTED_LOCATION "${CPLEX_CP_LIBRARY}"
    IMPORTED_LOCATION_DEBUG "${CPLEX_CP_LIBRARY_DEBUG}"
    INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_CP_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "cplex-concert;${CPLEX_CP_EXTRA_LIBRARIES}")
endif ()

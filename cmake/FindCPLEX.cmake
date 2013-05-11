# Try to find the CPLEX, Concert, IloCplex and CP Optimizer libraries.
#
# Once done this will define
#
#  CPLEX_FOUND - System has CPLEX
#  CPLEX_INCLUDE_DIRS - The CPLEX include directories
#  CPLEX_LIBRARIES - The libraries needed to use CPLEX
#
#  CPLEX_ILOCPLEX_FOUND - System has IloCplex
#  CPLEX_ILOCPLEX_INCLUDE_DIRS - The IloCplex include directories
#  CPLEX_ILOCPLEX_LIBRARIES - The libraries needed to use IloCplex
#  CPLEX_ILOCPLEX_DEFINITIONS - Compiler switches required for using IloCplex
#
#  CPLEX_CONCERT_FOUND - System has Concert
#  CPLEX_CONCERT_INCLUDE_DIRS - The Concert include directories
#  CPLEX_CONCERT_LIBRARIES - The libraries needed to use Concert
#  CPLEX_CONCERT_DEFINITIONS - Compiler switches required for using Concert
#
#  CPLEX_CP_FOUND - System has CP Optimizer
#  CPLEX_CP_INCLUDE_DIRS - The CP Optimizer include directories
#  CPLEX_CP_LIBRARIES - The libraries needed to use CP Optimizer

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
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_darwin9_gcc4.0/static_pic)
  else ()
    set(CPLEX_LIB_PATH_SUFFIXES lib/${CPLEX_ARCH}_sles10_4.1/static_pic)
  endif ()
else ()
  set(CPLEX_ILOG_DIRS "C:/Program Files/IBM/ILOG")
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CPLEX_ARCH x64)
  else ()
    set(CPLEX_ARCH x86)
  endif ()
  if (MSVC10)
    set(CPLEX_LIB_PATH_SUFFIXES
      lib/${CPLEX_ARCH}_windows_vs2010/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
      lib/${CPLEX_ARCH}_windows_vs2010/stat_mdd)
  elseif (MSVC9)
    set(CPLEX_LIB_PATH_SUFFIXES
      lib/${CPLEX_ARCH}_windows_vs2008/stat_mda)
    set(CPLEX_LIB_PATH_SUFFIXES_DEBUG
      lib/${CPLEX_ARCH}_windows_vs2008/stat_mdd)
  endif ()
endif ()
if (NOT CPLEX_STUDIO_DIR)
  foreach (dir ${CPLEX_ILOG_DIRS})
    file(GLOB CPLEX_STUDIO_DIRS "${dir}/CPLEX_Studio*")
    if (CPLEX_STUDIO_DIRS)
      list(GET CPLEX_STUDIO_DIRS 0 CPLEX_STUDIO_DIR_)
      message(STATUS "Found CPLEX Studio: ${CPLEX_STUDIO_DIR_}")
      break ()
    endif ()
  endforeach ()
  if (NOT CPLEX_STUDIO_DIR_)
    set(CPLEX_STUDIO_DIR_ CPLEX_STUDIO_DIR-NOTFOUND)
  endif ()
  set(CPLEX_STUDIO_DIR ${CPLEX_STUDIO_DIR_} CACHE PATH
    "Path to the CPLEX Studio directory")
endif ()

find_package(Threads)

# ----------------------------------------------------------------------------
# CPLEX

set(CPLEX_DIR ${CPLEX_STUDIO_DIR}/cplex)

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
  set(CPLEX_LIBRARY ${CPLEX_LIB} CACHE FILEPATH "Path to the CPLEX library")
  find_win_cplex_library(CPLEX_LIB "${CPLEX_LIB_PATH_SUFFIXES_DEBUG}")
  set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIB} CACHE
    FILEPATH "Path to the debug CPLEX library")
endif ()

set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
set(CPLEX_LIBRARIES
  optimized ${CPLEX_LIBRARY} debug ${CPLEX_LIBRARY_DEBUG}
  ${CMAKE_THREAD_LIBS_INIT})
check_library_exists(m floor "" HAVE_LIBM)
if (HAVE_LIBM)
  set(CPLEX_LIBRARIES ${CPLEX_LIBRARIES} m)
endif ()

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR)

mark_as_advanced(CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR)

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

# Require standard compliance.
set(CPLEX_CONCERT_DEFINITIONS -DIL_STD)

# Find the Concert include directory.
find_path(CPLEX_CONCERT_INCLUDE_DIR ilconcert/ilosys.h
  PATHS ${CPLEX_CONCERT_DIR}/include)

# Find the Concert library.
find_cplex_library(CPLEX_CONCERT_LIBRARY concert ${CPLEX_CONCERT_DIR})

set(CPLEX_CONCERT_INCLUDE_DIRS ${CPLEX_CONCERT_INCLUDE_DIR})
set(CPLEX_CONCERT_LIBRARIES
  optimized ${CPLEX_CONCERT_LIBRARY} debug ${CPLEX_CONCERT_LIBRARY_DEBUG}
  ${CMAKE_THREAD_LIBS_INIT})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CONCERT_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CONCERT DEFAULT_MSG CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_DEBUG
  CPLEX_CONCERT_INCLUDE_DIR)

mark_as_advanced(CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_DEBUG
  CPLEX_CONCERT_INCLUDE_DIR)

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

set(CPLEX_ILOCPLEX_INCLUDE_DIRS
  ${CPLEX_ILOCPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIRS})
set(CPLEX_ILOCPLEX_LIBRARIES
  optimized ${CPLEX_ILOCPLEX_LIBRARY} debug ${CPLEX_ILOCPLEX_LIBRARY_DEBUG}
  ${CPLEX_CONCERT_LIBRARIES} ${CPLEX_LIBRARIES})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_ILOCPLEX_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_ILOCPLEX DEFAULT_MSG
  CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_DEBUG
  CPLEX_ILOCPLEX_INCLUDE_DIR CPLEX_FOUND CPLEX_CONCERT_FOUND)

mark_as_advanced(CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_DEBUG
  CPLEX_ILOCPLEX_INCLUDE_DIR)

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

set(CPLEX_CP_INCLUDE_DIRS
  ${CPLEX_CP_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIRS})
set(CPLEX_CP_LIBRARIES
  optimized ${CPLEX_CP_LIBRARY} debug ${CPLEX_CP_LIBRARY_DEBUG}
  ${CPLEX_CONCERT_LIBRARIES} ${CPLEX_CP_EXTRA_LIBRARIES})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CP_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CP DEFAULT_MSG CPLEX_CP_LIBRARY CPLEX_CP_LIBRARY_DEBUG
  CPLEX_CP_INCLUDE_DIR CPLEX_CONCERT_FOUND)

mark_as_advanced(CPLEX_CP_LIBRARY CPLEX_CP_LIBRARY_DEBUG CPLEX_CP_INCLUDE_DIR)

# Try to find the CPLEX, Concert, IloCplex and CP Optimizer libraries.
#
# Once done this will define
#
#  CPLEX_FOUND - System has CPLEX
#  CPLEX_INCLUDE_DIRS - The CPLEX include directories
#  CPLEX_LIBRARIES - The libraries needed to use CPLEX
#  CPLEX_LIBRARIES_DEBUG - The debug libraries needed to use CPLEX
#
#  CPLEX_ILOCPLEX_FOUND - System has IloCplex
#  CPLEX_ILOCPLEX_INCLUDE_DIRS - The IloCplex include directories
#  CPLEX_ILOCPLEX_LIBRARIES - The libraries needed to use IloCplex
#  CPLEX_ILOCPLEX_LIBRARIES_DEBUG - The debug libraries needed to use IloCplex
#  CPLEX_ILOCPLEX_DEFINITIONS - Compiler switches required for using IloCplex
#
#  CPLEX_CONCERT_FOUND - System has Concert
#  CPLEX_CONCERT_INCLUDE_DIRS - The Concert include directories
#  CPLEX_CONCERT_LIBRARIES - The libraries needed to use Concert
#  CPLEX_CONCERT_LIBRARIES_DEBUG - The debug libraries needed to use Concert
#  CPLEX_CONCERT_DEFINITIONS - Compiler switches required for using Concert
#
#  CPLEX_CP_FOUND - System has CP Optimizer
#  CPLEX_CP_INCLUDE_DIRS - The CP Optimizer include directories
#  CPLEX_CP_LIBRARIES - The libraries needed to use CP Optimizer
#  CPLEX_CP_LIBRARIES_DEBUG - The debug libraries needed to use CP Optimizer

include(FindPackageHandleStandardArgs)

# Recent versions of CPLEX Studio are installed in the following locations:
#   /opt/ibm/ILOG/CPLEX_Studio<version> - Linux (checked version 12.4)
#   C:\ILOG\CPLEX_Studio<version> - Windows (checked version 12.2)
#   C:\Program Files\IBM\ILOG\CPLEX_Studio<version> - Windows (README)
if (UNIX)
  set(CPLEX_ILOG_DIRS
      /opt/ibm/ILOG /opt/IBM/ILOG $ENV{HOME}/ILOG $ENV{HOME}/ilog)
  if (APPLE)
    set(CPLEX_LIB_PATH_SUFFIXES lib/x86_darwin9_gcc4.0/static_pic)
  else ()
    set(CPLEX_LIB_PATH_SUFFIXES
      lib/x86_sles10_4.1/static_pic
      lib/x86-64_sles10_4.1/static_pic)
  endif ()
else ()
  set(CPLEX_ILOG_DIRS C:/ILOG "C:/Program Files/IBM/ILOG")
  if (MSVC10)
    set(CPLEX_LIB_PATH_SUFFIXES lib/x86_windows_vs2010/stat_mda)
    set(CPLEX_DEBUG_LIB_PATH_SUFFIXES lib/x86_windows_vs2010/stat_mdd)
  elseif (MSVC9)
    set(CPLEX_LIB_PATH_SUFFIXES lib/x86_windows_vs2008/stat_mda)
    set(CPLEX_DEBUG_LIB_PATH_SUFFIXES lib/x86_windows_vs2008/stat_mdd)
  endif ()
endif ()
foreach (d ${CPLEX_ILOG_DIRS})
  if (EXISTS ${d})
    set(CPLEX_STUDIO_PATH ${d})
    break ()
  endif ()
endforeach ()

find_package(Threads)

# ----------------------------------------------------------------------------
# CPLEX

file(GLOB CPLEX_PATHS "${CPLEX_STUDIO_PATH}/*/cplex")

# Find the CPLEX include directory.
find_path(CPLEX_INCLUDE_DIR ilcplex/cplex.h
  PATHS ${CPLEX_PATHS} PATH_SUFFIXES include)

# Find the CPLEX library.
macro(find_cplex_library var path_suffixes)
  if (UNIX)
    find_library(${var} NAMES cplex
      PATHS ${CPLEX_PATHS} PATH_SUFFIXES ${path_suffixes})
  elseif (NOT ${var})
    # On Windows the version is appended to the library name which cannot be
    # handled by find_library, so search manually.
    foreach (p ${CPLEX_PATHS})
      foreach (s ${path_suffixes})
        file(GLOB CPLEX_LIBRARY_CANDIDATES "${p}/${s}/cplex*.lib")
        if (CPLEX_LIBRARY_CANDIDATES)
          list(GET CPLEX_LIBRARY_CANDIDATES 0 CPLEX_LIB)
          set(${var} ${CPLEX_LIB}
            CACHE FILEPATH "Path to the CPLEX library")
          break ()
        endif ()
      endforeach ()
      if (${var})
        break ()
      endif ()
    endforeach ()
  endif ()
endmacro ()

find_cplex_library(CPLEX_LIBRARY "${CPLEX_LIB_PATH_SUFFIXES}")
if (UNIX)
  set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIBRARY})
else ()
  find_cplex_library(CPLEX_LIBRARY_DEBUG ${CPLEX_DEBUG_LIB_PATH_SUFFIXES})
endif ()

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG CPLEX_INCLUDE_DIR)

set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
set(CPLEX_LIBRARIES ${CPLEX_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
set(CPLEX_LIBRARIES_DEBUG ${CPLEX_LIBRARY_DEBUG} ${CMAKE_THREAD_LIBS_INIT})

mark_as_advanced(CPLEX_INCLUDE_DIR CPLEX_LIBRARY CPLEX_LIBRARY_DEBUG)

# ----------------------------------------------------------------------------
# Concert

macro(find_cplex_cxx_library var name paths)
  find_library(${var} NAMES ${name}
    PATHS ${paths} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})
  if (UNIX)
    set(${var}_DEBUG ${${var}})
  else ()
    find_library(${var}_DEBUG NAMES ${name}
      PATHS ${paths}
      PATH_SUFFIXES ${CPLEX_DEBUG_LIB_PATH_SUFFIXES})
  endif ()
endmacro()

file(GLOB CPLEX_CONCERT_PATHS "${CPLEX_STUDIO_PATH}/*/concert")

# Require standard compliance.
set(CPLEX_CONCERT_DEFINITIONS -DIL_STD)

# Find the Concert include directory.
find_path(CPLEX_CONCERT_INCLUDE_DIR ilconcert/ilosys.h
  PATHS ${CPLEX_CONCERT_PATHS} PATH_SUFFIXES include)

# Find the Concert library.
find_cplex_cxx_library(CPLEX_CONCERT_LIBRARY concert "${CPLEX_CONCERT_PATHS}")

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CONCERT_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CONCERT DEFAULT_MSG CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_DEBUG
  CPLEX_CONCERT_INCLUDE_DIR)

set(CPLEX_CONCERT_INCLUDE_DIRS ${CPLEX_CONCERT_INCLUDE_DIR})
set(CPLEX_CONCERT_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
set(CPLEX_CONCERT_LIBRARIES_DEBUG ${CPLEX_CONCERT_LIBRARY_DEBUG}
  ${CMAKE_THREAD_LIBS_INIT})

mark_as_advanced(CPLEX_CONCERT_INCLUDE_DIR CPLEX_CONCERT_LIBRARY
  CPLEX_CONCERT_LIBRARY_DEBUG)

# ----------------------------------------------------------------------------
# IloCplex - depends on CPLEX and Concert

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-long-long HAS_WNO_LONG_LONG_FLAG)
if (HAS_WNO_LONG_LONG_FLAG)
  # Required if -pedantic is used.
  set(CPLEX_ILOCPLEX_DEFINITIONS -Wno-long-long)
endif ()

# Find the IloCplex include directory - normally the same as the one for CPLEX
# but check if ilocplex.h is there anyway.
find_path(CPLEX_ILOCPLEX_INCLUDE_DIR ilcplex/ilocplex.h
  PATHS ${CPLEX_INCLUDE_DIR})

# Find the IloCplex library.
find_cplex_cxx_library(CPLEX_ILOCPLEX_LIBRARY ilocplex "${CPLEX_PATHS}")

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_ILOCPLEX_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_ILOCPLEX DEFAULT_MSG
  CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_DEBUG
  CPLEX_ILOCPLEX_INCLUDE_DIR CPLEX_FOUND CPLEX_CONCERT_FOUND)

set(CPLEX_ILOCPLEX_INCLUDE_DIRS
  ${CPLEX_ILOCPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIRS})
set(CPLEX_ILOCPLEX_LIBRARIES
  ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_CONCERT_LIBRARIES} ${CPLEX_LIBRARIES})
set(CPLEX_ILOCPLEX_LIBRARIES_DEBUG
  ${CPLEX_ILOCPLEX_LIBRARY_DEBUG} ${CPLEX_CONCERT_LIBRARIES_DEBUG}
  ${CPLEX_LIBRARIES})

mark_as_advanced(CPLEX_ILOCPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY
  CPLEX_ILOCPLEX_LIBRARY_DEBUG)

# ----------------------------------------------------------------------------
# CP Optimizer - depends on Concert

file(GLOB CPLEX_CP_PATHS "${CPLEX_STUDIO_PATH}/*/cpoptimizer")

# Find the CP Optimizer include directory.
find_path(CPLEX_CP_INCLUDE_DIR ilcp/cp.h
  PATHS ${CPLEX_CP_PATHS} PATH_SUFFIXES include)

# Find the CP Optimizer library.
find_cplex_cxx_library(CPLEX_CP_LIBRARY cp "${CPLEX_CP_PATHS}")

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CP_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CP DEFAULT_MSG CPLEX_CP_LIBRARY CPLEX_CP_LIBRARY_DEBUG
  CPLEX_CP_INCLUDE_DIR CPLEX_CONCERT_FOUND)

if (WIN32)
  set(CPLEX_CP_EXTRA_LIBRARIES Ws2_32.lib)
endif ()

set(CPLEX_CP_INCLUDE_DIRS
  ${CPLEX_CP_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIRS})
set(CPLEX_CP_LIBRARIES
  ${CPLEX_CP_LIBRARY} ${CPLEX_CONCERT_LIBRARIES}
  ${CPLEX_CP_EXTRA_LIBRARIES})
set(CPLEX_CP_LIBRARIES_DEBUG
  ${CPLEX_CP_LIBRARY_DEBUG} ${CPLEX_CONCERT_LIBRARIES_DEBUG}
  ${CPLEX_CP_EXTRA_LIBRARIES})

mark_as_advanced(CPLEX_CP_INCLUDE_DIR CPLEX_CP_LIBRARY CPLEX_CP_LIBRARY_DEBUG)

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
    set(CPLEX_LIB_PATH_SUFFIXES
      lib/x86_windows_vs2010/stat_mdd
      lib/x86_windows_vs2010/stat_mda)
  elseif (MSVC9)
    set(CPLEX_LIB_PATH_SUFFIXES
      lib/x86_windows_vs2008/stat_mda
      lib/x86_windows_vs2008/stat_mta)
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
if (NOT WIN32)
  find_library(CPLEX_LIBRARY NAMES cplex
    PATHS ${CPLEX_PATHS} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})
elseif (NOT CPLEX_LIBRARY)
  # On Windows the version is appended to the library name which cannot be
  # handled by find_library, so search manually.
  foreach (p ${CPLEX_PATHS})
    foreach (s ${CPLEX_LIB_PATH_SUFFIXES})
      file(GLOB CPLEX_LIBRARY_CANDIDATES "${p}/${s}/cplex*.lib")
      if (CPLEX_LIBRARY_CANDIDATES)
        list(GET CPLEX_LIBRARY_CANDIDATES 0 CPLEX_LIB)
        set(CPLEX_LIBRARY ${CPLEX_LIB}
          CACHE FILEPATH "Path to the CPLEX library")
        break ()
      endif ()
    endforeach ()
    if (CPLEX_LIBRARY)
      break ()
    endif ()
  endforeach ()
endif ()

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
set(CPLEX_LIBRARIES ${CPLEX_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

mark_as_advanced(CPLEX_INCLUDE_DIR CPLEX_LIBRARY)

# ----------------------------------------------------------------------------
# Concert

file(GLOB CPLEX_CONCERT_PATHS "${CPLEX_STUDIO_PATH}/*/concert")

# Require standard compliance.
set(CPLEX_CONCERT_DEFINITIONS -DIL_STD)

# Find the Concert include directory.
find_path(CPLEX_CONCERT_INCLUDE_DIR ilconcert/ilosys.h
  PATHS ${CPLEX_CONCERT_PATHS} PATH_SUFFIXES include)

# Find the Concert library.
find_library(CPLEX_CONCERT_LIBRARY NAMES concert
  PATHS ${CPLEX_CONCERT_PATHS} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CONCERT_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CONCERT DEFAULT_MSG CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_INCLUDE_DIR)

set(CPLEX_CONCERT_INCLUDE_DIRS ${CPLEX_CONCERT_INCLUDE_DIR})
set(CPLEX_CONCERT_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

mark_as_advanced(CPLEX_CONCERT_INCLUDE_DIR CPLEX_CONCERT_LIBRARY)

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
find_library(CPLEX_ILOCPLEX_LIBRARY NAMES ilocplex
  PATHS ${CPLEX_PATHS} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_ILOCPLEX_FOUND to
# TRUE if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_ILOCPLEX DEFAULT_MSG
  CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_INCLUDE_DIR
  CPLEX_FOUND CPLEX_CONCERT_FOUND)

set(CPLEX_ILOCPLEX_INCLUDE_DIRS
  ${CPLEX_ILOCPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIRS})
set(CPLEX_ILOCPLEX_LIBRARIES
  ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_CONCERT_LIBRARIES} ${CPLEX_LIBRARIES})

mark_as_advanced(CPLEX_ILOCPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY)

# ----------------------------------------------------------------------------
# CP Optimizer - depends on Concert

file(GLOB CPLEX_CP_PATHS "${CPLEX_STUDIO_PATH}/*/cpoptimizer")

# Find the CP Optimizer include directory.
find_path(CPLEX_CP_INCLUDE_DIR ilcp/cp.h
  PATHS ${CPLEX_CP_PATHS} PATH_SUFFIXES include)

# Find the CP Optimizer library.
find_library(CPLEX_CP_LIBRARY NAMES cp
  PATHS ${CPLEX_CP_PATHS} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_CP_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX_CP DEFAULT_MSG
  CPLEX_CP_LIBRARY CPLEX_CP_INCLUDE_DIR CPLEX_CONCERT_FOUND)

set(CPLEX_CP_INCLUDE_DIRS
  ${CPLEX_CP_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIRS})
set(CPLEX_CP_LIBRARIES
  ${CPLEX_CP_LIBRARY} ${CPLEX_CONCERT_LIBRARIES})
if (WIN32)
  set(CPLEX_CP_LIBRARIES ${CPLEX_CP_LIBRARIES} Ws2_32.lib)
endif ()

mark_as_advanced(CPLEX_CP_INCLUDE_DIR CPLEX_CP_LIBRARY)

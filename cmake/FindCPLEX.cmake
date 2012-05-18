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
#
#  CPLEX_CONCERT_FOUND - System has Concert
#  CPLEX_CONCERT_INCLUDE_DIRS - The Concert include directories
#  CPLEX_CONCERT_LIBRARIES - The libraries needed to use Concert
#
#  CPLEX_CP_FOUND - System has CP Optimizer
#  CPLEX_CP_INCLUDE_DIRS - The CP Optimizer include directories
#  CPLEX_CP_LIBRARIES - The libraries needed to use CP Optimizer

include(FindPackageHandleStandardArgs)

set(CPLEX_STUDIO_PATH /opt/ibm/ILOG)
set(CPLEX_LIB_PATH_SUFFIXES lib/x86-64_sles10_4.1/static_pic)

# ----------------------------------------------------------------------------
# CPLEX

file(GLOB CPLEX_PATHS "${CPLEX_STUDIO_PATH}/*/cplex")

# Find the CPLEX include directory.
find_path(CPLEX_INCLUDE_DIR ilcplex/cplex.h
  PATHS ${CPLEX_PATHS} PATH_SUFFIXES include)

# Find the CPLEX library.
find_library(CPLEX_LIBRARY NAMES cplex
  PATHS ${CPLEX_PATHS} PATH_SUFFIXES ${CPLEX_LIB_PATH_SUFFIXES})

# Handle the QUIETLY and REQUIRED arguments and set CPLEX_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
set(CPLEX_LIBRARIES ${CPLEX_LIBRARY})

mark_as_advanced(CPLEX_INCLUDE_DIR CPLEX_LIBRARY)

# ----------------------------------------------------------------------------
# Concert

file(GLOB CPLEX_CONCERT_PATHS "${CPLEX_STUDIO_PATH}/*/concert")

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
set(CPLEX_CONCERT_LIBRARIES ${CPLEX_CONCERT_LIBRARY})

mark_as_advanced(CPLEX_CONCERT_INCLUDE_DIR CPLEX_CONCERT_LIBRARY)

# ----------------------------------------------------------------------------
# IloCplex - depends on CPLEX and Concert

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

mark_as_advanced(CPLEX_CP_INCLUDE_DIR CPLEX_CP_LIBRARY)

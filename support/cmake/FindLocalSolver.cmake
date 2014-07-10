# Try to find the LocalSolver libraries.
#
# Once done this will define
#
#  LOCALSOLVER_FOUND - System has LocalSolver
#  LOCALSOLVER_INCLUDE_DIRS - The LocalSolver include directories
#  LOCALSOLVER_LIBRARIES - The libraries needed to use LocalSolver

if (UNIX)
  set(LOCALSOLVER_DIR /opt)
else ()
  set(LOCALSOLVER_DIR "C:")
endif ()

file(GLOB LOCALSOLVER_DIRS "${LOCALSOLVER_DIR}/localsolver*")
if (LOCALSOLVER_DIRS)
  list(GET LOCALSOLVER_DIRS 0 LOCALSOLVER_DIR)
  message(STATUS "Found LocalSolver directory: ${LOCALSOLVER_DIR}")
endif ()

find_path(LOCALSOLVER_INCLUDE_DIR localsolver.h PATHS ${LOCALSOLVER_DIR}/include)
find_library(LOCALSOLVER_LIBRARY localsolver PATHS ${LOCALSOLVER_DIR}/bin)

set(LOCALSOLVER_INCLUDE_DIRS ${LOCALSOLVER_INCLUDE_DIR})
set(LOCALSOLVER_LIBRARIES ${LOCALSOLVER_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set ODBC_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(LOCALSOLVER DEFAULT_MSG
  LOCALSOLVER_LIBRARY LOCALSOLVER_INCLUDE_DIR)

mark_as_advanced(LOCALSOLVER_LIBRARY LOCALSOLVER_INCLUDE_DIR)

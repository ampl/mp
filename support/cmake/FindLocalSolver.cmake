# Try to find the LocalSolver libraries.
#
# Once done this will define
#
#  LocalSolver_FOUND - System has LocalSolver
#  LocalSolver_INCLUDE_DIRS - The LocalSolver include directories
#  LocalSolver_LIBRARIES - The libraries needed to use LocalSolver

if (UNIX)
  set(LocalSolver_DIR /opt)
else ()
  set(LocalSolver_DIR "C:")
endif ()

file(GLOB LocalSolver_DIRS "${LocalSolver_DIR}/localsolver*")
if (LocalSolver_DIRS)
  list(GET LocalSolver_DIRS 0 LocalSolver_DIR)
  message(STATUS "Found LocalSolver directory: ${LocalSolver_DIR}")
endif ()

find_path(LocalSolver_INCLUDE_DIR localsolver.h PATHS ${LocalSolver_DIR}/include)
find_library(LocalSolver_LIBRARY localsolver PATHS ${LocalSolver_DIR}/bin)

set(LocalSolver_INCLUDE_DIRS ${LocalSolver_INCLUDE_DIR})
set(LocalSolver_LIBRARIES ${LocalSolver_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set ODBC_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(LocalSolver DEFAULT_MSG
  LocalSolver_LIBRARY LocalSolver_INCLUDE_DIR)

mark_as_advanced(LocalSolver_LIBRARY LocalSolver_INCLUDE_DIR)

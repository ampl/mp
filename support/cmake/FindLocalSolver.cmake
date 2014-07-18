# Try to find the LocalSolver libraries.
#
# Once done this will define
#
#  LOCALSOLVER_FOUND - System has LocalSolver
#  LOCALSOLVER_INCLUDE_DIRS - The LocalSolver include directories
#  LOCALSOLVER_LIBRARIES - The libraries needed to use LocalSolver

if (UNIX)
  set(LOCALSOLVER_DIR /opt)
  set(LOCALSOLVER_LIB_NAME localsolver)
else ()
  set(PROGRAM_FILES_DIR "C:/Program Files")
  set(PROGRAM_FILES_X86_DIR "${PROGRAM_FILES_DIR} (x86)")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4 AND EXISTS ${PROGRAM_FILES_X86_DIR})
    set(PROGRAM_FILES_DIR ${PROGRAM_FILES_X86_DIR})
  endif ()
  set(LOCALSOLVER_DIR ${PROGRAM_FILES_DIR})
  set(LOCALSOLVER_LIB_NAME localsolver.dll)
endif ()

file(GLOB LOCALSOLVER_DIRS "${LOCALSOLVER_DIR}/localsolver*")
if (LOCALSOLVER_DIRS)
  list(GET LOCALSOLVER_DIRS 0 LOCALSOLVER_DIR)
  message(STATUS "Found LocalSolver directory: ${LOCALSOLVER_DIR}")
endif ()

find_path(LOCALSOLVER_INCLUDE_DIR
  localsolver.h PATHS ${LOCALSOLVER_DIR}/include)
find_library(LOCALSOLVER_LIBRARY
  ${LOCALSOLVER_LIB_NAME} PATHS ${LOCALSOLVER_DIR}/bin)

set(LOCALSOLVER_INCLUDE_DIRS ${LOCALSOLVER_INCLUDE_DIR})
set(LOCALSOLVER_LIBRARIES ${LOCALSOLVER_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set LOCALSOLVER_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(LocalSolver DEFAULT_MSG
  LOCALSOLVER_LIBRARY LOCALSOLVER_INCLUDE_DIR)

mark_as_advanced(LOCALSOLVER_LIBRARY LOCALSOLVER_INCLUDE_DIR)

# Try to find the LocalSolver library.
#
# Once done this will add the following imported targets:
#
#  localsolver-library - the LocalSolver library

if (UNIX)
  set(INSTALL_DIRS /opt)
  set(LOCALSOLVER_LIB_NAME localsolver)
else ()
  set(PROGRAM_FILES_DIR "C:/Program Files")
  set(PROGRAM_FILES_X86_DIR "${PROGRAM_FILES_DIR} (x86)")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4 AND EXISTS ${PROGRAM_FILES_X86_DIR})
    set(PROGRAM_FILES_DIR ${PROGRAM_FILES_X86_DIR})
  endif ()
  set(INSTALL_DIRS ${PROGRAM_FILES_DIR} "C:")
  set(LOCALSOLVER_LIB_NAME localsolver.dll)
endif ()

foreach (dir ${INSTALL_DIRS})
  file(GLOB LOCALSOLVER_DIRS "${dir}/localsolver*")
  if (LOCALSOLVER_DIRS)
    list(GET LOCALSOLVER_DIRS 0 LOCALSOLVER_DIR)
    message(STATUS "Found LocalSolver directory: ${LOCALSOLVER_DIR}")
    if (WIN32)
      set(LOCALSOLVER_DLL ${LOCALSOLVER_DIR}/bin/localsolver.dll
        CACHE PATH "Path to the LocalSolver DLL.")
    endif ()
    break ()
  endif ()
endforeach ()

find_path(LOCALSOLVER_INCLUDE_DIR
  localsolver.h PATHS ${LOCALSOLVER_DIR}/include)
find_library(LOCALSOLVER_LIBRARY
  ${LOCALSOLVER_LIB_NAME} PATHS ${LOCALSOLVER_DIR}/bin)

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set LOCALSOLVER_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(LocalSolver DEFAULT_MSG
  LOCALSOLVER_LIBRARY LOCALSOLVER_INCLUDE_DIR)

mark_as_advanced(LOCALSOLVER_LIBRARY LOCALSOLVER_INCLUDE_DIR)

if (LOCALSOLVER_FOUND AND NOT TARGET localsolver-library)
  add_library(localsolver-library STATIC IMPORTED GLOBAL)
  set_target_properties(localsolver-library PROPERTIES
    IMPORTED_LOCATION "${LOCALSOLVER_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LOCALSOLVER_INCLUDE_DIR}")
endif ()

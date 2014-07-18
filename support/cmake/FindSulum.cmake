# Try to find the Sulum libraries.
#
# Once done this will define
#
#  SULUM_FOUND - System has Sulum
#  SULUM_INCLUDE_DIRS - The Sulum include directories
#  SULUM_LIBRARIES - The libraries needed to use Sulum

if (UNIX)
  set(SULUM_SYS linux)
  set(SULUM_DIR /opt/sulum)
else ()
  set(SULUM_SYS win)
  set(PROGRAM_FILES_DIR "C:/Program Files")
  set(PROGRAM_FILES_X86_DIR "${PROGRAM_FILES_DIR} (x86)")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4 AND EXISTS ${PROGRAM_FILES_X86_DIR})
    set(PROGRAM_FILES_DIR ${PROGRAM_FILES_X86_DIR})
  endif ()
  set(SULUM_DIR ${PROGRAM_FILES_DIR}/sulum)
endif ()

file(GLOB SULUM_DIRS "${SULUM_DIR}/*")
if (SULUM_DIRS)
  list(GET SULUM_DIRS 0 SULUM_DIR)
  message(STATUS "Found Sulum directory: ${SULUM_DIR}")
endif ()

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(SULUM_BITS 64)
else ()
  set(SULUM_BITS 32)
endif ()

set(SULUM_BIN_DIR ${SULUM_DIR}/${SULUM_SYS}${SULUM_BITS}/bin)

find_path(SULUM_INCLUDE_DIR sulumcpp.h PATHS ${SULUM_DIR}/header)
find_path(SULUM_AMPL_INCLUDE_DIR optsulum.ampl PATHS ${SULUM_BIN_DIR})
find_library(SULUM_LIBRARY sulum20 PATHS ${SULUM_BIN_DIR})

set(SULUM_INCLUDE_DIRS ${SULUM_INCLUDE_DIR} ${SULUM_AMPL_INCLUDE_DIR})
set(SULUM_LIBRARIES ${SULUM_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set SULUM_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(Sulum DEFAULT_MSG
  SULUM_LIBRARY SULUM_INCLUDE_DIR SULUM_AMPL_INCLUDE_DIR)

mark_as_advanced(SULUM_LIBRARY SULUM_INCLUDE_DIR SULUM_AMPL_INCLUDE_DIR)

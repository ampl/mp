# Try to find the Sulum libraries.
#
# Once done this will define
#
#  Sulum_FOUND - System has Sulum
#  Sulum_INCLUDE_DIRS - The Sulum include directories
#  Sulum_LIBRARIES - The libraries needed to use Sulum

if (UNIX)
  set(Sulum_DIR /opt/sulum)
  set(Sulum_SYS linux)
else ()
  set(Sulum_DIR "C:/Program Files/sulum")
  set(Sulum_SYS win)
endif ()

file(GLOB Sulum_DIRS "${Sulum_DIR}/*")
if (Sulum_DIRS)
  list(GET Sulum_DIRS 0 Sulum_DIR)
  message(STATUS "Found Sulum directory: ${Sulum_DIR}")
endif ()

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(Sulum_BITS 64)
else ()
  set(Sulum_BITS 32)
endif ()

set(Sulum_BIN_DIR ${Sulum_DIR}/${Sulum_SYS}${Sulum_BITS}/bin)

find_path(Sulum_INCLUDE_DIR sulumcpp.h PATHS ${Sulum_DIR}/header)
find_path(Sulum_AMPL_INCLUDE_DIR optsulum.ampl PATHS ${Sulum_BIN_DIR})
find_library(Sulum_LIBRARY sulum20 PATHS ${Sulum_BIN_DIR})

set(Sulum_INCLUDE_DIRS ${Sulum_INCLUDE_DIR} ${Sulum_AMPL_INCLUDE_DIR})
set(Sulum_LIBRARIES ${Sulum_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set LocalSolver_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(SULUM DEFAULT_MSG
  Sulum_LIBRARY Sulum_INCLUDE_DIR Sulum_AMPL_INCLUDE_DIR)

mark_as_advanced(Sulum_LIBRARY Sulum_INCLUDE_DIR Sulum_AMPL_INCLUDE_DIR)

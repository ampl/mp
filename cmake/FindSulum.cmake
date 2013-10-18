# Try to find the Sulum libraries.
#
# Once done this will define
#
#  SULUM_FOUND - System has Sulum
#  SULUM_INCLUDE_DIRS - The Sulum include directories
#  SULUM_LIBRARIES - The libraries needed to use Sulum

if (WIN32)
  set(SULUM_DIR "C:/Program Files/sulum")
else ()
  set(SULUM_DIR /opt/sulum)
endif ()

file(GLOB SULUM_DIRS "${SULUM_DIR}/*")
if (SULUM_DIRS)
  list(GET SULUM_DIRS 0 SULUM_DIR)
  message(STATUS "Found CPLEX Studio: ${SULUM_DIR}")
endif ()

find_path(SULUM_INCLUDE_DIR sulumcpp.h PATHS ${SULUM_DIR}/header)
find_library(SULUM_LIBRARY sulum20 PATHS ${SULUM_DIR}/linux64/bin)

set(SULUM_INCLUDE_DIRS ${SULUM_INCLUDE_DIR})
set(SULUM_LIBRARIES ${SULUM_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set ODBC_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(SULUM DEFAULT_MSG SULUM_LIBRARY SULUM_INCLUDE_DIR)

mark_as_advanced(SULUM_LIBRARY SULUM_INCLUDE_DIR)

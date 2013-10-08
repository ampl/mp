# Try to find the Sulum libraries.
#
# Once done this will define
#
#  SULUM_FOUND - System has Sulum
#  SULUM_INCLUDE_DIRS - The Sulum include directories
#  SULUM_LIBRARIES - The libraries needed to use Sulum

set(SULUM_PATH /opt/sulum/RedEric)
find_path(SULUM_INCLUDE_DIR sulumcpp.h PATHS ${SULUM_PATH}/header)
find_library(SULUM_LIBRARY sulum20 PATHS ${SULUM_PATH}/linux64/bin)

set(SULUM_INCLUDE_DIRS ${SULUM_INCLUDE_DIR})
set(SULUM_LIBRARIES ${SULUM_LIBRARY})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set ODBC_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(SULUM DEFAULT_MSG SULUM_LIBRARY SULUM_INCLUDE_DIR)

mark_as_advanced(SULUM_LIBRARY SULUM_INCLUDE_DIR)

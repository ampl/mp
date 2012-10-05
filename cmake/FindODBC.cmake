# Try to find the ODBC libraries.
#
# Once done this will define
#
#  ODBC_FOUND - System has ODBC
#  ODBC_INCLUDE_DIRS - The ODBC include directories
#  ODBC_LIBRARIES - The libraries needed to use ODBC

find_path(ODBC_INCLUDE_DIR sql.h)
find_library(ODBC_LIBRARY odbc)

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set ODBC_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  ODBC DEFAULT_MSG ODBC_LIBRARY ODBC_INCLUDE_DIR)

set(ODBC_INCLUDE_DIRS ${ODBC_INCLUDE_DIR})
set(ODBC_LIBRARIES ${ODBC_LIBRARY})

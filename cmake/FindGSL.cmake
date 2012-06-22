# Try to find GNU Scientific Library (GSL).
#
# Once done this will define
#
#  GSL_FOUND - System has GSL
#  GSL_INCLUDE_DIRS - The GSL include directories
#  GSL_LIBRARIES - The libraries needed to use GSL

find_library(GSL_LIBRARY NAMES gsl)

set(GSL_LIBRARIES ${GSL_LIBRARY} gslcblas)
set(GSL_INCLUDE_DIRS "")

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARY)

mark_as_advanced(GSL_LIBRARY)

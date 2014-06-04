# Try to find the Gecode libraries.
#
# Once done this will define
#
#  GECODE_FOUND - System has Gecode
#  GECODE_INCLUDE_DIRS - The Gecode include directories
#  GECODE_LIBRARIES - The libraries needed to use Gecode

find_path(GECODE_INCLUDE_DIR gecode/search.hh)
set(GECODE_INCLUDE_DIRS ${Gecode_INCLUDE_DIR})

foreach (module flatzinc driver gist search minimodel set int kernel support)
  set(var GECODE_${module}_LIBRARY)
  find_library(${var} gecode${module})
  set(GECODE_LIBRARIES ${GECODE_LIBRARIES} ${${var}})
  set(GECODE_LIBRARY_VARS ${GECODE_LIBRARY_VARS} ${var})
endforeach ()

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set GECODE_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(
  Gecode DEFAULT_MSG ${GECODE_LIBRARY_VARS} GECODE_INCLUDE_DIR)

mark_as_advanced(${GECODE_LIBRARY_VARS} GECODE_INCLUDE_DIR)

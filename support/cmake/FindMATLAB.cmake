# Try to find MATLAB.
#
# Once done this will define the following variable:
#
#  MATLAB_FOUND - found MATLAB
#  MATLAB_MEX   - path to the MATLAB mex executable

if (UNIX)
  set(MATLAB_DIR /opt/MATLAB)
else ()
  set(PROGRAM_FILES_DIR "C:/Program Files")
  set(PROGRAM_FILES_X86_DIR "${PROGRAM_FILES_DIR} (x86)")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4 AND EXISTS ${PROGRAM_FILES_X86_DIR})
    set(PROGRAM_FILES_DIR ${PROGRAM_FILES_X86_DIR})
  endif ()
  set(MATLAB_DIR "${PROGRAM_FILES_DIR}/MATLAB")
endif ()

file(GLOB MATLAB_DIRS "${MATLAB_DIR}/*")

find_program(MATLAB_MEX mex PATHS ${MATLAB_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set MATLAB_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(MATLAB DEFAULT_MSG MATLAB_MEX)

mark_as_advanced(MATLAB_MEX)

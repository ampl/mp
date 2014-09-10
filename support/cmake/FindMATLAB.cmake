# Try to find MATLAB.
#
# Once done this will define the following variable:
#
#  MATLAB_FOUND - found MATLAB
#  MATLAB_MEX   - path to the MATLAB mex executable

if (APPLE)
  set(MATLAB_GLOB "/Applications/MATLAB*")
  set(MATLAB_MEX_SUFFIX mac)
elseif (UNIX)
  set(MATLAB_GLOB "/opt/MATLAB/*")
  set(MATLAB_MEX_SUFFIX a64)
else ()
  set(PROGRAM_FILES_DIR "C:/Program Files")
  set(PROGRAM_FILES_X86_DIR "${PROGRAM_FILES_DIR} (x86)")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(MATLAB_MEX_SUFFIX w32)
    if (EXISTS ${PROGRAM_FILES_X86_DIR})
      set(PROGRAM_FILES_DIR ${PROGRAM_FILES_X86_DIR})
    endif ()
  else ()
    set(MATLAB_MEX_SUFFIX w64)
  endif ()
  set(MATLAB_GLOB "${PROGRAM_FILES_DIR}/MATLAB/*")
  set(MEX_EXT .bat)
endif ()

file(GLOB MATLAB_DIRS "${MATLAB_GLOB}")

find_program(MATLAB_MEX mex${MEX_EXT}
  PATHS ${MATLAB_DIRS} PATH_SUFFIXES bin NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set MATLAB_FOUND to TRUE
# if all listed variables are TRUE.
find_package_handle_standard_args(MATLAB DEFAULT_MSG MATLAB_MEX)

mark_as_advanced(MATLAB_MEX)

function (add_mex name)
  cmake_parse_arguments(add_mex "" "" "COMPILE_FLAGS;LIBRARIES" ${ARGN})
  set(filename ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name}.mex${MATLAB_MEX_SUFFIX})
  set(sources ${add_mex_UNPARSED_ARGUMENTS})
  set(libs)
  foreach (lib ${add_mex_LIBRARIES})
    set(libs ${libs} $<TARGET_FILE:${lib}>)
  endforeach ()
  add_custom_command(OUTPUT ${filename}
    COMMAND ${MATLAB_MEX} ${add_mex_COMPILE_FLAGS}
      ${sources} ${libs} -output ${filename}
    DEPENDS ${sources} ${add_mex_LIBRARIES})
  add_custom_target(${name} ALL SOURCES ${filename})
endfunction ()

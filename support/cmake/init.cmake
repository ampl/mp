# CMake initialization code that should be run before the project command.

include(CMakeParseArguments)

# Joins arguments and sets the result to <var>.
# Usage:
#   join(<var> [<arg>...])
function (join var)
  unset(result)
  foreach (arg ${ARGN})
    if (DEFINED result)
      set(result "${result} ${arg}")
    else ()
      set(result "${arg}")
    endif ()
  endforeach ()
  set(${var} "${result}" PARENT_SCOPE)
endfunction ()

# Sets cache variable <var> to the value <value>. The arguments
# following <type> are joined into a single docstring which allows
# breaking long documentation into smaller strings.
# Usage:
#   set_cache(<var> <value> <type> docstring... [FORCE])
function (set_cache var value type)
  cmake_parse_arguments(set_cache FORCE "" "" ${ARGN})
  unset(force)
  if (set_cache_FORCE)
    set(force FORCE)
  endif ()
  join(docstring ${set_cache_UNPARSED_ARGUMENTS})
  set(${var} ${value} CACHE ${type} "${docstring}" ${force})
endfunction ()

if (NOT CMAKE_BUILD_TYPE)
  # Set the default CMAKE_BUILD_TYPE to Release.
  # This should be done before the project command since the latter sets
  # CMAKE_BUILD_TYPE itself.
  set_cache(CMAKE_BUILD_TYPE Release STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or"
    "CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif ()

function (override var file)
  if (EXISTS "${file}")
    message("!!! override ${var} ${file}")
    set(${var} ${file} PARENT_SCOPE)
  endif ()
endfunction ()

# Use static MSVC runtime.
# This should be done before the project command.
override(CMAKE_USER_MAKE_RULES_OVERRIDE
  ${CMAKE_CURRENT_LIST_DIR}/c_flag_overrides.cmake)
override(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
  ${CMAKE_CURRENT_LIST_DIR}/cxx_flag_overrides.cmake)

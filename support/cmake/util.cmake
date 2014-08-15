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

# This CMake script generates the config header ${CONFIG_HEADER}.

# Gets the date of the last revision.
execute_process(
  COMMAND git log -1 --date=iso --pretty=format:%cd
  OUTPUT_VARIABLE DATE)
if (DATE MATCHES "([0-9]+)-([0-9]+)-([0-9]+).*")
  set(DATE ${CMAKE_MATCH_1}${CMAKE_MATCH_2}${CMAKE_MATCH_3})
else ()
  set(DATE 0)
endif ()

set(CONFIG
"#define MP_DATE ${DATE}
#define MP_SYSINFO \"${SYSINFO}\"
")

if (EXISTS ${CONFIG_HEADER})
  file(READ ${CONFIG_HEADER} OLD_CONFIG)
endif ()
if (NOT CONFIG STREQUAL OLD_CONFIG)
  message(STATUS "Writing ${CONFIG_HEADER}")
  file(WRITE ${CONFIG_HEADER} "${CONFIG}")
endif ()

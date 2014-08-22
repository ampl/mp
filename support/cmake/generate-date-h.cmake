# This CMake script generates the header ${DATE_HEADER} with the
# date of the last revision.

execute_process(
  COMMAND git log -1 --date=iso --pretty=format:%cd
  OUTPUT_VARIABLE DATE)
if (DATE MATCHES "([0-9]+)-([0-9]+)-([0-9]+).*")
  set(DATE ${CMAKE_MATCH_1}${CMAKE_MATCH_2}${CMAKE_MATCH_3})
else ()
  set(DATE 0)
endif ()

set(CONTENT "#define MP_DATE ${DATE}\n")

if (EXISTS ${DATE_HEADER})
  file(READ ${DATE_HEADER} OLD_CONTENT)
endif ()
if (NOT CONTENT STREQUAL OLD_CONTENT)
  message(STATUS "Writing ${DATE_HEADER}")
  file(WRITE ${DATE_HEADER} "${CONTENT}")
endif ()

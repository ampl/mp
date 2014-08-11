include_directories(../thirdparty)

add_test(support-test nose2 -s support)

# Add gtest library with some extra features.
find_package(Threads)
add_library(gtest-extra
  ../thirdparty/gtest/gtest-all.cc util.h util.cc test-main.cc)
target_link_libraries(gtest-extra mp ${CMAKE_THREAD_LIBS_INIT})

# Adds an MP test.
# Usage:
#   add_mp_test(name [PREFIX] sources... [LIBS libraries])
function(add_mp_test name)
  set(srcs )
  set(libs )
  set(library FALSE)
  foreach (arg ${ARGN})
    if (library)
      set (libs ${libs} ${arg})
    elseif (arg STREQUAL LIBS)
      set (library TRUE)
    elseif (arg STREQUAL PREFIX)
      set(prefix mp)
    else ()
      set (srcs ${srcs} ${arg})
    endif ()
  endforeach ()
  set(target_name "${prefix}${name}")
  add_executable(${target_name} ${srcs})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} gtest-extra ${libs})
  add_test(NAME ${name} COMMAND $<TARGET_FILE:${target_name}>)
endfunction()

add_mp_test(clock-test clock-test.cc)
add_mp_test(error-test error-test.cc)
add_mp_test(nl-test nl-test.cc)

add_executable(test-helper test-helper.cc)
target_link_libraries(test-helper mp)

add_mp_test(os-test os-test.cc)
add_dependencies(os-test test-helper)

add_mp_test(rstparser-test rstparser-test.cc)
add_mp_test(safeint-test safeint-test.cc)

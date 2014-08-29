include_directories(../thirdparty)

add_test(cmake-test
  ${CMAKE_COMMAND} -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/init-test.cmake)

find_program(NOSE2 nose2)
if (NOSE2)
  add_test(support-test ${NOSE2} -s ${CMAKE_CURRENT_SOURCE_DIR}/support)
endif ()

# Add gtest library with some extra features.
set(GTEST_SRC ../thirdparty/gtest/gtest-all.cc util.h util.cc test-main.cc)
add_library(gtest-extra ${GTEST_SRC})
find_package(Threads)
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
  target_compile_definitions(${target_name} PRIVATE
    MP_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/test/data")
  add_test(NAME ${name} COMMAND $<TARGET_FILE:${target_name}>)
endfunction()

add_ampl_library(testlib testlib.cc)
add_library(function function.cc function.h)
target_link_libraries(function asl)

add_mp_test(function-test function-test.cc LIBS function)

if (TARGET amplgsl)
  include_directories(${GSL_BINARY_DIR})
  add_mp_test(gsl-test gsl-test.cc LIBS function gsl gslcblas)
endif ()

add_mp_test(cp-test cp-test.cc LIBS function)
add_dependencies(cp-test cp)

add_subdirectory(asl)
add_subdirectory(solvers)
add_subdirectory(tables)

add_mp_test(util-test util-test.cc)

add_mp_test(clock-test clock-test.cc)
add_mp_test(error-test error-test.cc)
add_mp_test(nl-test nl-test.cc)

add_executable(test-helper test-helper.cc)
target_link_libraries(test-helper mp)

add_mp_test(os-test os-test.cc)
add_dependencies(os-test test-helper)

add_mp_test(rstparser-test rstparser-test.cc)
add_mp_test(safeint-test safeint-test.cc)

add_library(ampltestsolver SHARED testsolver.cc ../../src/solver-c.cc)
target_link_libraries(ampltestsolver asl)

add_mp_test(solver-c-test solver-c-test.cc LIBS ampltestsolver)

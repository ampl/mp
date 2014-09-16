include_directories(../thirdparty)

add_test(cmake-test
  ${CMAKE_COMMAND} -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/init-test.cmake)

find_program(NOSE2 nose2)
if (NOSE2)
  add_test(support-test ${NOSE2} -s ${CMAKE_CURRENT_SOURCE_DIR}/support)
endif ()

# Add gtest library with some extra features.
set(GTEST_SRC ../thirdparty/gtest/gtest-all.cc
  gtest-extra.cc gtest-extra.h util.cc util.h test-main.cc)
add_library(gtest-extra STATIC ${GTEST_SRC})
target_compile_definitions(gtest-extra PUBLIC FMT_USE_FILE_DESCRIPTORS)
find_package(Threads)
target_link_libraries(gtest-extra ${CMAKE_THREAD_LIBS_INIT})

# Adds an MP test.
# Usage:
#   add_mp_test(name [PREFIX] sources... [LIBS libraries])
function(add_mp_test name)
  cmake_parse_arguments(add_mp_test PREFIX "" LIBS ${ARGN})
  set(target_name ${name})
  if (add_mp_test_PREFIX)
    set(target_name mp${name})
  endif ()
  add_executable(${target_name} ${add_mp_test_UNPARSED_ARGUMENTS})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  if (MINGW)
    set_target_properties(${target_name} PROPERTIES
      LINK_FLAGS "-static-libgcc -static-libstdc++")
  endif ()
  if (NOT add_mp_test_LIBS)
    set(add_mp_test_LIBS mp)
  endif ()
  target_link_libraries(${target_name} gtest-extra ${add_mp_test_LIBS})
  target_compile_definitions(${target_name} PRIVATE
    MP_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/test/data")
  add_test(NAME ${name} COMMAND $<TARGET_FILE:${target_name}>)
endfunction()

add_ampl_library(testlib PRIVATE testlib.cc)
add_library(function STATIC function.cc function.h)
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
if (MINGW)
  set_target_properties(test-helper PROPERTIES
    LINK_FLAGS "-static-libgcc -static-libstdc++")
endif ()

add_mp_test(os-test os-test.cc)
add_dependencies(os-test test-helper)

find_program(LSOF lsof)
if (LSOF)
  target_compile_definitions(os-test PRIVATE HAVE_LSOF=1)
endif ()

add_mp_test(problem-builder-test problem-builder-test.cc)

add_mp_test(rstparser-test rstparser-test.cc)
add_mp_test(safeint-test safeint-test.cc)

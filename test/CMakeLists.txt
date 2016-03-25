include_directories(../thirdparty)

if (MSVC)
  # Fix error C1128: number of sections exceeded object file format limit.
  add_compile_options(/bigobj)
endif ()

add_test(cmake-test
  ${CMAKE_COMMAND} -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/init-test.cmake)

find_program(NOSE2 nose2)
if (NOSE2)
  add_test(support-test ${NOSE2} -s ${CMAKE_CURRENT_SOURCE_DIR}/support)
endif ()

add_prefix(GTEST_SRC ../thirdparty/gmock/ gmock-gtest-all.cc gmock.h)

# Add gtest library with some extra features.
set(GTEST_SRC ${GTEST_SRC}
              gtest-extra.cc gtest-extra.h util.cc util.h test-main.cc)
add_library(gtest-extra STATIC ${GTEST_SRC})
target_compile_definitions(gtest-extra PUBLIC FMT_USE_FILE_DESCRIPTORS)
find_package(Threads)
target_link_libraries(gtest-extra ${CMAKE_THREAD_LIBS_INIT})

# Workaround GTest bug https://github.com/google/googletest/issues/705.
check_cxx_compiler_flag(
  -fno-delete-null-pointer-checks HAVE_FNO_DELETE_NULL_POINTER_CHECKS)
if (HAVE_FNO_DELETE_NULL_POINTER_CHECKS)
  target_compile_options(gtest-extra PUBLIC -fno-delete-null-pointer-checks)
endif ()

# Workaround a bug in implementation of variadic templates in MSVC11.
if (MSVC)
  target_compile_definitions(gtest-extra PUBLIC _VARIADIC_MAX=10)
endif ()

check_cxx_source_compiles("
  #include <initializer_list>
  int main() {}" MP_INITIALIZER_LISTS)

if (NOT MP_VARIADIC_TEMPLATES OR NOT MP_INITIALIZER_LISTS)
  add_definitions(-DGTEST_LANG_CXX11=0)
endif ()

# Adds an MP test.
# Usage:
#   add_mp_test(name sources... [LIBS libraries])
function(add_mp_test name)
  cmake_parse_arguments(add_mp_test "" "" LIBS ${ARGN})
  add_executable(${name} ${add_mp_test_UNPARSED_ARGUMENTS})
  set_target_properties(${name} PROPERTIES OUTPUT_NAME ${name})
  if (MINGW)
    set_target_properties(${name} PROPERTIES
      LINK_FLAGS "-static-libgcc -static-libstdc++")
  endif ()
  if (NOT add_mp_test_LIBS)
    set(add_mp_test_LIBS mp)
  endif ()
  target_link_libraries(${name} gtest-extra ${add_mp_test_LIBS})
  target_compile_definitions(${name} PRIVATE
    MP_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/test/data")
  add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}>)
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

add_executable(test-helper test-helper.cc)
target_link_libraries(test-helper mp)
if (MINGW)
  set_target_properties(test-helper PROPERTIES
    LINK_FLAGS "-static-libgcc -static-libstdc++")
endif ()

add_mp_test(assert-test assert-test.cc)
add_mp_test(clock-test clock-test.cc)
add_mp_test(common-test common-test.cc)
add_mp_test(error-test error-test.cc)
add_mp_test(expr-test expr-test.cc mock-allocator.h test-assert.h)
add_mp_test(expr-visitor-test expr-visitor-test.cc test-assert.h)
add_mp_test(expr-writer-test expr-writer-test.cc)
add_mp_test(nl-reader-test nl-reader-test.cc mock-file.h mock-problem-builder.h)
add_mp_test(option-test option-test.cc)
add_mp_test(os-test os-test.cc mock-file.h)
add_dependencies(os-test test-helper)
add_mp_test(problem-test problem-test.cc)

add_mp_test(solver-test
  solver-test.cc mock-problem-builder.h solution-handler.h)
target_compile_definitions(solver-test
  PRIVATE MP_SYSINFO="${MP_SYSINFO}" MP_DATE=${MP_DATE})

add_mp_test(suffix-test suffix-test.cc)

find_program(LSOF lsof)
if (LSOF)
  target_compile_definitions(os-test PRIVATE HAVE_LSOF=1)
endif ()

add_mp_test(problem-builder-test problem-builder-test.cc)

add_mp_test(rstparser-test rstparser-test.cc)
add_mp_test(safeint-test safeint-test.cc)

if (TARGET benchmark)
  add_executable(nl-microbenchmark nl-benchmark.cc)
  target_link_libraries(nl-microbenchmark benchmark)
endif ()

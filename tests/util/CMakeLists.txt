cxx_test_with_flags(clock-test "" "gtest;util" ../test-main.cc clock-test.cc)

cxx_test_with_flags(error-test "" "gtest;util" ../test-main.cc error-test.cc)

cxx_test_with_flags(expr-test "" "gtest;util"
  ../test-main.cc ../expr-builder.cc expr-test.cc no-unordered-map-test.cc)

cxx_test_with_flags(expr-factory-test "" "gtest;util"
  ../test-main.cc ../util.cc expr-factory-test.cc)

cxx_test_with_flags(format-test "" "gtest;util" format-test.cc)

cxx_test_with_flags(nl-test "" "gtest;util"
  ../test-main.cc ../util.cc nl-test.cc)

add_executable(numberofmap-speed-test
  numberofmap-speed-test.cc ../expr-builder.cc)
target_link_libraries(numberofmap-speed-test util)

add_executable(test-helper test-helper.cc)
target_link_libraries(test-helper util)

cxx_test_with_flags(os-test "" "gtest;util"
  ../test-main.cc ../util.cc os-test.cc)
add_dependencies(os-test test-helper)

cxx_test_with_flags(problem-test "" "gtest;util"
  ../test-main.cc ../expr-builder.cc ../util.cc problem-test.cc)

cxx_test_with_flags(rstparser-test "" "gtest;util" rstparser-test.cc)

cxx_test_with_flags(solver-test "" "gtest;util"
  ../test-main.cc ../args.cc ../util.cc solver-test.cc)

add_library(ampltestsolver SHARED testsolver.cc ../../solvers/util/solver-c.cc)
target_link_libraries(ampltestsolver util)

cxx_test_with_flags(solver-c-test "" "gtest;ampltestsolver;util"
  ../test-main.cc solver-c-test.cc)

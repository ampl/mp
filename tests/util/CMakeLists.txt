add_definitions(-DFMT_USE_FILE_DESCRIPTORS=1)

add_ampl_test(clock-test clock-test.cc)
add_ampl_test(error-test error-test.cc)

if (HAVE_AUTO)
  add_ampl_test(expr-test
    ../expr-builder.cc expr-test.cc no-unordered-map-test.cc)
endif ()

add_ampl_test(aslbuilder-test aslbuilder-test.cc)
add_ampl_test(format-test format-test.cc gtest-extra.cc util.h util.cc)
add_ampl_test(gtest-extra-test gtest-extra-test.cc gtest-extra.cc)
add_ampl_test(nl-test nl-test.cc)

add_executable(numberofmap-speed-test
  numberofmap-speed-test.cc ../expr-builder.cc)
target_link_libraries(numberofmap-speed-test util)

add_executable(test-helper test-helper.cc)
target_link_libraries(test-helper util)

add_ampl_test(os-test os-test.cc)
add_dependencies(os-test test-helper)

add_ampl_test(problem-test ../expr-builder.cc problem-test.cc)
add_ampl_test(rstparser-test rstparser-test.cc)
add_ampl_test(solver-test ../args.cc solver-test.cc)

add_library(ampltestsolver SHARED testsolver.cc ../../solvers/util/solver-c.cc)
target_link_libraries(ampltestsolver util)

add_ampl_test(solver-c-test solver-c-test.cc LIBS ampltestsolver)

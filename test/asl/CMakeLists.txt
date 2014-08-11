if (HAVE_AUTO)
  add_ampl_test(expr-test expr-test.cc no-unordered-map-test.cc)
endif ()

add_ampl_test(aslbuilder-test aslbuilder-test.cc)

add_executable(numberofmap-speed-test numberofmap-speed-test.cc)
target_link_libraries(numberofmap-speed-test util)

add_ampl_test(problem-test problem-test.cc)
add_ampl_test(solver-test ../args.cc solver-test.cc)

add_library(ampltestsolver SHARED testsolver.cc ../../solvers/util/solver-c.cc)
target_link_libraries(ampltestsolver util)

add_ampl_test(solver-c-test solver-c-test.cc LIBS ampltestsolver)

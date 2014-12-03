add_mp_test(aslexpr-test aslexpr-test.cc no-unordered-map-test.cc)
target_link_libraries(aslexpr-test asl)

add_mp_test(aslbuilder-test aslbuilder-test.cc LIBS asl)

add_executable(numberofmap-speed-test numberofmap-speed-test.cc)
target_link_libraries(numberofmap-speed-test asl)

add_mp_test(aslproblem-test aslproblem-test.cc stderr-redirect.h LIBS asl)
add_mp_test(aslsolver-test aslsolver-test.cc LIBS asl)

add_library(ampltestsolver SHARED
  testsolver.cc ${PROJECT_SOURCE_DIR}/src/solver-c.cc)
target_link_libraries(ampltestsolver asl)

add_mp_test(solver-c-test solver-c-test.cc LIBS ampltestsolver)

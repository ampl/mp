if (HAVE_AUTO)
  add_mp_test(expr-test expr-test.cc no-unordered-map-test.cc)
endif ()

add_mp_test(aslbuilder-test aslbuilder-test.cc LIBS asl)

add_executable(numberofmap-speed-test numberofmap-speed-test.cc)
target_link_libraries(numberofmap-speed-test asl)

add_mp_test(problem-test problem-test.cc stderr-redirect.h LIBS asl)
add_mp_test(solver-test PREFIX ../args.cc solver-test.cc LIBS asl)
target_compile_definitions(mpsolver-test PRIVATE
  MP_SYSINFO="${MP_SYSINFO}" MP_DATE=${MP_DATE})

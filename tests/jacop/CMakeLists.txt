include_directories(${JNI_INCLUDE_DIRS})
add_custom_command(OUTPUT ${JACOP_JAR}
  COMMAND ${CMAKE_COMMAND} -E copy
    ../../solvers/jacop/lib/${JACOP_JAR} ${JACOP_JAR})
if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
endif ()
cxx_test_with_flags(jacop_test "" "gtest;ampljacop" ${JACOP_JAR}
  ../test_main.cc ../args.cc ../expr_builder.cc
  ../solver_test.cc ../util.cc jacop_test.cc)

include_directories(${JNI_INCLUDE_DIRS})
add_custom_command(OUTPUT ${JACOP_JAR}
  COMMAND ${CMAKE_COMMAND} -E copy
    ../../solvers/jacop/lib/${JACOP_JAR} ${JACOP_JAR})
add_custom_command(OUTPUT jacoplisteners.jar
  COMMAND ${CMAKE_COMMAND} -E copy
    ../../solvers/jacop/jacoplisteners.jar jacoplisteners.jar)
cxx_test_with_flags(jacop_test "" "gtest;ampljacop" ${JACOP_JAR} jacoplisteners.jar
  ../test_main.cc ../args.cc ../expr_builder.cc
  ../solver_test.cc ../util.cc jacop_test.cc)
if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
  target_link_libraries(jacop_test delayimp)
endif ()

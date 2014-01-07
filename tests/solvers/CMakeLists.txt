cxx_test_with_flags(cp-test "" "gtest;amplsolver;function"
  ../test-main.cc cp-test.cc)
add_dependencies(cp-test cp)

if (HAVE_GECODE)
  include_directories(${GECODE_SOURCE_DIR})
  cxx_test_with_flags(amplgecode-test "" "gtest;amplgecode-static"
    ../test-main.cc ../args.cc ../expr-builder.cc
    ../solver-test.cc ../util.cc gecode-test.cc)
  set_target_properties(amplgecode-test PROPERTIES OUTPUT_NAME gecode-test)
endif ()

if (HAVE_JACOP)
  include_directories(${JNI_INCLUDE_DIRS})
  add_custom_command(OUTPUT ${JACOP_JAR}
    COMMAND ${CMAKE_COMMAND} -E copy
      ../../solvers/jacop/lib/${JACOP_JAR} ${JACOP_JAR})
  add_custom_command(OUTPUT ampljacop.jar DEPENDS ampljacop
    COMMAND ${CMAKE_COMMAND} -E copy
      ../../solvers/jacop/ampljacop.jar ampljacop.jar)
  cxx_test_with_flags(jacop-test "" "gtest;ampljacop-static" ${JACOP_JAR} ampljacop.jar
    ../test-main.cc ../args.cc ../expr-builder.cc
    ../solver-test.cc ../util.cc jacop-test.cc)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

cxx_test_with_flags(smpswriter-test "" "gtest;amplsmpswriter"
  ../test-main.cc ../args.cc ../util.cc smpswriter-test.cc)

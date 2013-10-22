cxx_test_with_flags(cp_test "" "gtest;amplsolver;function"
  ../test_main.cc cp_test.cc)
add_dependencies(cp_test cp)

if (HAVE_GECODE)
  include_directories(${GECODE_SOURCE_DIR})
  cxx_test_with_flags(gecode_test "" "gtest;amplgecode"
    ../test_main.cc ../args.cc ../expr_builder.cc
    ../solver_test.cc ../util.cc gecode_test.cc)
endif ()

if (HAVE_JACOP)
  include_directories(${JNI_INCLUDE_DIRS})
  add_custom_command(OUTPUT ${JACOP_JAR}
    COMMAND ${CMAKE_COMMAND} -E copy
      ../../solvers/jacop/lib/${JACOP_JAR} ${JACOP_JAR})
  add_custom_command(OUTPUT ampljacop.jar DEPENDS ampljacop_jar
    COMMAND ${CMAKE_COMMAND} -E copy
      ../../solvers/jacop/ampljacop.jar ampljacop.jar)
  cxx_test_with_flags(jacop_test "" "gtest;ampljacop" ${JACOP_JAR} ampljacop.jar
    ../test_main.cc ../args.cc ../expr_builder.cc
    ../solver_test.cc ../util.cc jacop_test.cc)
  if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop_test delayimp)
  endif ()
endif ()

cxx_test_with_flags(smpswriter_test "" "gtest;amplsmpswriter"
  ../test_main.cc ../args.cc ../util.cc smpswriter_test.cc)

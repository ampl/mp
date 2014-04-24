add_ampl_test(cp-test cp-test.cc LIBS amplsolver function)
add_dependencies(cp-test cp)

if (HAVE_GECODE)
  include_directories(${GECODE_SOURCE_DIR})
  add_ampl_test(gecode-test PREFIX
    ../args.cc ../expr-builder.cc ../solver-test.cc gecode-test.cc
    LIBS amplgecode-static)
endif ()

if (HAVE_JACOP)
  include_directories(${JNI_INCLUDE_DIRS})
  add_custom_command(OUTPUT ${JACOP_JAR}
    COMMAND ${CMAKE_COMMAND} -E copy
      ../../solvers/jacop/lib/${JACOP_JAR} ${JACOP_JAR})
  add_custom_command(OUTPUT ampljacop.jar DEPENDS ampljacop
    COMMAND ${CMAKE_COMMAND} -E copy
      ../../solvers/jacop/ampljacop.jar ampljacop.jar)
  add_ampl_test(jacop-test ${JACOP_JAR} ampljacop.jar
    ../args.cc ../expr-builder.cc ../solver-test.cc jacop-test.cc
    LIBS ampljacop-static)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

add_ampl_test(smpswriter-test ../args.cc smpswriter-test.cc LIBS amplsmpswriter)

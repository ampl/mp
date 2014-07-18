add_ampl_test(cp-test cp-test.cc LIBS amplsolver function)
add_dependencies(cp-test cp)

if (TARGET gecode)
  add_ampl_test(gecode-test PREFIX
    ../args.cc ../expr-builder.cc ../solver-test.cc gecode-test.cc
    LIBS amplgecode-static)
endif ()

if (TARGET ilogcp)
  include_directories(../gtest)
  # Include the gtest source in test instead of linking with the library to
  # avoid a problem when gtest is built with a different runtime than ilogcp.
  add_executable(ilogcp-test
    ../test-main.cc ../args.cc ../expr-builder.cc ../solver-test.cc
    ../gtest/src/gtest-all.cc ilogcp-test.cc)
  target_link_libraries(ilogcp-test amplilogcp-static)
  add_test(NAME ilogcp-test COMMAND $<TARGET_FILE:ilogcp-test>)
endif ()

if (TARGET jacop)
  add_ampl_test(jacop-test
    ../args.cc ../expr-builder.cc ../solver-test.cc jacop-test.cc
    LIBS ampljacop-static)
  add_custom_command(TARGET jacop-test
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ../../solvers/jacop/lib/${JACOP_JAR} $<TARGET_FILE_DIR:jacop-test>)
  add_custom_command(TARGET jacop-test DEPENDS ampljacop
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ../../solvers/jacop/ampljacop.jar $<TARGET_FILE_DIR:jacop-test>)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

if (TARGET amplsmpswriter)
  add_ampl_test(smpswriter-test ../args.cc smpswriter-test.cc LIBS amplsmpswriter)
endif ()

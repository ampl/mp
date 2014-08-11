if (TARGET gecode)
  add_mp_test(gecode-test PREFIX
    ../args.cc ../solver-test.cc gecode-test.cc LIBS amplgecode-static)
endif ()

if (TARGET ilogcp)
  include_directories(../gtest)
  # Include the gtest source in test instead of linking with the library to
  # avoid a problem when gtest is built with a different runtime than ilogcp.
  add_executable(ilogcp-test
    ../test-main.cc ../args.cc ../solver-test.cc
    ../gtest/src/gtest-all.cc ilogcp-test.cc)
  target_link_libraries(ilogcp-test amplilogcp-static)
  add_test(NAME ilogcp-test COMMAND $<TARGET_FILE:ilogcp-test>)
endif ()

if (TARGET jacop)
  add_mp_test(jacop-test
    ../args.cc ../solver-test.cc jacop-test.cc LIBS ampljacop-static)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

if (TARGET amplsmpswriter)
  add_mp_test(smpswriter-test ../args.cc smpswriter-test.cc LIBS amplsmpswriter)
endif ()

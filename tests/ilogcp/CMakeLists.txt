find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  include_directories(../gtest
    ${CPLEX_ILOCPLEX_INCLUDE_DIRS} ${CPLEX_CP_INCLUDE_DIRS})
  add_definitions(${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})

  use_dll_runtime()

  # Include the gtest source in test instead of linking with the library to
  # avoid a problem when gtest is built with a different runtime than ilogcp.
  add_executable(ilogcp-test
    ../test-main.cc ../args.cc ../expr-builder.cc ../solver-test.cc
    ../gtest/src/gtest-all.cc ilogcp-test.cc)
  target_link_libraries(ilogcp-test amplilogcp-static)
  add_test(NAME ilogcp-test COMMAND $<TARGET_FILE:ilogcp-test>)
  copy_cplex_dll(ilogcp-test)
endif ()

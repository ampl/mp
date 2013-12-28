cxx_test_with_flags(basic_solver_test "" "gtest;util"
  ../test_main.cc ../args.cc ../util.cc basic_solver_test.cc)

cxx_test_with_flags(clock_test "" "gtest;util" ../test_main.cc clock_test.cc)

cxx_test_with_flags(error_test "" "gtest;util" ../test_main.cc error_test.cc)

cxx_test_with_flags(expr_test "" "gtest;util"
  ../test_main.cc ../expr_builder.cc expr_test.cc no_unordered_map_test.cc)

cxx_test_with_flags(format_test "" "gtest;util" format_test.cc)

add_executable(numberofmap_speed_test
  numberofmap_speed_test.cc ../expr_builder.cc)
target_link_libraries(numberofmap_speed_test util)

add_executable(print-executable-path print-executable-path.cc)
target_link_libraries(print-executable-path util)

cxx_test_with_flags(os_test "" "gtest;util"
  ../test_main.cc ../util.cc os_test.cc)
add_dependencies(os_test print-executable-path)

cxx_test_with_flags(problem_test "" "gtest;util"
  ../test_main.cc ../expr_builder.cc ../util.cc problem_test.cc)

cxx_test_with_flags(rstparser_test "" "gtest;util" rstparser_test.cc)

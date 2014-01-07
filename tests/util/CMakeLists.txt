cxx_test_with_flags(basic-solver-test "" "gtest;util"
  ../test_main.cc ../args.cc ../util.cc basic-solver-test.cc)

cxx_test_with_flags(clock-test "" "gtest;util" ../test_main.cc clock-test.cc)

cxx_test_with_flags(error-test "" "gtest;util" ../test_main.cc error-test.cc)

cxx_test_with_flags(expr-test "" "gtest;util"
  ../test_main.cc ../expr_builder.cc expr-test.cc no-unordered-map-test.cc)

cxx_test_with_flags(format-test "" "gtest;util" format-test.cc)

add_executable(numberofmap-speed-test
  numberofmap-speed-test.cc ../expr_builder.cc)
target_link_libraries(numberofmap-speed-test util)

add_executable(print-executable-path print-executable-path.cc)
target_link_libraries(print-executable-path util)

cxx_test_with_flags(os-test "" "gtest;util"
  ../test_main.cc ../util.cc os-test.cc)
add_dependencies(os-test print-executable-path)

cxx_test_with_flags(problem-test "" "gtest;util"
  ../test_main.cc ../expr_builder.cc ../util.cc problem-test.cc)

cxx_test_with_flags(rstparser-test "" "gtest;util" rstparser-test.cc)

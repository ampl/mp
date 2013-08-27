cxx_test_with_flags(gecode_test "" "gtest;amplgecode"
  ../test_main.cc ../args.cc ../expr_builder.cc
  ../solver_test.cc ../util.cc gecode_test.cc)

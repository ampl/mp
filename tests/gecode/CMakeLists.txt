# Disable -pedantic when building Gecode because it causes errors
# when compiling on Mac with older GCC.
string(REPLACE "-pedantic" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

include_directories(${GECODE_SOURCE_DIR})
cxx_test_with_flags(gecode_test "" "gtest;amplgecode"
  ../test_main.cc ../args.cc ../expr_builder.cc
  ../solver_test.cc ../util.cc gecode_test.cc)

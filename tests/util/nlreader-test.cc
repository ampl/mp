/*
 An .nl reader tests.

 Copyright (C) 2013 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "gtest/gtest.h"
#include "solvers/util/nlreader.h"
#include "solvers/util/problem.h"
#include "tests/util.h"

using ampl::NLReader;

namespace {

TEST(NLReaderTest, InvalidFormat) {
  EXPECT_THROW_MSG(NLReader().ReadString("x"),
      ampl::ParseError, "(input):1:1: invalid format 'x'");
}

// The problem below represents the following AMPL problem in NL format
// without the first line (options):
//   var x >= 0;
//   minimize o: x;
const char *TEST_PROBLEM_NO_OPTIONS =
  " 1 1 0\n"
  " 0 0\n"
  " 0 0\n"
  " 0 0 0\n"
  " 0 0 0 1\n"
  " 0 0 0 0 0\n"
  " 0 0\n"
  " 0 0\n"
  " 0 0 0 0 0\n";

struct TestNLHandler : ampl::NLHandler {
  ampl::NLHeader header;
  void HandleHeader(const ampl::NLHeader &h) { header = h; }
};

ampl::NLHeader ReadOptions(const char *options) {
  fmt::Writer w;
  w << options << '\n' << TEST_PROBLEM_NO_OPTIONS;
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(w.c_str());
  return handler.header;
}

TEST(NLReaderTest, InvalidNumOptions) {
  EXPECT_EQ(0, ReadOptions("ga").num_options);
  EXPECT_EQ(0, ReadOptions("g-1").num_options);
  EXPECT_THROW_MSG(NLReader().ReadString("g10"),
      ampl::ParseError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(NLReader().ReadString(
      str(fmt::Format("g{}") << static_cast<unsigned>(INT_MAX) + 1)),
      ampl::ParseError, "(input):1:2: number is too big");
}

void CheckReadOptions(size_t num_options,
    size_t num_options_to_write, const int *options) {
  fmt::Writer w;
  w << 'g' << num_options;
  for (size_t i = 0; i < num_options_to_write; ++i)
    w << ' ' << options[i];
  ampl::NLHeader header = ReadOptions(w.c_str());
  ASSERT_EQ(num_options, header.num_options);
  size_t min_num_options = std::min(num_options, num_options_to_write);
  for (size_t i = 0; i < min_num_options; ++i)
    EXPECT_EQ(options[i], header.options[i]);
  for (size_t i = min_num_options; i < num_options_to_write; ++i)
    EXPECT_EQ(0, header.options[i]);
}

TEST(NLReaderTest, ReadOptions) {
  const int options[ampl::MAX_NL_OPTIONS + 1] = {
      3, 5, 7, 11, 13, 17, 19, 23, 29, 31
  };
  for (size_t i = 0; i < ampl::MAX_NL_OPTIONS; ++i) {
    for (size_t j = 0; j < ampl::MAX_NL_OPTIONS + 1; ++j)
      CheckReadOptions(i, j, options);
  }
  EXPECT_EQ(0, ReadOptions("g").num_options);
}

TEST(NLReaderTest, ReadAMPLVBTol) {
  EXPECT_EQ(4.2, ReadOptions("g2 0 3 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadOptions("g2 0 0 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadOptions("g2 0 3").ampl_vbtol);
}

TEST(NLReaderTest, NumVars) {
  EXPECT_EQ(1, ReadOptions("g").num_vars);
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(
    "g\n"
    " 42 1 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0\n"
    " 0 0 0 1\n"
    " 0 0 0 0 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0 0 0\n");
  EXPECT_EQ(42, handler.header.num_vars);
}

TEST(NLReaderTest, NumCons) {
  EXPECT_EQ(1, ReadOptions("g").num_cons);
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(
    "g\n"
    " 1 42 1\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0\n"
    " 0 0 0 1\n"
    " 0 0 0 0 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0 0 0\n");
  EXPECT_EQ(42, handler.header.num_cons);
}

TEST(NLReaderTest, NumObjs) {
  EXPECT_EQ(0, ReadOptions("g").num_objs);
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(
    "g\n"
    " 1 0 42\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0\n"
    " 0 0 0 1\n"
    " 0 0 0 0 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0 0 0\n"
    "O0 0\n"
    "n0\n"
    "k0\n");
  EXPECT_EQ(42, handler.header.num_objs);
}

// TODO: more tests
}

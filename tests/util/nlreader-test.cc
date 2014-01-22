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

TEST(NLReaderTest, InvalidNumOptions) {
  EXPECT_THROW_MSG(NLReader().ReadString("ga"),
      ampl::ParseError, "(input):1:2: expected integer");
  EXPECT_THROW_MSG(NLReader().ReadString("g-1"),
      ampl::ParseError, "(input):1:2: expected integer");
  EXPECT_THROW_MSG(NLReader().ReadString("g10"),
      ampl::ParseError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(NLReader().ReadString(
      str(fmt::Format("g{}") << static_cast<unsigned>(INT_MAX) + 1)),
      ampl::ParseError, "(input):1:2: number is too big");
}

// The problem below represents the following AMPL problem in NL format
// without the first line (options):
//   var x >= 0;
//   minimize o: x;
const char *TEST_PROBLEM_NO_OPTIONS =
  " 1 0 1 0 0\n"
  " 0 0\n"
  " 0 0\n"
  " 0 0 0\n"
  " 0 0 0 1\n"
  " 0 0 0 0 0\n"
  " 0 1\n"
  " 0 0\n"
  " 0 0 0 0 0\n"
  "O0 0\n"
  "n0\n"
  "b\n"
  "2 0\n"
  "k0\n"
  "G0 1\n"
  "0 1\n";

struct TestNLHandler : ampl::NLHandler {
  ampl::NLHeader header;
  void HandleHeader(const ampl::NLHeader &h) { header = h; }
};

void CheckReadOptions(size_t num_options,
    size_t num_options_to_write, const int *options) {
  TestNLHandler handler;
  NLReader reader(&handler);
  fmt::Writer w;
  w << 'g' << num_options;
  for (size_t i = 0; i < num_options_to_write; ++i)
    w << ' ' << options[i];
  w << "\n" << TEST_PROBLEM_NO_OPTIONS;
  reader.ReadString(w.c_str());
  const ampl::NLHeader &header = handler.header;
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
  // TODO: test optional num_options
}

TEST(NLReaderTest, AMPLVBTol) {
  fmt::Writer w;
  w << "g2 0 3 4.2\n" << TEST_PROBLEM_NO_OPTIONS;
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(w.c_str());
  EXPECT_EQ(4.2, handler.header.ampl_vbtol);
}

TEST(NLReaderTest, IgnoreAMPLVBTol) {
  fmt::Writer w;
  w << "g2 0 0 4.2\n" << TEST_PROBLEM_NO_OPTIONS;
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(w.c_str());
  EXPECT_EQ(0, handler.header.ampl_vbtol);
}

TEST(NLReaderTest, OptionalAMPLVBTol) {
  fmt::Writer w;
  w << "g2 0 3\n" << TEST_PROBLEM_NO_OPTIONS;
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(w.c_str());
  EXPECT_EQ(0, handler.header.ampl_vbtol);
}

// TODO: more tests
}

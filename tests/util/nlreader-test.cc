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

void CheckReadOptions(size_t num_options, const int *options) {
  struct TestNLHandler : ampl::NLHandler {
    ampl::NLHeader header;
    void HandleHeader(const ampl::NLHeader &h) { header = h; }
  } handler;
  NLReader reader(&handler);
  // The problem below represents the following AMPL problem in NL format:
  //   var x >= 0;
  //   minimize o: x;
  fmt::Writer w;
  w << 'g' << num_options;
  for (size_t i = 0; i < num_options; ++i)
    w << ' ' << options[i];
  w << "\n"
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
  reader.ReadString(w.c_str());
  const ampl::NLHeader &header = handler.header;
  ASSERT_EQ(num_options, header.num_options);
  for (size_t i = 0; i < num_options; ++i)
    EXPECT_EQ(options[i], header.options[i]);
}

TEST(NLReaderTest, ReadOptions) {
  const int options[] = {3, 5, 7, 11, 13, 17, 19, 23, 29};
  size_t num_options = sizeof(options) / sizeof(*options);
  for (size_t i = 1; i < num_options; ++i)
    CheckReadOptions(i, options);
  // TODO: test the number of options different from actual
}

// TODO: more tests
}

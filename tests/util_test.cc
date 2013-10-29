/*
 Tests of test utilities

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

#include "tests/util.h"
#include "gtest/gtest.h"

namespace {

TEST(UtilTest, Split) {
  auto result = Split("abc", ' ');
  EXPECT_EQ(1, result.size());
  EXPECT_EQ("abc", result[0]);
  result = Split("a b c", ' ');
  EXPECT_EQ(3, result.size());
  EXPECT_EQ("a", result[0]);
  EXPECT_EQ("b", result[1]);
  EXPECT_EQ("c", result[2]);
  result = Split("abc ", ' ');
  EXPECT_EQ(2, result.size());
  EXPECT_EQ("abc", result[0]);
  EXPECT_EQ("", result[1]);
}
}

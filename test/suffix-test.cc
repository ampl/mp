/*
 Expression tests

 Copyright (C) 2015 AMPL Optimization Inc

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

#include <gtest/gtest.h>
#include "mp/suffix.h"

TEST(SuffixTest, Suffix) {
  mp::Suffix s;
  EXPECT_TRUE(s == 0);
}

TEST(SuffixTest, AddSuffix) {
  mp::SuffixSet suffixes;
  mp::Suffix s = suffixes.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
}

TEST(SuffixTest, ConversionToSuffix) {
  // Test that BasicSuffix is not convertible to Suffix&. If it was
  // convertible there would be an error because of an ambigous call.
  // The conversion is forbidden because it compromises type safety
  // as illustrated in the following example:
  //   auto i = suffixes.Add<int>("a", 0, 1);
  //   auto d = suffixes.Add<double>("b", 0, 1);
  //   mp::Suffix &s = i;
  //   s = d;
  struct Test {
    static void f(mp::Suffix) {}
    static void f(mp::Suffix &) {}
  };
  mp::SuffixSet suffixes;
  auto s = suffixes.Add<int>("a", 0, 1);
  Test::f(s);
}

// TODO: test VisitValues, BasicSuffix, ...

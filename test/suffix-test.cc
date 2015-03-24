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

#include "gmock/gmock.h"
#include "test-assert.h"
#include "mp/suffix.h"

class SuffixTest : public testing::Test {
 protected:
  mp::SuffixSet suffixes_;
};

TEST_F(SuffixTest, Suffix) {
  mp::Suffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
}

TEST_F(SuffixTest, IntSuffix) {
  mp::IntSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
}

TEST_F(SuffixTest, IntSuffixValue) {
  auto s = suffixes_.Add<int>("test", 11, 3);
  s.set_value(0, 42);
  EXPECT_EQ(42, s.value(0));
  s.set_value(2, 0);
  EXPECT_EQ(0, s.value(2));
  EXPECT_ASSERT(s.value(-1), "index out of bounds");
  EXPECT_ASSERT(s.value(3), "index out of bounds");
  EXPECT_ASSERT(s.set_value(-1, 0), "index out of bounds");
  EXPECT_ASSERT(s.set_value(3, 0), "index out of bounds");
}

TEST_F(SuffixTest, DoubleSuffix) {
  mp::DoubleSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<double>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
  s.set_value(0, 4.2);
  EXPECT_EQ(4.2, s.value(0));
  // TODO: test index checks in value() and set_value()
}

struct MockIntValueVisitor {
  MOCK_METHOD2(Visit, void (int index, int value));
};

TEST_F(SuffixTest, VisitIntSuffixValues) {
  mp::IntSuffix s = suffixes_.Add<int>("test", 0, 3);
  s.set_value(0, 42);
  s.set_value(1, 0);
  s.set_value(2, 11);
  MockIntValueVisitor v;
  EXPECT_CALL(v, Visit(0, 42));
  EXPECT_CALL(v, Visit(2, 11));
  s.VisitValues(v);
}

struct MockDoubleValueVisitor {
  MOCK_METHOD2(Visit, void (int index, double value));
};

TEST_F(SuffixTest, VisitDoubleSuffixValues) {
  mp::DoubleSuffix s = suffixes_.Add<double>("test", 0, 3);
  s.set_value(0, 4.2);
  s.set_value(1, 0);
  s.set_value(2, 1.1);
  MockDoubleValueVisitor v;
  EXPECT_CALL(v, Visit(0, 4.2));
  EXPECT_CALL(v, Visit(2, 1.1));
  s.VisitValues(v);
}

TEST_F(SuffixTest, ConversionToSuffix) {
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
  auto s = suffixes_.Add<int>("a", 0, 1);
  Test::f(s);
}

// TODO: test VisitValues, Cast, SuffixSet, SuffixManager

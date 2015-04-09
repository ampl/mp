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
#include "mock-allocator.h"
#include "test-assert.h"
#include "mp/error.h"
#include "mp/suffix.h"

using testing::_;
using testing::Matcher;
using testing::Return;

using mp::Suffix;
using mp::MutSuffix;

class SuffixTest : public testing::Test {
 protected:
  mp::SuffixSet suffixes_;
};

TEST_F(SuffixTest, Suffix) {
  Suffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
}

TEST_F(SuffixTest, MutSuffix) {
  MutSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
  Suffix csuf = s;
  csuf = s;
}

template <typename T>
struct IsInt { enum { VALUE = 0 }; };

template <>
struct IsInt<int> { enum { VALUE = 1 }; };

TEST_F(SuffixTest, IntSuffix) {
  mp::IntSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
  EXPECT_TRUE(IsInt<mp::IntSuffix::Type>::VALUE);
}

TEST_F(SuffixTest, MutIntSuffix) {
  mp::MutIntSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<int>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11, s.kind());
  EXPECT_EQ(222, s.num_values());
  EXPECT_TRUE(IsInt<mp::IntSuffix::Type>::VALUE);
  MutSuffix csuf = s;
  csuf = s;
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

template <typename T>
struct IsDouble { enum { VALUE = 0 }; };

template <>
struct IsDouble<double> { enum { VALUE = 1 }; };

TEST_F(SuffixTest, DoubleSuffix) {
  mp::DoubleSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<double>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11 | mp::suf::FLOAT, s.kind());
  EXPECT_EQ(222, s.num_values());
  EXPECT_TRUE(IsDouble<mp::DoubleSuffix::Type>::VALUE);
}

TEST_F(SuffixTest, MutDoubleSuffix) {
  mp::MutDoubleSuffix s;
  EXPECT_TRUE(s == 0);
  s = suffixes_.Add<double>("test", 11, 222);
  EXPECT_STREQ("test", s.name());
  EXPECT_EQ(11 | mp::suf::FLOAT, s.kind());
  EXPECT_EQ(222, s.num_values());
  EXPECT_TRUE(IsDouble<mp::DoubleSuffix::Type>::VALUE);
  MutSuffix csuf = s;
  csuf = s;
}

TEST_F(SuffixTest, DoubleSuffixValue) {
  auto s = suffixes_.Add<double>("test", 11, 3);
  s.set_value(0, 4.2);
  EXPECT_EQ(4.2, s.value(0));
  s.set_value(2, 0);
  EXPECT_EQ(0, s.value(2));
  EXPECT_ASSERT(s.value(-1), "index out of bounds");
  EXPECT_ASSERT(s.value(3), "index out of bounds");
  EXPECT_ASSERT(s.set_value(-1, 0), "index out of bounds");
  EXPECT_ASSERT(s.set_value(3, 0), "index out of bounds");
}

struct MockValueVisitor {
  MOCK_METHOD2(Visit, void (int index, int value));
  MOCK_METHOD2(Visit, void (int index, double value));
};

TEST_F(SuffixTest, VisitIntSuffixValues) {
  mp::MutIntSuffix s = suffixes_.Add<int>("test", 0, 3);
  s.set_value(0, 42);
  s.set_value(1, 0);
  s.set_value(2, 11);
  MockValueVisitor v;
  EXPECT_CALL(v, Visit(0, Matcher<int>(42)));
  EXPECT_CALL(v, Visit(2, Matcher<int>(11)));
  s.VisitValues(v);
}

TEST_F(SuffixTest, VisitDoubleSuffixValues) {
  mp::MutDoubleSuffix s = suffixes_.Add<double>("test", 0, 3);
  s.set_value(0, 4.2);
  s.set_value(1, 0);
  s.set_value(2, 1.1);
  MockValueVisitor v;
  EXPECT_CALL(v, Visit(0, Matcher<double>(4.2)));
  EXPECT_CALL(v, Visit(2, Matcher<double>(1.1)));
  s.VisitValues(v);
}

TEST_F(SuffixTest, VisitSuffixValues) {
  Suffix s;
  {
    auto is = suffixes_.Add<int>("is", 0, 3);
    is.set_value(0, 42);
    is.set_value(1, 0);
    is.set_value(2, 11);
    s = is;
  }
  {
    MockValueVisitor v;
    EXPECT_CALL(v, Visit(0, Matcher<int>(42)));
    EXPECT_CALL(v, Visit(2, Matcher<int>(11)));
    s.VisitValues(v);
  }
  {
    auto ds = suffixes_.Add<double>("ds", 0, 3);
    ds.set_value(0, 4.2);
    ds.set_value(1, 0);
    ds.set_value(2, 1.1);
    s = ds;
  }
  {
    MockValueVisitor v;
    EXPECT_CALL(v, Visit(0, Matcher<double>(4.2)));
    EXPECT_CALL(v, Visit(2, Matcher<double>(1.1)));
    s.VisitValues(v);
  }
}

TEST_F(SuffixTest, ConversionToSuffix) {
  // Test that BasicSuffix is not convertible to Suffix&. If it was
  // convertible there would be an error because of an ambiguous call.
  // The conversion is forbidden because it compromises type safety
  // as illustrated in the following example:
  //   auto i = suffixes.Add<int>("a", 0, 1);
  //   auto d = suffixes.Add<double>("b", 0, 1);
  //   Suffix &s = i;
  //   s = d;
  struct Test {
    static void f(Suffix) {}
    static void f(Suffix &) {}
  };
  auto s = suffixes_.Add<int>("a", 0, 1);
  Test::f(s);
  struct MutTest {
    static void f(mp::MutSuffix) {}
    static void f(mp::MutSuffix &) {}
  };
  MutTest::f(s);
}

TEST_F(SuffixTest, Is) {
  Suffix s = suffixes_.Add<int>("a", 0, 1);
  EXPECT_TRUE(mp::internal::Is<mp::IntSuffix>(s));
  EXPECT_FALSE(mp::internal::Is<mp::DoubleSuffix>(s));
}

TEST_F(SuffixTest, Cast) {
  mp::IntSuffix is = suffixes_.Add<int>("a", 0, 1);
  Suffix s = is;
  EXPECT_EQ(is, mp::Cast<mp::IntSuffix>(s));
  EXPECT_EQ(Suffix(), mp::Cast<mp::DoubleSuffix>(s));
}

TEST_F(SuffixTest, SuffixKindAgreesWithType) {
  EXPECT_EQ(0, suffixes_.Add<int>("a", 0, 1).kind());
  EXPECT_ASSERT(suffixes_.Add<int>("b", mp::suf::FLOAT, 1),
                "invalid suffix kind");
  EXPECT_EQ(mp::suf::FLOAT, suffixes_.Add<double>("c", 0, 1).kind());
  EXPECT_EQ(mp::suf::FLOAT,
            suffixes_.Add<double>("d", mp::suf::FLOAT, 1).kind());
}

TEST(SuffixSetTest, Empty) {
  mp::SuffixSet s;
  EXPECT_EQ(s.begin(), s.end());
  EXPECT_EQ(MutSuffix(), s.Find("test"));
}

TEST(SuffixSetTest, Add) {
  mp::SuffixSet s;
  MutSuffix suf = s.Add<int>("test", 0, 0);
  EXPECT_EQ(suf, s.Find("test"));
}

TEST(SuffixSetTest, NonNulTerminatedSuffixName) {
  mp::SuffixSet s;
  EXPECT_STREQ("foo", s.Add<int>(fmt::StringRef("foobar", 3), 0, 1).name());
}

TEST(SuffixSetTest, AddDuplicateSuffix) {
  mp::SuffixSet s;
  s.Add<int>("foo", 0, 1);
  EXPECT_THROW_MSG(s.Add<int>(fmt::StringRef("foobar", 3), 0, 1),
                   mp::Error, "duplicate suffix 'foo'");
}

TEST(SuffixSetTest, Find) {
  mp::SuffixSet s;
  auto suf = s.Add<int>("test", 0, 0);
  struct Test {
    static bool ismut(Suffix) { return false; }
    static bool ismut(MutSuffix) { return true; }
  };
  EXPECT_EQ(suf, s.Find("test"));
  EXPECT_TRUE(Test::ismut(s.Find("test")));
  const mp::SuffixSet &cs = s;
  EXPECT_EQ(suf, cs.Find("test"));
  EXPECT_FALSE(Test::ismut(cs.Find("test")));
}

TEST(SuffixSetTest, Iterator) {
  enum {NUM_SUFFIXES = 3};
  mp::SuffixSet s;
  Suffix suffixes[NUM_SUFFIXES] = {
    s.Add<int>("a", 0, 0), s.Add<int>("b", 0, 0), s.Add<int>("c", 0, 0)
  };
  mp::SuffixSet::iterator i = s.begin();
  EXPECT_EQ(suffixes[0], *i);
  EXPECT_STREQ("a", i->name());
  EXPECT_EQ(mp::SuffixSet::iterator(s.begin()), i);
  auto j = i;
  EXPECT_TRUE(i == j);
  j = i++;
  EXPECT_TRUE(i != j);
  EXPECT_EQ(suffixes[0], *j);
  EXPECT_EQ(suffixes[1], *i);
  j = ++i;
  EXPECT_EQ(j, i);
  EXPECT_EQ(suffixes[2], *i);
  ++i;
  EXPECT_EQ(s.end(), i);
}

TEST(ExprFactoryTest, ValueMemoryAllocation) {
  typedef testing::StrictMock<MockAllocator> Alloc;
  Alloc alloc;
  mp::BasicSuffixSet< AllocatorRef<Alloc> > s((AllocatorRef<Alloc>(&alloc)));
  char buffer[100];
  EXPECT_CALL(alloc, allocate(_)).WillOnce(Return(buffer));
  s.Add<int>("test", 0, 1);
  EXPECT_CALL(alloc, deallocate(buffer, _));
}

// TODO: test deallocation of suffix names
// TODO: test SuffixManager

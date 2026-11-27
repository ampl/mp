/*
 Tests of error classes and functions.

 Copyright (C) 2012 AMPL Optimization Inc

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
#include "mp/error.h"
#include "mp/nl-reader.h"

namespace {

TEST(ErrorTest, Error) {
  mp::Error e(fmt::CStringRef("test"));
  EXPECT_STREQ("test", e.what());
}

TEST(ErrorTest, ThrowError) {
  mp::Error error("");
  try {
    throw mp::Error(fmt::format("test {} message", "error"));
  } catch (const mp::Error &e) {
    error = e;
  }
  EXPECT_STREQ("test error message", error.what());
}

TEST(ErrorTest, UnwindsLocalObjects) {
  bool destroyed = false;
  struct Guard {
    bool& destroyed;
    ~Guard() { destroyed = true; }
  };
  try {
    Guard guard{destroyed};
    throw mp::Error("unwind");
  } catch (const mp::Error&) {
    EXPECT_TRUE(destroyed);
  }
}

TEST(ErrorTest, SolveCode) {
  mp::Error formatted(fmt::format("answer {}", 42));
  EXPECT_STREQ("answer 42", formatted.what());
  EXPECT_EQ(-1, formatted.exit_code());
  mp::Error coded(fmt::format("answer {}", 42), mp::SolveCode(200));
  EXPECT_STREQ("answer 42", coded.what());
  EXPECT_EQ(200, coded.exit_code());
  try {
    MP_INFEAS("test");
  } catch (const mp::Error& e) {
    EXPECT_EQ(200, e.exit_code());
  }
}

// An error without a solve result code must not get one:
// drivers write a code >= 0 to the .sol file as is
// (1 would read as "solved").
TEST(ErrorTest, NoCodeByDefault) {
  EXPECT_EQ(-1, mp::Error("plain").exit_code());
  EXPECT_EQ(-1, mp::UnsupportedError("feature {}", 1).exit_code());
  EXPECT_EQ(-1, mp::MakeUnsupportedError("feature {}", 1).exit_code());
  EXPECT_EQ(-1, mp::OptionError("option").exit_code());
  EXPECT_EQ(-1, mp::ReadError("f.nl", 1, 2, "bad {}", 3).exit_code());
  EXPECT_EQ(-1, mp::BinaryReadError("f.nl", 4, "bad").exit_code());
}

// The deprecated constructors keep their behavior.
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
# pragma warning(push)
# pragma warning(disable: 4996)
#endif
TEST(ErrorTest, DeprecatedConstructors) {
  mp::Error formatted("test {} message", "error");
  EXPECT_STREQ("test error message", formatted.what());
  EXPECT_EQ(-1, formatted.exit_code());
  mp::Error coded("literal {}", 200);    // the int is the code
  EXPECT_STREQ("literal {}", coded.what());
  EXPECT_EQ(200, coded.exit_code());
}
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic pop
#elif defined(_MSC_VER)
# pragma warning(pop)
#endif
}

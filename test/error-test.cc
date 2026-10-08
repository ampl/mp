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

namespace {

TEST(ErrorTest, Error) {
  mp::Error e(fmt::CStringRef("test"));
  EXPECT_STREQ("test", e.what());
}

TEST(ErrorTest, ThrowError) {
  mp::Error error("");
  try {
    throw mp::Error("test {} message", "error");
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

TEST(ErrorTest, IntegerMessageArgumentAndExplicitExitCode) {
  mp::Error formatted(fmt::format("answer {}", 42));
  EXPECT_STREQ("answer 42", formatted.what());
  mp::Error coded("literal {}", 200);
  EXPECT_STREQ("literal {}", coded.what());
  EXPECT_EQ(200, coded.exit_code());
  mp::Error both(fmt::format("answer {}", 42), 200);
  EXPECT_STREQ("answer 42", both.what());
  EXPECT_EQ(200, both.exit_code());
}
}

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

#include "gtest/gtest.h"
#include "solvers/util/error.h"

namespace {

TEST(ErrorTest, Error) {
  ampl::Error e(fmt::StringRef("test"));
  EXPECT_STREQ("test", e.what());
}

TEST(ErrorTest, ThrowError) {
  ampl::Error error("");
  try {
    ampl::ThrowError("test {} message") << "error";
  } catch (const ampl::Error &e) {
    error = e;
  }
  EXPECT_STREQ("test error message", error.what());
}
}

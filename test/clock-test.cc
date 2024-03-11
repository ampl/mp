/*
 Clock tests.

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
#include <chrono>
#include "mp/utils-clock.h"

#include <gtest/gtest.h>

using mp::GetTimeAndReset;
using std::chrono::duration;
using std::chrono::steady_clock;
using std::chrono::time_point;

namespace {

TEST(ClockTest, GetTimeAndReset) {
  steady_clock::time_point start;
  steady_clock::time_point t = start;
  double time = GetTimeAndReset(t);
  EXPECT_DOUBLE_EQ(
      std::chrono::duration_cast< duration<double> >(t - start).count(), time);
  EXPECT_GT(t, start);
}
}

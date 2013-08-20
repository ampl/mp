/*
 Clock tests.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "gtest/gtest.h"
#include "solvers/util/clock.h"

using ampl::ratio;
using ampl::duration;
using ampl::steady_clock;
using ampl::time_point;

TEST(ClockTest, Ratio) {
  typedef ratio<42, 77> Ratio;
  EXPECT_EQ(42, Ratio::num);
  EXPECT_EQ(77, Ratio::den);
}

TEST(ClockTest, Duration) {
  struct Rep {};
  struct Period {};
  typedef duration<Rep, Period> Duration;
  static_cast<Rep>(Duration::rep());
  static_cast<Period>(Duration::period());
  static_cast< ratio<1> >(duration<Rep>::period());
  EXPECT_EQ(0, duration<int>().count());
  EXPECT_EQ(42, duration<int>(42).count());
  EXPECT_EQ(10, (duration<int>(15) - duration<int>(5)).count());
}

TEST(ClockTest, DurationCast) {
  EXPECT_EQ(4200,
      ampl::duration_cast< duration<double> >(
          duration<int, ratio<100> >(42)).count());
  EXPECT_EQ(6,
      ampl::duration_cast< duration<double> >(
          duration<int, ratio<3, 5> >(10)).count());
}

TEST(ClockTest, TimePoint) {
  struct Duration {};
  static_cast<Duration>(time_point<void, Duration>::duration());
  struct Clock {
    typedef Duration duration;
  };
  static_cast<Duration>(time_point<Clock>::duration());
  typedef time_point<Clock, duration<int> > TimePoint;
  EXPECT_EQ(0, TimePoint().time_since_epoch().count());
  EXPECT_EQ(42, TimePoint(duration<int>(42)).time_since_epoch().count());
  EXPECT_EQ(12,
      (TimePoint(duration<int>(16)) - TimePoint(duration<int>(4))).count());
}

TEST(ClockTest, Now) {
  EXPECT_EQ(42, steady_clock::rep(42));
  static_cast< duration<steady_clock::rep, steady_clock::period> >(
      steady_clock::duration());
  static_cast< time_point<steady_clock> >(steady_clock::time_point());
  steady_clock::time_point start = steady_clock::now();
  steady_clock::time_point finish = steady_clock::now();
  EXPECT_GT(finish.time_since_epoch().count(),
      start.time_since_epoch().count());
}

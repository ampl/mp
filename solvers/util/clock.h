/*
 Clock classes similar to the ones from the C++11 time library (chrono),
 but with a reduced API. The goal is to eventually migrate to C++11.

 Copyright (C) 2013 AMPL Optimization LLC

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

#ifndef SOLVERS_UTIL_CLOCK_H_
#define SOLVERS_UTIL_CLOCK_H_

namespace ampl {

template <int N, int D = 1>
struct ratio {
  static const int num = N;
  static const int den = D;
};

// A time span.
template <typename Rep, typename Period = ratio<1> >
class duration {
 private:
  Rep rep_;

 public:
  typedef Rep rep;
  typedef Period period;

  duration() : rep_() {}
  explicit duration(const Rep &r) : rep_(r) {}

  rep count() const { return rep_; }

  friend duration operator-(const duration &lhs, const duration &rhs) {
    return duration(lhs.rep_ - rhs.rep_);
  }
};

template <typename DurationTo, typename Rep, typename Period>
DurationTo duration_cast(const duration<Rep, Period> &d) {
  // Only cast to duration<double> is supported so far.
  return duration<double>(
      static_cast<double>(d.count()) * Period::num / Period::den);
}

// A point in time relative to a clock's epoch.
template <typename Clock, typename Duration = typename Clock::duration>
class time_point {
 private:
  Duration d_;

 public:
  typedef Duration duration;

  time_point() {}
  explicit time_point(const duration &d) : d_(d) {}

  duration time_since_epoch() const { return d_; }

  friend duration operator-(const time_point &lhs, const time_point &rhs) {
    return duration(lhs.d_ - rhs.d_);
  }
};

// A clock for calculating time differences.
struct steady_clock {
  typedef long long rep;
  typedef ratio<1, 1000000000> period;
  typedef ampl::duration<rep, period> duration;
  typedef ampl::time_point<steady_clock> time_point;

  static time_point now() noexcept;
};
}

#endif  // SOLVERS_UTIL_CLOCK_H_

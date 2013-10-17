/*
 Clock classes similar to the ones from the C++11 time library (chrono),
 but with a reduced API. The goal is to eventually migrate to C++11.

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

#include "solvers/util/clock.h"

#if defined(__APPLE__)
# include <mach/mach_time.h>
#elif defined(_WIN32)
# include <windows.h>
#else
# include <time.h>
# include <sys/time.h>
#endif

#include <cassert>

namespace ampl {

#if defined(__APPLE__)

steady_clock::time_point steady_clock::now() {
  mach_timebase_info_data_t info;
  if (mach_timebase_info(&info)) {
    assert(0 && "mach_timebase_info failed");
    return time_point();
  }
  return time_point(duration(static_cast<rep>(
      static_cast<double>(mach_absolute_time()) * info.numer / info.denom)));
}

#elif defined(_WIN32)

double GetNanosecondsPerCount() {
  LARGE_INTEGER freq;
  typedef steady_clock::period period;
  if (!QueryPerformanceFrequency(&freq)) {
    assert(0 && "QueryPerformanceFrequency failed");
    return 0;
  }
  return static_cast<double>(period::den) / period::num / freq.QuadPart;
}

steady_clock::time_point steady_clock::now() {
  static const double NS_PER_COUNT = GetNanosecondsPerCount();
  if (!NS_PER_COUNT)
    return time_point();
  LARGE_INTEGER count;
  if (!QueryPerformanceCounter(&count)) {
    assert(0 && "QueryPerformanceCounter failed");
    return time_point();
  }
  return time_point(duration(static_cast<rep>(count.QuadPart * NS_PER_COUNT)));
}

#else

steady_clock::time_point steady_clock::now() {
  timespec ts;
#ifdef USE_CLOCK_GETTIME
  if (clock_gettime(CLOCK_MONOTONIC, &ts))
    assert(0 && "clock_gettime failed");
  return time_point(duration(
      static_cast<rep>(ts.tv_sec) * 1000000000 + ts.tv_nsec));
#else
  timeval tv;
  if (gettimeofday(&tv, 0))
    assert(0 && "gettimeofday failed");
  return time_point(duration(
      static_cast<rep>(ts.tv_sec) * 1000000000 + tv.tv_usec * 1000));
#endif
}

#endif

double GetTimeAndReset(steady_clock::time_point &t) {
  steady_clock::time_point now = steady_clock::now();
  double seconds = duration_cast< duration<double> >(now - t).count();
  t = now;
  return seconds;
}
}

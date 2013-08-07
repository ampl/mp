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

#include "solvers/util/clock.h"

#include <cassert>

#if defined(__APPLE__)
# include <mach/mach_time.h>
#elif defined(WIN32)
# include <windows.h>
#else
# include <time.h>
#endif

namespace ampl {

#if defined(__APPLE__)

steady_clock::time_point steady_clock::now() {
  mach_timebase_info_data_t info;
  if (mach_timebase_info(&info)) {
    assert(0 && "mach_timebase_info failed");
    return time_point();
  }
  return static_cast<steady_clock::rep>(
      static_cast<double>(mach_absolute_time()) * info.numer / info.denom);
}

#elif defined(WIN32)

double GetNanosecondsPerCount() {
  LARGE_INTEGER freq;
  typedef steady_clock::period period;
  return QueryPerformanceFrequency(&freq) ?
      static_cast<double>(period::den) / period::num / freq.QuadPart : 0;
}

steady_clock::time_point steady_clock::now() {
  static const double NS_PER_COUNT = GetNanosecondsPerCount();
  LARGE_INTEGER count;
  return time_point(NS_PER_COUNT && QueryPerformanceCounter(&count) ?
      static_cast<rep>(count.QuadPart * NS_PER_COUNT) : 0);
}

#else

steady_clock::time_point steady_clock::now() {
  timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts))
    assert(0 && "clock_gettime failed");
  return time_point(duration(
      static_cast<rep>(ts.tv_sec) * 1000000000 + ts.tv_nsec));
}

#endif
}

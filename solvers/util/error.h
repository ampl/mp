/*
 Exception classes.

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

#ifndef SOLVERS_UTIL_ERROR_H_
#define SOLVERS_UTIL_ERROR_H_

#include <stdexcept>
#include "solvers/util/format.h"

namespace ampl {

// A general error.
class Error : public std::runtime_error {
 public:
  explicit Error(fmt::StringRef message) : std::runtime_error(message) {}
  ~Error() throw() {}
};

struct Throw {
  void operator()(const fmt::Writer &w) const {
    throw Error(fmt::StringRef(w.c_str(), w.size()));
  }
};

// Throws Error with a formatted message.
inline fmt::Formatter<Throw> ThrowError(fmt::StringRef format) {
  fmt::Formatter<Throw> f(format);
  return f;
}

// Reports a system error without throwing an exception.
// Can be used to report errors from destructors.
void ReportSystemError(int error_code, const char *message) FMT_NOEXCEPT(true);

#ifdef _WIN32
// Reports a Windows error without throwing an exception.
// Can be used to report errors from destructors.
void ReportWinError(int error_code, const char *message) FMT_NOEXCEPT(true);
#endif
}

#endif  // SOLVERS_UTIL_ERROR_H_

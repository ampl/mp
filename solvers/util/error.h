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
};

struct Throw {
  void operator()(const fmt::Writer &w) const {
    throw Error(fmt::StringRef(w.c_str(), w.size()));
  }
};

// Throws Error with a formatted message.
inline fmt::Formatter<Throw> ThrowError(fmt::StringRef format) {
  return fmt::Formatter<Throw>(format);
}

// An error returned by the operating system or the language runtime,
// for example a file opening error.
class SystemError : public Error {
 private:
  int error_code_;

 public:
  SystemError(fmt::StringRef message, int error_code)
  : Error(message), error_code_(error_code) {}

  int error_code() const { return error_code_; }
};

class SystemThrow {
 private:
  int error_code_;

 public:
  explicit SystemThrow(int error_code) : error_code_(error_code) {}

  void operator()(const fmt::Writer &w) const;
};

// Throws SystemError with a code and a formatted message.
inline fmt::Formatter<SystemThrow> ThrowSystemError(
    int error_code, fmt::StringRef format) {
  return fmt::Formatter<SystemThrow>(format, SystemThrow(error_code));
}

// Reports a system error without throwing an exception.
// Can be used to report errors from destructors.
void ReportSystemError(int error_code, const char *message) FMT_NOEXCEPT(true);
}

#endif  // SOLVERS_UTIL_ERROR_H_

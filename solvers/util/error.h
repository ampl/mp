/*
 Exception classes.

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
  void operator()(const fmt::BasicFormatter<char> &f) const {
    throw Error(f.str());
  }
};

// Throws Error with formatted error message.
inline fmt::TempFormatter<Throw> ThrowError(fmt::StringRef format) {
  return fmt::TempFormatter<Throw>(format);
}
}

#endif  // SOLVERS_UTIL_ERROR_H_

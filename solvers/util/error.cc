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

#include "solvers/util/error.h"
#include "solvers/util/os.h"
#include <string.h>

void ampl::SystemThrow::operator()(const fmt::Writer &w) const {
#ifndef _WIN32
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  const char *message = 0;
  for (;;) {
    errno = 0;
#ifdef _GNU_SOURCE
    message = strerror_r(error_code_, &buffer[0], buffer.size());
#else
    strerror_r(error_code_, message = &buffer[0], buffer.size());
#endif
    if (errno == 0)
      break;
    if (errno != ERANGE) {
      // Can't get error message, print error code instead.
      throw SystemError(fmt::Format("{}: error code = {}")
          << w.c_str() << error_code_, error_code_);
    }
    buffer.resize(buffer.size() * 2);
  }
  throw SystemError(fmt::Format("{}: {}")
    << w.c_str() << message, error_code_);
#else
  // TODO: use FormatMessageW instead of strerror on Windows
  // and convert to UTF-8
  throw SystemError(fmt::Format("{}: error code = {}")
    << w.c_str() << error_code_, error_code_);
#endif
}

void ampl::LogSystemError(
    int error_code, const char *message) FMT_NOEXCEPT(true) {
  // TODO: log error without throwing an exception
}

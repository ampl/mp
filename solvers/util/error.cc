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

#include <errno.h>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32
# include <windows.h>
#endif

void ampl::SystemThrow::operator()(const fmt::Writer &w) const {
#ifndef _WIN32
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  const char *message = 0;
  for (;;) {
    errno = 0;
# ifdef _GNU_SOURCE
    message = strerror_r(error_code_, &buffer[0], buffer.size());
# else
    strerror_r(error_code_, message = &buffer[0], buffer.size());
# endif
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
  class String {
   private:
    LPWSTR str_;

   public:
    String() : str_() {}
    ~String() { LocalFree(str_); }
    LPWSTR *ptr() const { return &str_; }
    LPWSTR c_str() const { return str_; }
  };
  String message;
  if (!FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER |
      FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, 0,
      error_code_, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
      reinterpret_cast<LPWSTR>(message.ptr()), 0, 0)) {
    // Can't get error message, print error code instead.
    throw SystemError(fmt::Format("{}: error code = {}")
        << w.c_str() << error_code_, error_code_);
  }
  throw SystemError(fmt::Format("{}: {}")
    << w.c_str() << UTF16ToUTF8(message.c_str()), error_code_);
#endif
}

void ampl::LogSystemError(
    int error_code, const char *message) FMT_NOEXCEPT(true) {
  // TODO: log error without throwing an exception
}

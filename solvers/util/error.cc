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
#include <string.h>

#ifdef _WIN32
# include <windows.h>
#endif

namespace {

void FormatSystemErrorMessage(
    fmt::Writer &w, int error_code, fmt::StringRef message) {
#ifndef _WIN32
  fmt::internal::Array<char, ampl::BUFFER_SIZE> buffer;
  buffer.resize(ampl::BUFFER_SIZE);
  char *system_message = 0;
  for (;;) {
    errno = 0;
# ifdef _GNU_SOURCE
    system_message = strerror_r(error_code, &buffer[0], buffer.size());
# else
    strerror_r(error_code, system_message = &buffer[0], buffer.size());
# endif
    if (errno == 0)
      break;
    if (errno != ERANGE) {
      // Can't get error message, report error code instead.
      w.Format("{}: error code = {}") << message << error_code;
      return;
    }
    buffer.resize(buffer.size() * 2);
  }
  w.Format("{}: {}") << message << system_message;
#else
  class String {
   private:
    LPWSTR str_;

   public:
    String() : str_() {}
    ~String() { LocalFree(str_); }
    LPWSTR *ptr() { return &str_; }
    LPCWSTR c_str() const { return str_; }
  };
  String system_message;
  if (FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER |
      FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, 0,
      error_code, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
      reinterpret_cast<LPWSTR>(system_message.ptr()), 0, 0)) {
    ampl::UTF16ToUTF8 utf8_message;
    if (!utf8_message.Convert(system_message.c_str())) {
      w.Format("{}: {}") << message << utf8_message;
      return;
    }
  }
  // Can't get error message, report error code instead.
  w.Format("{}: error code = {}") << message << error_code;
#endif
}
}

void ampl::SystemThrow::operator()(const fmt::Writer &w) const {
  fmt::Writer message;
  FormatSystemErrorMessage(message, error_code_, w.c_str());
  throw SystemError(message.c_str(), error_code_);
}

void ampl::ReportSystemError(
    int error_code, const char *message) FMT_NOEXCEPT(true) {
  try {
    fmt::Writer full_message;
    FormatSystemErrorMessage(full_message, error_code, message);
    std::fwrite(full_message.c_str(), full_message.size(), 1, stderr);
    std::fputc('\n', stderr);
  } catch (...) {}
}

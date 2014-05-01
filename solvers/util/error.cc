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

namespace {
typedef void (*FormatFunc)(fmt::Writer &, int , fmt::StringRef);

void ReportError(
    FormatFunc func, int error_code, const char *message) FMT_NOEXCEPT(true) {
  try {
    fmt::Writer full_message;
    func(full_message, error_code, message);
    std::fwrite(full_message.c_str(), full_message.size(), 1, stderr);
    std::fputc('\n', stderr);
  } catch (...) {}
}
}

void ampl::ReportSystemError(
    int error_code, const char *message) FMT_NOEXCEPT(true) {
  ReportError(fmt::internal::FormatSystemErrorMessage, error_code, message);
}

#ifdef _WIN32
void ampl::ReportWinError(
    int error_code, const char *message) FMT_NOEXCEPT(true) {
  ReportError(fmt::internal::FormatWinErrorMessage, error_code, message);
}
#endif

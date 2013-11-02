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
#include <cstring>

void ampl::SystemThrow::operator()(const fmt::Writer &w) const {
#ifndef _WIN32
  // TODO: use strerror_r
  throw SystemError(fmt::Format("{}: {}")
    << w.c_str() << std::strerror(error_code_), error_code_);
#else
  // TODO: use FormatMessageW instead of strerror on Windows
  // and convert to UTF-8
  throw SystemError(fmt::Format("{}: error code = {}")
    << w.c_str() << error_code_, error_code_);
#endif
}

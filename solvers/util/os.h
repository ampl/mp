/*
 Operating system dependent functionality.

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

#ifndef SOLVERS_UTIL_OS_H_
#define SOLVERS_UTIL_OS_H_

#include <string>

#include "solvers/util/format.h"
#include "solvers/util/noncopyable.h"

namespace ampl {

// This class provides a subset of boost filesystem's path API.
class path {
 private:
  std::string str_;
  
  inline std::size_t FindLastSep() const {
#ifdef _WIN32
    const char *sep = "/\\";
#else
    const char sep = '/';
#endif
    return str_.find_last_of(sep);
  }

 public:
#ifdef _WIN32
  static const char preferred_separator = '\\';
#else
  static const char preferred_separator = '/';
#endif

  path() {}
  explicit path(const std::string &s): str_(s) {}
  path(const char *begin, const char *end) : str_(begin, end) {}

  const std::string &string() const { return str_; }
  
  path filename() const {
    size_t last_sep = FindLastSep();
    return last_sep == std::string::npos ?
        *this : path(str_.substr(last_sep + 1));
  }

  path &remove_filename() {
    size_t last_sep = FindLastSep();
    if (last_sep == 0)
      ++last_sep;
    str_.resize(last_sep != std::string::npos ? last_sep : 0);
    return *this;
  }

  // Returns a path to the system-specific temporary directory.
  static path temp_directory_path();
};

// Returns the path to the currently running executable file.
// Throws Error on error.
path GetExecutablePath();

class MemoryMappedFile : Noncopyable {
 private:
  char *start_;
  unsigned long long size_;

 public:
  explicit MemoryMappedFile(fmt::StringRef filename);
  ~MemoryMappedFile();

  const char *start() const { return start_; }
  unsigned long long size() const { return size_; }
};

// The default buffer size.
enum { BUFFER_SIZE = 500 };
}

#endif  // SOLVERS_UTIL_OS_H_

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

#include "solvers/util/noncopyable.h"

namespace ampl {

class path {
 private:
  std::string str_;

 public:
  path() {}
  path(const char *begin, const char *end) : str_(begin, end) {}

  const std::string &string() const { return str_; }

  path &remove_filename() {
    size_t last_sep = str_.find_last_of('/');
    str_.resize(last_sep != std::string::npos ? last_sep : 0);
    return *this;
  }
};

// Returns the path to the currently running executable file.
// Throws Error on error.
path GetExecutablePath();

class MemoryMappedFile : Noncopyable {
 private:
  char *start_;
  std::size_t length_;

 public:
  explicit MemoryMappedFile(const char *filename);
  ~MemoryMappedFile();

  const char *start() const { return start_; }
};

#ifdef _WIN32
// A converter from UTF-8 to UTF-16.
class UTF8ToUTF16 {
 private:
  fmt::internal::Array<wchar_t, BUFFER_SIZE> buffer_;

 public:
  explicit UTF8ToUTF16(const char *s);
  operator const wchar_t*() const { return &buffer_[0]; }
};

// A converter from UTF-16 to UTF-8.
class UTF16ToUTF8 {
 private:
  fmt::internal::Array<char, BUFFER_SIZE> buffer_;

 public:
  explicit UTF16ToUTF8(const wchar_t *s);
  operator const char*() const { return &buffer_[0]; }
  size_t size() const { return buffer_.size(); }
};
#endif
}

#endif  // SOLVERS_UTIL_OS_H_

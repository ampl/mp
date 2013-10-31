/*
 Test utilities

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

#ifndef TESTS_UTIL_H_
#define TESTS_UTIL_H_

#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>

#include "solvers/util/format.h"

std::string ReadFile(fmt::StringRef name);
void WriteFile(fmt::StringRef name, fmt::StringRef data);

inline std::string FixPath(fmt::StringRef s) {
#ifdef _WIN32
  std::string fixed = s;
  std::replace(fixed.begin(), fixed.end(), '/', '\\');
  return fixed;
#else
  return s;
#endif
}

// Redirects Stderr to a file.
class StderrRedirect {
 private:
  std::FILE *saved_stderr;

 public:
  explicit StderrRedirect(const char *filename);
  ~StderrRedirect();
};

// Changes the current working directory. Throws Error on error.
void ChangeDirectory(fmt::StringRef path);

// Executes a shell command. Throws Error on error.
void ExecuteShellCommand(fmt::StringRef command);

// Splits the string into an array of substrings.
std::vector<std::string> Split(const std::string &s, char sep);

#endif  // TESTS_UTIL_H_

#ifndef TESTS_UTIL_H_
#define TESTS_UTIL_H_

#include <cstdio>
#include <algorithm>
#include <string>

#include "solvers/util/format.h"

std::string ReadFile(const char *name);
void WriteFile(const char *name, const char *data);

inline std::string FixPath(fmt::StringRef s) {
#ifdef WIN32
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

#endif  // TESTS_UTIL_H_

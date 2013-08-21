/*
 Functions and classes for working with files and directories.

 Copyright (C) 2013 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "solvers/util/filesystem.h"

#include "solvers/util/error.h"

namespace {
enum { BUFFER_SIZE = 500 };
}

#if defined(__APPLE__)

#include <mach-o/dyld.h>

ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  uint32_t size = BUFFER_SIZE;
  buffer.resize(size);
  if (_NSGetExecutablePath(&buffer[0], &size) != 0) {
    buffer.resize(size);
    if (_NSGetExecutablePath(&buffer[0], &size) != 0)
      assert(0 && "_NSGetExecutablePath failed");
  }
  if (size == BUFFER_SIZE)
    size = std::strlen(&buffer[0]);
  const char *s = &buffer[0];
  return path(s, s + size);
}

#elif defined(WIN32)

#include <windows.h>

ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  DWORD size = 0;
  for (;;) {
    size = GetModuleFileNameA(0, &buffer[0], buffer.size());
    if (size == 0)
      ThrowError("GetModuleFileName failed, error code = {}") << GetLastError();
    if (size != buffer.size()) break;
    buffer.resize(2 * buffer.size());
  }
  std::replace(&buffer[0], &buffer[0] + size, '\\', '/');
  const char *s = &buffer[0];
  return path(s, s + size);
}

#else

#include <errno.h>
#include <unistd.h>

ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  ssize_t size = 0;
  for (;;) {
    size = readlink("/proc/self/exe", &buffer[0], buffer.size());
    if (size < 0) {
      ThrowError("readlink failed for /proc/self/exe, error code = {}")
          << errno;
    }
    if (static_cast<std::size_t>(size) != buffer.size()) break;
    buffer.resize(2 * buffer.size());
  }
  const char *s = &buffer[0];
  return path(s, s + size);
}

#endif

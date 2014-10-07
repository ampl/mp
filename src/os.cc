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

#include "mp/os.h"
#include "mp/error.h"
#include "mp/posix.h"

#include <cerrno>
#include <cstdlib>
#include <algorithm>

#ifndef _WIN32
# include <sys/mman.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <fcntl.h>
#endif

#if defined(__APPLE__)
# include <mach-o/dyld.h>
#elif defined(_WIN32)
# include <windows.h>
# include <io.h>
# undef min
#else
# include <unistd.h>
#endif

#undef getenv

using std::size_t;
using fmt::SystemError;
using mp::path;

// Workaround for a bug in MSVC.
// http://connect.microsoft.com/VisualStudio/feedback/details/
// 786583/in-class-static-const-member-initialization-and-lnk2005
#ifndef _MSC_EXTENSIONS
const char path::preferred_separator;
#endif

#ifndef _WIN32

#ifdef __APPLE__

// Mac OS X implementation.
path mp::GetExecutablePath() {
  fmt::internal::MemoryBuffer<char, BUFFER_SIZE> buffer;
  uint32_t size = BUFFER_SIZE;
  buffer.resize(size);
  if (_NSGetExecutablePath(&buffer[0], &size) != 0) {
    buffer.resize(size);
    if (_NSGetExecutablePath(&buffer[0], &size) != 0)
      throw SystemError(errno, "cannot get executable path");
  }
  if (size == BUFFER_SIZE)
    size = std::strlen(&buffer[0]);
  const char *s = &buffer[0];
  return path(s, s + size);
}

#else

// Linux implementation.
path mp::GetExecutablePath() {
  fmt::internal::MemoryBuffer<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  ssize_t size = 0;
  for (;;) {
    size = readlink("/proc/self/exe", &buffer[0], buffer.size());
    if (size < 0)
      throw SystemError(errno, "cannot get executable path");
    if (static_cast<size_t>(size) != buffer.size()) break;
    buffer.resize(2 * buffer.size());
  }
  const char *s = &buffer[0];
  return path(s, s + size);
}

#endif

// POSIX implementation.

path path::temp_directory_path() {
  const char *dir = std::getenv("TMPDIR");
  if (!dir) {
#ifdef P_tmpdir
    dir = P_tmpdir;
#else
    dir = "/tmp";
#endif
  }
  return path(dir);
}

mp::internal::MemoryMappedFileBase::MemoryMappedFileBase(
    int fd, std::size_t size)
: start_(), size_(size) {
  start_ = reinterpret_cast<char*>(
      mmap(0, size_, PROT_READ, MAP_FILE | MAP_PRIVATE, fd, 0));
  if (start_ == MAP_FAILED)
    throw SystemError(errno, "cannot map file");
}

mp::internal::MemoryMappedFileBase::~MemoryMappedFileBase() {
  if (munmap(start_, size_) == -1)
    fmt::report_system_error(errno, "cannot unmap file");
}

#else

// Windows implementation.

using fmt::WindowsError;

path path::temp_directory_path() {
  enum { BUFFER_SIZE = MAX_PATH + 1 };
  wchar_t buffer[BUFFER_SIZE];
  DWORD result = GetTempPathW(BUFFER_SIZE, &buffer[0]);
  if (result == 0) {
    throw WindowsError(
      GetLastError(), "cannot get path to the temporary directory");
  }
  assert(result <= BUFFER_SIZE);
  buffer[BUFFER_SIZE - 1] = L'\0';
  fmt::internal::UTF16ToUTF8 utf8_str(buffer);
  const char *s = utf8_str.c_str();
  return path(s, s + utf8_str.size());
}

path mp::GetExecutablePath() {
  fmt::internal::MemoryBuffer<wchar_t, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  DWORD size = 0;
  for (;;) {
    size = GetModuleFileNameW(0, &buffer[0], static_cast<DWORD>(buffer.size()));
    if (size == 0)
      throw WindowsError(GetLastError(), "cannot get executable path");
    if (size < buffer.size()) break;
    buffer.resize(2 * buffer.size());
  }
  fmt::internal::UTF16ToUTF8 utf8_str(&buffer[0]);
  const char *s = utf8_str.c_str();
  return path(s, s + utf8_str.size());
}

mp::internal::MemoryMappedFileBase::MemoryMappedFileBase(
    int fd, std::size_t size)
: start_(), size_(size) {
  class Handle {
    HANDLE handle_;
    Handle(const Handle &);
    void operator=(const Handle &);
   public:
    explicit Handle(HANDLE h) : handle_(h) {}
    ~Handle() { CloseHandle(handle_); }
    operator HANDLE() const { return handle_; }
  };
  HANDLE handle = reinterpret_cast<HANDLE>(_get_osfhandle(fd));
  if (handle == INVALID_HANDLE_VALUE)
    throw SystemError(errno, "cannot get file handle");
  Handle mapping(CreateFileMappingW(handle, 0, PAGE_READONLY, 0, 0, 0));
  if (!mapping)
    throw WindowsError(GetLastError(), "cannot create file mapping");
  start_ = reinterpret_cast<char*>(
      MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0));
  if (!start_)
    throw WindowsError(GetLastError(), "cannot map file");
}

mp::internal::MemoryMappedFileBase::~MemoryMappedFileBase() {
  if (!UnmapViewOfFile(start_))
    throw WindowsError(GetLastError(), "cannot unmap file");
}

#endif

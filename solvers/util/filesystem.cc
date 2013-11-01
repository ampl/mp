/*
 Functions and classes for working with files and directories.

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

#include "solvers/util/filesystem.h"

#include <cerrno>

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
#else
# include <unistd.h>
#endif

#include "solvers/util/error.h"

namespace {
enum { BUFFER_SIZE = 500 };
}

#if defined(__APPLE__)

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

#elif defined(_WIN32)

ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  DWORD size = 0;
  for (;;) {
    size = GetModuleFileNameA(0, &buffer[0], static_cast<DWORD>(buffer.size()));
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

#ifndef _WIN32

ampl::MemoryMappedFile::MemoryMappedFile(const char *filename)
: start_(), length_() {
  class File : Noncopyable {
    int fd_;
   public:
    explicit File(const char *filename) : fd_(open(filename, O_RDONLY)) {
      if (fd_ == -1)
        ThrowSystemError(errno, "cannot open file {}") << filename;
    }
    ~File() { close(fd_); }
    operator int() const { return fd_; }
  };
  File file(filename);
  struct stat file_stat = {};
  if (fstat(file, &file_stat) == -1)
    ThrowSystemError(errno, "cannot get attributes of file {}") << filename;
  length_ = file_stat.st_size;
  long pagesize = sysconf(_SC_PAGESIZE);
  std::size_t extra_bytes = length_ % pagesize;
  if (extra_bytes == 0) {
    // TODO: don' use mmap
  }
  // Round length up to a multiple of the memory page size.
  length_ += pagesize - extra_bytes;
  start_ = reinterpret_cast<char*>(
      mmap(0, length_, PROT_READ, MAP_FILE | MAP_PRIVATE, file, 0));
  if (start_ == MAP_FAILED)
    ThrowSystemError(errno, "cannot map file {}") << filename;
}

ampl::MemoryMappedFile::~MemoryMappedFile() {
  munmap(start_, length_);
}

#else

namespace {

// A converter from UTF-8 to UTF-16.
class UTF8ToUTF16 {
 private:
  fmt::internal::Array<WCHAR, 500> buffer_;

 public:
  explicit UTF8ToUTF16(const char *s);
  operator const WCHAR*() const { return &buffer_[0]; }
};

UTF8ToUTF16::UTF8ToUTF16(const char *s) {
  int length = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s, -1, 0, 0);
  if (length == 0) {
    ThrowSystemError(GetLastError(),
        "Cannot convert string from UTF-8 to UTF-16");
  }
  buffer_.resize(length);
  length = MultiByteToWideChar(
    CP_UTF8, MB_ERR_INVALID_CHARS, s, -1, &buffer_[0], length);
  if (length == 0) {
    ThrowSystemError(GetLastError(),
        "Cannot convert string from UTF-8 to UTF-16");
  }
}
}

ampl::MemoryMappedFile::MemoryMappedFile(const char *filename)
: start_(), length_() {
  class Handle : Noncopyable {
    HANDLE handle_;
   public:
    explicit Handle(HANDLE h) : handle_(h) {}
    ~Handle() { CloseHandle(handle_); }
    operator HANDLE() const { return handle_; }
  };
  Handle file(CreateFileW(UTF8ToUTF16(filename), GENERIC_READ,
      FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0));
  // TODO: get file size and check if it is not a multiple of the page size
  Handle mapping(CreateFileMappingW(file, 0, PAGE_READONLY, 0, 0, 0));
  if (!mapping)
    ThrowSystemError(GetLastError(), "cannot map file {}") << filename;
}

ampl::MemoryMappedFile::~MemoryMappedFile() {
  // TODO: unmap
}

#endif

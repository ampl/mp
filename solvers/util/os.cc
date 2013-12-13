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

#include "solvers/util/os.h"

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
# undef ERROR
#else
# include <unistd.h>
#endif

#include "solvers/util/error.h"

using std::size_t;

namespace {

// Round n up to a multiple of page_size.
size_t RoundUpToMultipleOf(size_t n, size_t page_size) {
  size_t extra_bytes = n % page_size;
  if (extra_bytes == 0) {
    // TODO: don' use mmap
  }
  return n += page_size - extra_bytes;
}
}

#ifndef _WIN32

#ifdef __APPLE__

// Mac OS X implementation.
ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  uint32_t size = BUFFER_SIZE;
  buffer.resize(size);
  if (_NSGetExecutablePath(&buffer[0], &size) != 0) {
    buffer.resize(size);
    if (_NSGetExecutablePath(&buffer[0], &size) != 0)
      ThrowSystemError(errno, "cannot get executable path");
  }
  if (size == BUFFER_SIZE)
    size = std::strlen(&buffer[0]);
  const char *s = &buffer[0];
  return path(s, s + size);
}

#else

// Linux implementation.
ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<char, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  ssize_t size = 0;
  for (;;) {
    size = readlink("/proc/self/exe", &buffer[0], buffer.size());
    if (size < 0)
      ThrowSystemError(errno, "cannot get executable path");
    if (static_cast<size_t>(size) != buffer.size()) break;
    buffer.resize(2 * buffer.size());
  }
  const char *s = &buffer[0];
  return path(s, s + size);
}

#endif

// POSIX implementation of MemoryMappedFile.

ampl::MemoryMappedFile::MemoryMappedFile(fmt::StringRef filename)
: start_(), size_() {
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

  // Open file and check that its size is not a multiple of memory page size.
  File file(filename.c_str());
  struct stat file_stat = {};
  if (fstat(file, &file_stat) == -1)
    ThrowSystemError(errno, "cannot get attributes of file {}") << filename;
  size_ = file_stat.st_size;
  // TODO: don't use mmap if file size is a multiple of page size
  //size_t full_size = RoundUpToMultipleOf(size_, sysconf(_SC_PAGESIZE));

  // Map file to memory.
  start_ = reinterpret_cast<char*>(
      mmap(0, size_, PROT_READ, MAP_FILE | MAP_PRIVATE, file, 0));
  if (start_ == MAP_FAILED)
    ThrowSystemError(errno, "cannot map file {}") << filename;
}

ampl::MemoryMappedFile::~MemoryMappedFile() {
  if (munmap(start_, size_) == -1)
    ReportSystemError(errno, "cannot unmap file");
}

#else

// Windows implementation:

ampl::UTF8ToUTF16::UTF8ToUTF16(fmt::StringRef s) {
  int length = MultiByteToWideChar(
      CP_UTF8, MB_ERR_INVALID_CHARS, s.c_str(), -1, 0, 0);
  static const char ERROR[] = "cannot convert string from UTF-8 to UTF-16";
  if (length == 0)
    ThrowSystemError(GetLastError(), ERROR);
  buffer_.resize(length);
  length = MultiByteToWideChar(
    CP_UTF8, MB_ERR_INVALID_CHARS, s.c_str(), -1, &buffer_[0], length);
  if (length == 0)
    ThrowSystemError(GetLastError(), ERROR);
}

ampl::UTF16ToUTF8::UTF16ToUTF8(fmt::WStringRef s) {
  if (int error_code = Convert(s)) {
    ThrowSystemError(GetLastError(),
        "cannot convert string from UTF-16 to UTF-8");
  }
}

int ampl::UTF16ToUTF8::Convert(fmt::WStringRef s) {
  int length = WideCharToMultiByte(CP_UTF8, 0, s.c_str(), -1, 0, 0, 0, 0);
  if (length == 0)
    return GetLastError();
  buffer_.resize(length);
  length = WideCharToMultiByte(
    CP_UTF8, 0, s.c_str(), -1, &buffer_[0], length, 0, 0);
  if (length == 0)
    return GetLastError();
  return 0;
}

ampl::path ampl::GetExecutablePath() {
  fmt::internal::Array<wchar_t, BUFFER_SIZE> buffer;
  buffer.resize(BUFFER_SIZE);
  DWORD size = 0;
  for (;;) {
    size = GetModuleFileNameW(0, &buffer[0], static_cast<DWORD>(buffer.size()));
    if (size == 0)
      ThrowSystemError(GetLastError(), "cannot get executable path");
    if (size < buffer.size()) break;
    buffer.resize(2 * buffer.size());
  }
  std::replace(&buffer[0], &buffer[0] + size, L'\\', L'/');
  UTF16ToUTF8 utf16_str(&buffer[0]);
  const char *s = utf16_str;
  return path(s, s + utf16_str.size());
}

ampl::MemoryMappedFile::MemoryMappedFile(fmt::StringRef filename)
: start_(), size_() {
  class Handle : Noncopyable {
    HANDLE handle_;
   public:
    explicit Handle(HANDLE h) : handle_(h) {}
    ~Handle() { CloseHandle(handle_); }
    operator HANDLE() const { return handle_; }
  };

  // Open file.
  Handle file(CreateFileW(UTF8ToUTF16(filename.c_str()), GENERIC_READ,
      FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0));
  if (file == INVALID_HANDLE_VALUE)
    ThrowSystemError(GetLastError(), "cannot open file {}") << filename;

  // Get file size and check if it is not a multiple of a memory page size.
  LARGE_INTEGER size = {};
  if (!GetFileSizeEx(file, &size))
    ThrowSystemError(GetLastError(), "cannot get size of file {}") << filename;
  SYSTEM_INFO si = {};
  GetSystemInfo(&si);
  size_ = size.QuadPart;
  // TODO: don't use mmap if file size is a multiple of page size
  //size_ = RoundUpToMultipleOf(size.QuadPart, si.dwPageSize);

  // Map file to memory.
  Handle mapping(CreateFileMappingW(file, 0, PAGE_READONLY, 0, 0, 0));
  if (!mapping) {
    ThrowSystemError(GetLastError(),
        "cannot create file mapping for {}") << filename;
  }
  start_ = reinterpret_cast<char*>(
      MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0));
  if (!start_)
    ThrowSystemError(GetLastError(), "cannot map file {}");
}

ampl::MemoryMappedFile::~MemoryMappedFile() {
  if (!UnmapViewOfFile(start_))
    ReportSystemError(GetLastError(), "cannot unmap file");
}

#endif

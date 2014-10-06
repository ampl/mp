/*
 Tests of operating system dependent functionality.

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

#include <gtest/gtest.h>
#include "mp/error.h"
#include "mp/os.h"
#include "util.h"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifndef _WIN32
# include <unistd.h>
#else
# include <windows.h>
#endif

using fmt::File;
using mp::MemoryMappedFile;
using mp::path;

using std::string;

namespace {

class PathTest : public ::testing::TestWithParam<char> {};

TEST(PathTest, PreferredSeparator) {
#ifdef _WIN32
  EXPECT_EQ('\\', path::preferred_separator);
#else
  EXPECT_EQ('/', path::preferred_separator);
#endif
}

TEST_P(PathTest, PathCtor) {
  EXPECT_EQ("", path().string());
  std::string s = FixPath("/some/path", GetParam());
  path p(&s[0], &s[0] + s.size());
  EXPECT_EQ(s, p.string());
  EXPECT_EQ(s, path(s).string());
}

TEST_P(PathTest, RemoveFilename) {
  path p;
  EXPECT_EQ("", p.remove_filename().string());
  EXPECT_EQ("", p.string());
  std::string s = FixPath("/somewhere/out/in/space", GetParam());
  p = path(s);
  s = FixPath("/somewhere/out/in", GetParam());
  EXPECT_EQ(s, p.remove_filename().string());
  EXPECT_EQ(s, p.string());
  s = FixPath("/somewhere/out", GetParam());
  EXPECT_EQ(s, p.remove_filename().string());
  EXPECT_EQ(s, p.string());
  s = FixPath("/", GetParam());
  EXPECT_EQ(s, path(s).remove_filename().string());
  EXPECT_EQ("", path("test").remove_filename().string());
}

TEST_P(PathTest, Filename) {
  path p;
  EXPECT_EQ("", p.filename().string());
  EXPECT_EQ("", p.string());
  std::string s = FixPath("/somewhere/out/in/space", GetParam());
  p = path(s);
  EXPECT_EQ("space", p.filename().string());
  EXPECT_EQ(s, p.string());
  EXPECT_EQ("", path(FixPath("/", GetParam())).filename().string());
  EXPECT_EQ("test", path("test").filename().string());
}

INSTANTIATE_TEST_CASE_P(POSIX, PathTest, ::testing::Values('/'));

#ifdef _WIN32
INSTANTIATE_TEST_CASE_P(Win32, PathTest, ::testing::Values('\\'));
#endif

TEST(PathTest, TempDirectoryPath) {
#ifndef _WIN32
  const char *dir = std::getenv("TMPDIR");
#ifndef P_tmpdir
# define P_tmpdir "/tmp"
#endif
  EXPECT_EQ(dir ? dir : P_tmpdir, path::temp_directory_path().string());
#else
  wchar_t buffer[MAX_PATH + 1];
  DWORD result = GetTempPathW(MAX_PATH + 1, buffer);
  EXPECT_GT(result, 0u);
  EXPECT_LE(result, static_cast<DWORD>(MAX_PATH));
  EXPECT_EQ(fmt::internal::UTF16ToUTF8(buffer).str(),
      path::temp_directory_path().string());
#endif
}

TEST(OSTest, GetExecutablePath) {
  string filename = "os-test";
#ifdef _WIN32
  filename += ".exe";
#endif
  path p = mp::GetExecutablePath();
  EXPECT_EQ(filename, p.filename().string());
  p.remove_filename();
#if defined(_WIN32) && !defined(__MINGW32__)
  p.remove_filename();
#endif
  EXPECT_EQ("bin", p.filename().string());
}

// Creates a new link for a file or copies the file if the system doesn't
// support symlinks. Both filenames are UTF-8 encoded.
void LinkFile(fmt::StringRef filename, fmt::StringRef linkname) {
#ifndef _WIN32
  int result = link(filename.c_str(), linkname.c_str());
  if (result && errno != EEXIST) {
    throw fmt::SystemError(errno,
        "cannot create a symlink from {} to {}", filename, linkname);
  }
#else
  using fmt::internal::UTF8ToUTF16;
  if (!CopyFileW(UTF8ToUTF16(filename).c_str(),
      UTF8ToUTF16(linkname).c_str(), FALSE)) {
    throw fmt::WindowsError(GetLastError(), "cannot copy file {} to {}",
      filename.c_str(), linkname.c_str());
  }
#endif
}

TEST(OSTest, LinkFile) {
  WriteFile("out", "some content");
  std::remove("out-link");
  LinkFile("out", "out-link");
  EXPECT_EQ("some content", ReadFile("out-link"));
  // Test creating link if file exists.
  LinkFile("out", "out-link");
  EXPECT_EQ("some content", ReadFile("out-link"));
}

std::string GetExecutableDir() {
  return mp::GetExecutablePath().remove_filename().string();
}

TEST(OSTest, GetExecutablePathUnicode) {
  // Neither CMake nor NMake handle Unicode paths properly on Windows,
  // so copy test executable ourselves.
  std::string exe_dir = GetExecutableDir();
  std::string filename = exe_dir + "/test-helper";
  std::string linkname = exe_dir + mp::path::preferred_separator + "юникод";
#ifdef _WIN32
  filename += ".exe";
  linkname += ".exe";
#endif
  LinkFile(filename, linkname);
  ExecuteShellCommand(linkname + " > out", false);
  string path = ReadFile("out");
  string ending = linkname + "\n";
  EXPECT_EQ(ending, path.size() >= ending.size() ?
      path.substr(path.size() - ending.size()) : path);
}

TEST(MemoryMappedFileTest, MapZeroTerminated) {
  const char *content = "some content";
  std::string filename = GetExecutableDir() + "test";
  WriteFile(filename, content);
  File file(filename, File::RDONLY);
  MemoryMappedFile<> f(file, file.size());
  EXPECT_STREQ(content, f.start());
  EXPECT_EQ(std::strlen(content), f.size());
}

TEST(MemoryMappedFileTest, DtorUnmapsFile) {
  std::string filename = GetExecutableDir() + "test";
  WriteFile(filename, "abc");
  const volatile char *start = 0;
  {
    File file(filename, File::RDONLY);
    MemoryMappedFile<> f(file, file.size());
    start = f.start();
  }
  EXPECT_DEATH((void)*start, "");
}

#ifndef _WIN32
# ifdef HAVE_LSOF
TEST(MemoryMappedFileTest, CloseFile) {
  std::string path =
      mp::GetExecutablePath().remove_filename().string() + "/test";
  WriteFile(path, "abc");
  File file(path, File::RDONLY);
  MemoryMappedFile<> f(file, file.size());
  file.close();
  int exit_code = ExecuteShellCommand("lsof " + path + " > out", false);
  std::string out = ReadFile("out");
  std::vector<string> results = Split(out, '\n');
#  ifdef __APPLE__
  // For some reason lsof prints txt instead of mem for a mapped file on Mac.
  const char MEM[] = " txt ";
#  else
  const char MEM[] = " mem ";
#  endif
  // Check that lsof prints mem instead of a file descriptor.
  if (results.size() == 3 &&
    results[1].find(MEM) != string::npos && results[2] == "") {
    // Got expected output - do nothing.
  } else {
    // Running lsof with filename failed, so check the full lsof output instead.
    ExecuteShellCommand("lsof > out", false);
    std::string line;
    std::ifstream out("out");
    bool found = false;
    while (std::getline(out, line)) {
      if (line.find(path) != std::string::npos) {
        if (line.find(MEM) != string::npos)
          found = true;
        break;
      }
    }
    if (!found) {
      FAIL() << "Unexpected output from lsof:\n" << out
             << "\nExit code = " << exit_code << "\nPath = "<< path;
    }
  }
  EXPECT_STREQ("abc", f.start());
}
# endif

class StandardErrorHandling {};

#else
TEST(MemoryMappedFileTest, CloseFile) {
  WriteFile("test", "abc");
  DWORD handle_count_before = 0;
  ASSERT_TRUE(GetProcessHandleCount(
      GetCurrentProcess(), &handle_count_before) != 0);
  File file("test", File::RDONLY);
  MemoryMappedFile<> f(file, file.size());
  file.close();
  DWORD handle_count_after = 0;
  ASSERT_TRUE(GetProcessHandleCount(
      GetCurrentProcess(), &handle_count_after) != 0);
  EXPECT_EQ(handle_count_before, handle_count_after);
  EXPECT_STREQ("abc", f.start());
}

// Enables standard handling of invalid arguments in the current scope by
// setting _set_invalid_parameter_handler in ctor and restoring the old
// handler in dtor.
// Usage:
//   {
//     StandardErrorHandling seh;
//     // do something dangerous
//   }
class StandardErrorHandling {
 private:
  _invalid_parameter_handler saved_handler_;

  FMT_DISALLOW_COPY_AND_ASSIGN(StandardErrorHandling);

  static void handler(const wchar_t* expression, const wchar_t* function,
                      const wchar_t* file, unsigned int line,
                      uintptr_t pReserved) {}

 public:
  StandardErrorHandling()
    : saved_handler_(_set_invalid_parameter_handler(handler)) {}
  ~StandardErrorHandling() { _set_invalid_parameter_handler(saved_handler_); }
};

#endif

TEST(MemoryMappedFileTest, InvalidFile) {
  StandardErrorHandling seh;
  EXPECT_THROW(MemoryMappedFile<>(fmt::File(), 1), fmt::SystemError);
}
}

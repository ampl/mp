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

#include "solvers/util/error.h"
#include "solvers/util/os.h"
#include "tests/config.h"
#include "tests/util.h"
#include "gtest/gtest.h"

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

using ampl::MemoryMappedFile;

#ifdef _WIN32
using ampl::UTF8ToUTF16;
using ampl::UTF16ToUTF8;
#endif

using std::string;

namespace {

class PathTest : public ::testing::TestWithParam<char> {};

TEST(PathTest, PreferredSeparator) {
#ifdef _WIN32
  EXPECT_EQ('\\', ampl::path::preferred_separator);
#else
  EXPECT_EQ('/', ampl::path::preferred_separator);
#endif
}

TEST_P(PathTest, PathCtor) {
  EXPECT_EQ("", ampl::path().string());
  std::string s = FixPath("/some/path", GetParam());
  ampl::path p(&s[0], &s[0] + s.size());
  EXPECT_EQ(s, p.string());
  EXPECT_EQ(s, ampl::path(s).string());
}

TEST_P(PathTest, RemoveFilename) {
  ampl::path p;
  EXPECT_EQ("", p.remove_filename().string());
  EXPECT_EQ("", p.string());
  std::string s = FixPath("/somewhere/out/in/space", GetParam());
  p = ampl::path(s);
  s = FixPath("/somewhere/out/in", GetParam());
  EXPECT_EQ(s, p.remove_filename().string());
  EXPECT_EQ(s, p.string());
  s = FixPath("/somewhere/out", GetParam());
  EXPECT_EQ(s, p.remove_filename().string());
  EXPECT_EQ(s, p.string());
  s = FixPath("/", GetParam());
  EXPECT_EQ(s, ampl::path(s).remove_filename().string());
  EXPECT_EQ("", ampl::path("test").remove_filename().string());
}

TEST_P(PathTest, Filename) {
  ampl::path p;
  EXPECT_EQ("", p.filename().string());
  EXPECT_EQ("", p.string());
  std::string s = FixPath("/somewhere/out/in/space", GetParam());
  p = ampl::path(s);
  EXPECT_EQ("space", p.filename().string());
  EXPECT_EQ(s, p.string());
  EXPECT_EQ("", ampl::path(FixPath("/", GetParam())).filename().string());
  EXPECT_EQ("test", ampl::path("test").filename().string());
}

INSTANTIATE_TEST_CASE_P(POSIX, PathTest, ::testing::Values('/'));

#ifdef _WIN32
INSTANTIATE_TEST_CASE_P(Win32, PathTest, ::testing::Values('\\'));
#endif

TEST(PathTest, TempDirectoryPath) {
#ifndef _WIN32
  const char *dir = std::getenv("TMPDIR");
  EXPECT_EQ(dir ? dir : "/tmp", ampl::path::temp_directory_path().string());
#else
  wchar_t buffer[MAX_PATH + 1];
  DWORD result = GetTempPath(MAX_PATH + 1, buffer);
  EXPECT_GT(result, 0);
  EXPECT_LE(result, MAX_PATH);
  EXPECT_STREQ(UTF16ToUTF8(buffer), ampl::path::temp_directory_path().string());
#endif
}

TEST(OSTest, GetExecutablePath) {
  string path = ampl::GetExecutablePath().string();
  string ending = FixBinaryPath("/util/os-test");
#ifdef _WIN32
  ending += ".exe";
#endif
  EXPECT_EQ(ending, path.size() >= ending.size() ?
      path.substr(path.size() - ending.size()) : path);
}

// Creates a new link for a file or copies the file if the system doesn't
// support symlinks. Both filenames are UTF-8 encoded.
void LinkFile(fmt::StringRef filename, fmt::StringRef linkname) {
#ifndef _WIN32
  int result = link(filename.c_str(), linkname.c_str());
  if (result && errno != EEXIST) {
    ampl::ThrowSystemError(errno, "cannot create a symlink from {} to {}")
      << filename.c_str() << linkname.c_str();
  }
#else
  if (!CopyFileW(UTF8ToUTF16(filename), UTF8ToUTF16(linkname), FALSE)) {
    ampl::ThrowSystemError(GetLastError(), "cannot copy file {} to {}")
      << filename.c_str() << linkname.c_str();
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

TEST(OSTest, GetExecutablePathUnicode) {
  // Neither CMake nor NMake handle Unicode paths properly on Windows,
  // so copy test executable ourselves.
  std::string filename = FixBinaryPath("test-helper");
  std::string linkname = FixBinaryPath("юникод");
#ifdef _WIN32
  filename += ".exe";
  linkname += ".exe";
#endif
  LinkFile(filename, linkname);
  ExecuteShellCommand("./" + linkname + " > out", false);
  string path = ReadFile("out");
  string ending = "/" + linkname;
  EXPECT_EQ(ending, path.size() >= ending.size() ?
      path.substr(path.size() - ending.size()) : path);
}

#ifdef _WIN32
TEST(OSTest, UTF16ToUTF8) {
  std::string s = "ёжик";
  UTF16ToUTF8 u(L"\x0451\x0436\x0438\x043A");
  EXPECT_STREQ(s.c_str(), u);
  EXPECT_EQ(s.size(), u.size());
}

TEST(OSTest, UTF8ToUTF16) {
  std::string s = "лошадка";
  UTF8ToUTF16 u(s.c_str());
  EXPECT_STREQ(L"\x043B\x043E\x0448\x0430\x0434\x043A\x0430", u);
  EXPECT_EQ(7, u.size());
}
#endif  // _WIN32

TEST(MemoryMappedFileTest, MapZeroTerminated) {
  const char *content = "some content";
  WriteFile("test", content);
  MemoryMappedFile f("test");
  EXPECT_STREQ(content, f.start());
  EXPECT_EQ(std::strlen(content), f.size());
}

TEST(MemoryMappedFileTest, DtorUnmapsFile) {
  WriteFile("test", "abc");
  const volatile char *start = 0;
  {
    MemoryMappedFile f("test");
    start = f.start();
  }
  EXPECT_DEATH((void)*start, "");
}

#ifndef _WIN32
# ifdef HAVE_LSOF
TEST(MemoryMappedFileTest, CloseFile) {
  WriteFile("test", "abc");
  MemoryMappedFile f("test");
  std::string path =
    ampl::GetExecutablePath().remove_filename().string() + "/test";
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
}
# endif
#else
TEST(MemoryMappedFileTest, CloseFile) {
  WriteFile("test", "abc");
  DWORD handle_count_before = 0;
  ASSERT_TRUE(GetProcessHandleCount(
      GetCurrentProcess(), &handle_count_before) != 0);
  MemoryMappedFile f("test");
  DWORD handle_count_after = 0;
  ASSERT_TRUE(GetProcessHandleCount(
      GetCurrentProcess(), &handle_count_after) != 0);
  EXPECT_EQ(handle_count_before, handle_count_after);
  EXPECT_STREQ("abc", f.start());
}
#endif

TEST(MemoryMappedFileTest, NonexistentFile) {
  EXPECT_THROW(MemoryMappedFile("nonexistent"), ampl::SystemError);
}
}

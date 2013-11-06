/*
 Filesystem library tests.

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
#include "tests/util.h"
#include "gtest/gtest.h"

#include <cstring>
#include <string>
#include <vector>

using ampl::MemoryMappedFile;

using std::string;

namespace {

TEST(OSTest, EmptyPath) {
  ampl::path p;
  EXPECT_EQ("", p.string());
  EXPECT_EQ("", p.remove_filename().string());
  EXPECT_EQ("", p.string());
}

TEST(OSTest, NonemptyPath) {
  const char *s = "/somewhere/out/in/space";
  ampl::path p(s, s + std::strlen(s));
  EXPECT_EQ("/somewhere/out/in/space", p.string());
  EXPECT_EQ("/somewhere/out/in", p.remove_filename().string());
  EXPECT_EQ("/somewhere/out/in", p.string());
  EXPECT_EQ("/somewhere/out", p.remove_filename().string());
  EXPECT_EQ("/somewhere/out", p.string());
}

TEST(OSTest, GetExecutablePath) {
  string path = ampl::GetExecutablePath().string();
  string ending = "/util/os_test";
#ifdef _WIN32
  ending += ".exe";
#endif
  EXPECT_EQ(ending, path.size() >= ending.size() ?
      path.substr(path.size() - ending.size()) : path);
}

TEST(FilesystemTest, GetExecutablePathUnicode) {
  ExecuteShellCommand("./юникод > out");
  string path = ReadFile("out");
  string ending = "/юникод";
#ifdef _WIN32
  ending += ".exe";
#endif
  EXPECT_EQ(ending, path.size() >= ending.size() ?
      path.substr(path.size() - ending.size()) : path);
}

#ifdef _WIN32
TEST(OSTest, UTF16ToUTF8) {
  std::string s = "ёжик";
  ampl::UTF16ToUTF8 u(L"\x0451\x0436\x0438\x043A");
  EXPECT_EQ(s, u);
  EXPECT_EQ(s.size(), u.size());
}

TEST(OSTest, UTF8ToUTF16) {
  std::string s = "лошадка";
  ampl::UTF8ToUTF16 u(s.c_str());
  EXPECT_STREQ(L"\x043B\x043E\x0448\x0430\x0434\x043A\x0430", u);
  EXPECT_EQ(s.size(), u.size());
}
#endif  // _WIN32

TEST(MemoryMappedFileTest, MapZeroTerminated) {
  WriteFile("test", "some content");
  MemoryMappedFile f("test");
  EXPECT_STREQ("some content", f.start());
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

TEST(MemoryMappedFileTest, CloseFile) {
  WriteFile("test", "abc");
  MemoryMappedFile f("test");
#ifndef _WIN32
  ExecuteShellCommand("lsof test > out");
  std::string out = ReadFile("out");
  std::vector<string> results = Split(out, '\n');
#ifdef __APPLE__
  // For some reason lsof prints txt instead of mem for a mapped file on Mac.
  const char MEM[] = " txt ";
#else
  const char MEM[] = " mem ";
#endif
  // Check that lsof prints mem instead of a file descriptor.
  EXPECT_TRUE(results.size() == 3 &&
      results[1].find(MEM) != string::npos && results[2] == "")
    << "Unexpected output from lsof:\n" << out;
#else
  // TODO: windows test
#endif
}

TEST(MemoryMappedFileTest, NonexistentFile) {
  EXPECT_THROW(MemoryMappedFile("nonexistent"), ampl::SystemError);
}
}

// TODO: test UTF16ToUTF8 and UTF8ToUTF16

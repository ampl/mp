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
#include "solvers/util/filesystem.h"
#include "tests/util.h"
#include "gtest/gtest.h"

#include <cstring>
#include <string>
#include <vector>

using ampl::MemoryMappedFile;

using std::string;

namespace {

TEST(FilesystemTest, EmptyPath) {
  ampl::path p;
  EXPECT_EQ("", p.string());
  EXPECT_EQ("", p.remove_filename().string());
  EXPECT_EQ("", p.string());
}

TEST(FilesystemTest, NonemptyPath) {
  const char *s = "/somewhere/out/in/space";
  ampl::path p(s, s + std::strlen(s));
  EXPECT_EQ("/somewhere/out/in/space", p.string());
  EXPECT_EQ("/somewhere/out/in", p.remove_filename().string());
  EXPECT_EQ("/somewhere/out/in", p.string());
  EXPECT_EQ("/somewhere/out", p.remove_filename().string());
  EXPECT_EQ("/somewhere/out", p.string());
}

TEST(FilesystemTest, GetExecutablePath) {
  string path = ampl::GetExecutablePath().string();
  string ending = "/util/filesystem_test";
#ifdef WIN32
  ending += ".exe";
#endif
  EXPECT_EQ(ending, path.size() >= ending.size() ?
      path.substr(path.size() - ending.size()) : path);
}

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
#ifndef WIN32
  ExecuteShellCommand("lsof test > out");
  std::string out = ReadFile("out");
  std::vector<string> results = Split(out, '\n');
  ASSERT_EQ(3, results.size()) << "Unexpected output from lsof:\n" << out;
  // Check that lsof prints mem instead of a file descriptor.
  EXPECT_TRUE(results[1].find(" mem ") != string::npos);
  EXPECT_EQ("", results[2]);
#endif
}

TEST(MemoryMappedFileTest, NonexistentFile) {
  EXPECT_THROW(MemoryMappedFile("nonexistent"), ampl::SystemError);
}

// TODO: test Unicode path in GetExecutablePath
}

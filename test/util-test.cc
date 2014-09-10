/*
 Tests of test utilities

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

#include "gtest/gtest.h"
#include "mp/error.h"
#include "mp/os.h"
#include "util.h"
#include "gtest-extra.h"
#include <string>

namespace {

TEST(UtilTest, Split) {
  std::vector<std::string> result = Split("abc", ' ');
  EXPECT_EQ(1u, result.size());
  EXPECT_EQ("abc", result[0]);
  result = Split("a b c", ' ');
  EXPECT_EQ(3u, result.size());
  EXPECT_EQ("a", result[0]);
  EXPECT_EQ("b", result[1]);
  EXPECT_EQ("c", result[2]);
  result = Split("abc ", ' ');
  EXPECT_EQ(2u, result.size());
  EXPECT_EQ("abc", result[0]);
  EXPECT_EQ("", result[1]);
}

TEST(UtilTest, FixPath) {
  std::string path = "/somewhere/out/in/space";
  EXPECT_EQ("-somewhere-out-in-space", FixPath(path, '-'));
  std::replace(path.begin(), path.end(), '/', mp::path::preferred_separator);
  EXPECT_EQ(path, FixPath("/somewhere/out/in/space"));
}

TEST(UtilTest, ReplaceLine) {
  EXPECT_EQ("de", ReplaceLine("", 0, "de"));
  EXPECT_EQ("de", ReplaceLine("abc", 0, "de"));
  EXPECT_THROW(ReplaceLine("abc", 1, "de"), mp::Error);
  EXPECT_EQ("de\n", ReplaceLine("abc\n", 0, "de"));
  EXPECT_EQ("abc\nde", ReplaceLine("abc\n", 1, "de"));
  EXPECT_EQ("gh\ndef", ReplaceLine("abc\ndef", 0, "gh"));
  EXPECT_EQ("abc\ngh", ReplaceLine("abc\ndef", 1, "gh"));
  EXPECT_THROW(ReplaceLine("abc\ndef", 2, "gh"), mp::Error);
}

TEST(UtilTest, ExecuteShellCommand) {
  ExecuteShellCommand("cd .");
  std::string message;
  try {
    ExecuteShellCommand("dir/bad-command");
  } catch (const mp::Error &e) {
    message = e.what();
  }
  EXPECT_TRUE(message.find("process exited with code ") != std::string::npos)
    << message;
}

TEST(UtilTest, ExecuteShellThrowsOnNonzeroExitCode) {
  EXPECT_THROW_MSG(
    ExecuteShellCommand(GetExecutableDir() + "/test-helper > out"),
    mp::Error, "process exited with code 42");
}
}

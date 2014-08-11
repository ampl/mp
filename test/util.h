/*
 Test utilities

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

#ifndef TEST_UTIL_H_
#define TEST_UTIL_H_

#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>

#include "mp/format.h"
#include "mp/os.h"

std::string ReadFile(fmt::StringRef name);
void WriteFile(fmt::StringRef name, fmt::StringRef data);

// Replaces all occurrences of '/' in the path with sep. If sep is 0
// it is set to the system-specific directory separator.
std::string FixPath(
    fmt::StringRef path, char sep = mp::path::preferred_separator);

// Changes the current working directory. Throws Error on error.
void ChangeDirectory(fmt::StringRef path);

// Executes a shell command. Throws Error on error.
int ExecuteShellCommand(
  fmt::StringRef command, bool throw_on_nonzero_exit_code = true);

// Fix the path to a binary (executable or shared library) file by
// inserting a configuration directory (Debug or Release) which is
// deduced from the path of the currently running executable.
std::string FixBinaryPath(fmt::StringRef path);

// Splits the string into an array of substrings.
std::vector<std::string> Split(const std::string &s, char sep);

// Replace line at line_index in s with new_line.
std::string ReplaceLine(std::string s, int line_index, const char *new_line);

#define FORMAT_TEST_THROW_(statement, expected_exception, message, fail) \
  GTEST_AMBIGUOUS_ELSE_BLOCKER_ \
  if (::testing::internal::ConstCharPtr gtest_msg = "") { \
    bool gtest_caught_expected = false; \
    std::string gtest_actual_message; \
    try { \
      GTEST_SUPPRESS_UNREACHABLE_CODE_WARNING_BELOW_(statement); \
    } \
    catch (expected_exception const& e) { \
      gtest_caught_expected = true; \
      gtest_actual_message = e.what(); \
    } \
    catch (...) { \
      gtest_msg.value = \
          "Expected: " #statement " throws an exception of type " \
          #expected_exception ".\n  Actual: it throws a different type."; \
      goto GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__); \
    } \
    if (!gtest_caught_expected) { \
      gtest_msg.value = \
          "Expected: " #statement " throws an exception of type " \
          #expected_exception ".\n  Actual: it throws nothing."; \
      goto GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__); \
    } else { \
      EXPECT_EQ(message, gtest_actual_message); \
    } \
  } else \
    GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__): \
      fail(gtest_msg.value)

#define EXPECT_THROW_MSG(statement, expected_exception, expected_message) \
  FORMAT_TEST_THROW_(statement, expected_exception, expected_message, \
      GTEST_NONFATAL_FAILURE_)

#endif  // TEST_UTIL_H_

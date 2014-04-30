/*
 Tests of error classes and functions.

 Copyright (C) 2012 AMPL Optimization Inc

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
#include "solvers/util/error.h"
#include "solvers/util/os.h"

#include <errno.h>
#include <string.h>

#ifdef _WIN32
# include <windows.h>
#endif

namespace {

TEST(ErrorTest, Error) {
  ampl::Error e(fmt::StringRef("test"));
  EXPECT_STREQ("test", e.what());
}

TEST(ErrorTest, ThrowError) {
  ampl::Error error("");
  try {
    ampl::ThrowError("test {} message") << "error";
  } catch (const ampl::Error &e) {
    error = e;
  }
  EXPECT_STREQ("test error message", error.what());
}

TEST(ErrorTest, ReportSystemError) {
  const int TEST_ERROR = EDOM;
  EXPECT_EXIT({
    ampl::ReportSystemError(TEST_ERROR, "test error");
    std::fprintf(stderr, "end\n");
    std::exit(0);
  }, ::testing::ExitedWithCode(0),
      str(fmt::Format("test error: {}\nend\n") << strerror(TEST_ERROR)));
}

#ifdef _WIN32
TEST(ErrorTest, ReportWinError) {
  const int TEST_ERROR = ERROR_FILE_EXISTS;
  LPWSTR message = 0;
  FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER |
      FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, 0,
      ERROR_FILE_EXISTS, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
      reinterpret_cast<LPWSTR>(&message), 0, 0);
  fmt::internal::UTF16ToUTF8 utf8_message(message);
  LocalFree(message);
  EXPECT_EXIT({
    ampl::ReportSystemError(TEST_ERROR, "test error");
    std::fprintf(stderr, "end\n");
    std::exit(0);
  }, ::testing::ExitedWithCode(0),
      str(fmt::Format("test error: {}\nend\n") << fmt::str(utf8_message)));
}
#endif
}

/*
 Tests of the AMPL table proxy.

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
#include "util/format.h"
#include "tests/function.h"
#include "tests/util.h"

namespace {

class TableProxyTest : public ::testing::Test {
 protected:
  static fun::Library lib_;
  const fun::Handler *handler_;
  std::vector<std::string> strings_;

  static void SetUpTestCase() {
    lib_.Load();
  }

  static void TearDownTestCase() {
    lib_.Unload();
  }

  void SetUp() {
    handler_ = lib_.GetHandler("tableproxy");
    strings_.push_back("tableproxy");
    int bits = sizeof(void*) == 8 ? 64 : 32;
    std::string prog = FixBinaryPath(
        fmt::format("prog=../../bin/tableproxy{}", bits));
    strings_.push_back(prog);
    strings_.push_back(FixBinaryPath("lib=../../tables/fullbit.dll"));
    strings_.push_back("lib-tab");
  }
};

fun::Library TableProxyTest::lib_(FixBinaryPath("../../tables/ampltabl.dll"));

TEST_F(TableProxyTest, WriteTab) {
  fun::Table t("test", 1, 0, strings_);
  t = "N", 1, 2;
  std::remove("test.tab");
  handler_->Write(t);
  EXPECT_EQ("ampl.tab 1 0\nN\n1\n2\n", ReadFile("test.tab"));
}

TEST_F(TableProxyTest, WriteDuplicateStringNumKey) {
  fun::Table t("test", 2, 0, strings_);
  t = "S", "N",
      "a",  1;
  std::remove("test.tab");
  handler_->Write(t, fun::Handler::INOUT);
  handler_->Write(t, fun::Handler::INOUT);
  EXPECT_EQ("ampl.tab 2 0\nS\tN\na\t1\n", ReadFile("test.tab"));
}
}

/*
 Tests of the CSV ODBC connection.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "gtest/gtest.h"
#include "tests/function.h"
#include "tests/odbc.h"

using fun::Handler;
using fun::Table;

namespace {

class CSVTest : public ::testing::Test {
 protected:
  static fun::Library lib_;
  const Handler *handler_;

  static void SetUpTestCase() {
    lib_.Load();
  }

  void SetUp() {
    handler_ = lib_.GetHandler("odbc");
  }
};

fun::Library CSVTest::lib_("../tables/ampltabl.dll");

TEST_F(CSVTest, Read) {
  Table t("test", 1);
  t = "S";
  odbc::Env env;
  handler_->Read("DRIVER={" + env.FindDriver("*.csv") + "}; DBQ=data",
      &t, "SELECT * FROM test.csv");
  EXPECT_EQ(2u, t.num_rows());
  EXPECT_STREQ("abc", t(0, 0).string());
  EXPECT_STREQ("def", t(1, 0).string());
}
}

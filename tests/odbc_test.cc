/*
 Tests of the AMPL ODBC table handler.

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

#include <cctype>
#include <cstring>
#include <algorithm>
#include <deque>
#include <vector>

#undef VOID

#include "gtest/gtest.h"
#include "tests/function.h"
#include "tests/odbc.h"
#include "solvers/asl.h"
#include "tests/config.h"

#define SERVER "callisto.local"

using fun::Table;

namespace {

class ODBCTest : public ::testing::Test {
 protected:
  static fun::Library lib_;

  static void SetUpTestCase() {
    lib_.Load();
  }

  void Read(Table *t) {
    lib_.GetHandler("odbc")->Read(t);
  }
};

fun::Library ODBCTest::lib_("../tables/ampltabl.dll");

TEST_F(ODBCTest, Loaded) {
  EXPECT_EQ("", lib_.error());
  EXPECT_TRUE(lib_.GetHandler("odbc"));
}

TEST_F(ODBCTest, ReadMySQL) {
  std::string driver_name(odbc::Env().FindDriver("mysql"));
  if (driver_name.empty()) {
    std::cerr << "Skipping MySQL test." << std::endl;
    return;
  }
  std::string connection(
      "DRIVER={" + driver_name + "}; SERVER=" SERVER "; DATABASE=test;");
  Table t("", "ODBC", connection.c_str(), "SQL=SELECT VERSION();");
  t.AddCol("VERSION()");
  Read(&t);
  EXPECT_EQ(1, t.num_rows());
  EXPECT_TRUE(t.GetString(0) != nullptr);
}

// TODO(viz): more tests
}

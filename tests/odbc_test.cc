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

#include <cstring>
#include <algorithm>
#include <deque>
#include <vector>

#include "gtest/gtest.h"
#include "tests/function.h"
#include "solvers/asl.h"
#include "tests/config.h"

using fun::Table;

namespace {

class ODBCTest : public ::testing::Test {
 protected:
  static fun::Library lib_;

  static void SetUpTestCase() {
    lib_.Load();
  }
};

fun::Library ODBCTest::lib_("../tables/ampltabl.dll");

TEST_F(ODBCTest, Loaded) {
  EXPECT_EQ("", lib_.error());
  EXPECT_TRUE(lib_.GetHandler("odbc"));
}

TEST_F(ODBCTest, ReadMySQL) {
  const fun::Handler *handler = lib_.GetHandler("odbc");
  Table t("", "ODBC",
      "DRIVER={MySQL ODBC 5.2 Driver}; SOCKET={/var/run/mysqld/mysqld.sock};",
      "SQL=SELECT VERSION();");
  t.AddCol("VERSION()");
  handler->Read(t);
  EXPECT_EQ(1, t.num_rows());
  EXPECT_TRUE(t.GetString(0) != nullptr);
}

// TODO: more tests
}

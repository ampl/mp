/*
 Tests of the Access ODBC connection.

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
#include <climits>

#include "gtest/gtest.h"
#include "solvers/util/format.h"
#include "tests/function.h"
#include "tests/odbc.h"

using fun::Handler;
using fun::Table;

namespace {

class AccessTest : public ::testing::Test {
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

fun::Library AccessTest::lib_("../tables/ampltabl.dll");

TEST_F(AccessTest, Read) {
  Table t("Test", 1);
  t = "N";
  handler_->Read("data/test.accdb", &t);
  EXPECT_EQ(1u, t.num_rows());
  EXPECT_EQ(42, t(0, 0).number());
}

TEST_F(AccessTest, ReadNullField) {
  Table t("TableWithNullField", 1);
  t = "NullField";
  handler_->Read("data/test.accdb", &t);
  EXPECT_EQ(0u, t.num_rows());
}
}

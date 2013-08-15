/*
 Tests of the Excel ODBC connection.

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

class ExcelTest : public ::testing::Test {
 protected:
  static fun::Library lib_;
  const Handler *handler_;
  std::vector<std::string> strings_;

  static void SetUpTestCase() {
    lib_.Load();
  }

  void SetUp() {
    handler_ = lib_.GetHandler("odbc");
    strings_.push_back("odbc");
    strings_.push_back(
        "DRIVER={Microsoft Excel Driver (*.xls, *.xlsx, *.xlsm, *.xlsb)}; "
        "DBQ=data/test.xls; READONLY=FALSE;");
  }
};

fun::Library ExcelTest::lib_("../tables/ampltabl.dll");

TEST_F(ExcelTest, Read) {
  Table t("TableWith256CharCell", 0, 1, strings_);
  t = "S";
  handler_->Read(&t);
  EXPECT_EQ(1u, t.num_rows());
  EXPECT_STREQ(
      "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghij"
      "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghij"
      "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghij"
      "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghij"
      "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghij"
      "abcdef", t(0, 0).string());
}

// A DSN for "Microsoft Excel Driver (*.xls)" with "Read Only" checkbox
// unchecked may have to be created in ODBC Data Source Administrator
// for this test to pass.
TEST_F(ExcelTest, WriteMaxColumnsExcel2003) {
  // Excel 2003 and earlier support at most 256 columns
  // http://office.microsoft.com/en-us/excel-help/excel-specifications-and-limits-HP005199291.aspx
  int num_cols = 256;
  Table t("TableWithManyCols", 0, num_cols, strings_);
  for (int i = 1; i <= num_cols; ++i)
    t.Add(c_str(fmt::Format("c{}") << i));
  handler_->Write(t);
}

TEST_F(ExcelTest, AppendZeroRows) {
  Table t("AppendZeroRowsTable", 1, 0, strings_);
  t = "S",
       1;
  handler_->Write(t);
  t = "S";
  handler_->Write(t, Handler::INOUT);
}

#ifndef _WIN64
// The test fails on 64-bit Windows XP with 64-bit ODBC drivers.
// According to http://www.microsoft.com/en-us/download/details.aspx?id=13255
// only the 32-bit Access Database Engine (which includes the Excel driver) may
// be used on Windows XP Service Pack 3.
TEST_F(ExcelTest, WriteMaxColumnsExcel2007) {
  // Excel 2007 supports up to 16384 columns, but the ODBC driver only
  // allows 256.
  int num_cols = 256;
  strings_[1] =
      "DRIVER={Microsoft Excel Driver (*.xls, *.xlsx, *.xlsm, *.xlsb)}; "
      "DBQ=data/test.xlsx; READONLY=FALSE;";
  Table t("TableWithManyCols", 0, num_cols, strings_);
  for (int i = 1; i <= num_cols; ++i)
    t.Add(c_str(fmt::Format("c{}") << i));
  handler_->Write(t);
}
#endif
}

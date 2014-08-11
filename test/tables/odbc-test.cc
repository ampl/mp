/*
 Tests of the AMPL ODBC table handler.

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
#include "tests/function.h"
#include "tests/config.h"
#include "tests/tables/odbc.h"
#include "tests/util.h"

namespace {

TEST(ODBCTest, Load) {
  fun::Library lib(FixBinaryPath("../../bin/ampltabl.dll"));
  lib.Load();
  EXPECT_EQ("", lib.error());
  EXPECT_TRUE(lib.GetHandler("odbc") != nullptr);
}
}

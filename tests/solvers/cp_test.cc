/*
 CP function library tests.

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

using fun::Function;
using fun::FunctionInfo;
using fun::MakeArgs;

namespace {

class CPTest : public ::testing::Test {
 protected:
  static const FunctionInfo info;  // Default function info.
  static fun::Library lib_;

  static void SetUpTestCase() { lib_.Load(); }

  // Returns an AMPL function by name.
  static Function GetFunction(const char *name, const FunctionInfo &info) {
    const func_info *fi = lib_.GetFunction(name);
    if (!fi)
      throw std::runtime_error(std::string("function not found: ") + name);
    return Function(&lib_, fi, &info);
  }

  static Function GetFunction(const char *name) {
    return GetFunction(name, info);
  }
};

const FunctionInfo CPTest::info;
fun::Library CPTest::lib_("../../solvers/cp/cp.dll");

TEST_F(CPTest, NoLibErrors) {
  EXPECT_EQ("", lib_.error());
}

TEST_F(CPTest, Element) {
  Function element = GetFunction("element");
  EXPECT_THROW(element(fun::Tuple()), std::invalid_argument);
  EXPECT_THROW(element(42), std::invalid_argument);
  EXPECT_EQ(42, element(MakeArgs(42, 0)).value());
  EXPECT_EQ(11, element(MakeArgs(11, 22, 33, 0)).value());
  EXPECT_EQ(22, element(MakeArgs(11, 22, 33, 1)).value());
  EXPECT_EQ(33, element(MakeArgs(11, 22, 33, 2)).value());
  EXPECT_STREQ("invalid index -1", element(MakeArgs(11, 22, 33, -1)).error());
  EXPECT_STREQ("invalid index 3", element(MakeArgs(11, 22, 33, 3)).error());
  EXPECT_STREQ("derivatives are not provided",
      element(MakeArgs(42, 0), fun::DERIVS).error());
  EXPECT_STREQ("derivatives are not provided",
      element(MakeArgs(42, 0), fun::HES).error());
}

TEST_F(CPTest, InRelation) {
  Function in_relation = GetFunction("in_relation");
  EXPECT_THROW(in_relation(fun::Tuple()), std::invalid_argument);
  EXPECT_STREQ("can't evaluate in_relation", in_relation(MakeArgs(0)).error());
}
}

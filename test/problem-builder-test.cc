/*
 A problem builder reader tests.

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

#include "gtest-extra.h"
#include "mp/problem-builder.h"
#include "mp/nl.h"

class TestExpr {};

class TestProblemBuilder :
    public mp::ProblemBuilder<TestProblemBuilder, TestExpr> {};

// Test that ProblemBuilder can be used with ProblemBuilderToNLAdapter.
TEST(ProblemBuilderTest, UseWithProblemBuilderToNLAdapter) {
  TestProblemBuilder builder;
  mp::ProblemBuilderToNLAdapter<TestProblemBuilder> handler(builder);
  EXPECT_THROW_MSG(handler.OnNumericConstant(0), mp::Error,
                   "unsupported: numeric constant in nonlinear expression");
}

#define EXPECT_ERROR(call, construct) { \
  TestProblemBuilder builder; \
  EXPECT_THROW_MSG(builder.call, mp::Error, \
                   fmt::format("unsupported: {}", construct)); \
}

TEST(ProblemBuilderTest, ReportUnhandled) {
  EXPECT_ERROR(SetObj(0, mp::obj::MIN, TestExpr()), "objective");
  EXPECT_ERROR(SetCon(0, TestExpr()), "nonlinear constraint");
  EXPECT_ERROR(SetLogicalCon(0, TestExpr()), "logical constraint");
  EXPECT_ERROR(SetVar(0, TestExpr(), 0), "nonlinear defined variable");
  EXPECT_ERROR(SetComplement(0, 0, 0), "complementarity constraint");
  // TODO
}

// TODO: more tests

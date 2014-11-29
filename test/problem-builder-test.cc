/*
 Problem builder tests.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "gmock/gmock.h"
#include "gtest-extra.h"
#include "mp/problem-builder.h"
#include "mp/nl.h"

class TestExpr {};

struct TestProblemBuilder : mp::ProblemBuilder<TestProblemBuilder, TestExpr> {
  MOCK_METHOD1(ReportUnhandledConstruct, void (const std::string &name));
};

// Test that ProblemBuilder can be used with ProblemBuilderToNLAdapter.
TEST(ProblemBuilderTest, UseWithProblemBuilderToNLAdapter) {
  TestProblemBuilder builder;
  mp::ProblemBuilderToNLAdapter<TestProblemBuilder> handler(builder);
  EXPECT_CALL(builder, ReportUnhandledConstruct(::testing::_)).
      Times(::testing::Exactly(2));
  handler.OnNumericConstant(0);
  auto expr_builder = handler.BeginCommonExpr(0, 0);
  handler.EndCommonExpr(expr_builder, TestExpr(), 0);
}

// Check that handling problem info doesn't throw an exception.
TEST(ProblemBuilderTest, SetInfo) {
  TestProblemBuilder builder;
  builder.SetInfo(mp::ProblemInfo());
}

#define EXPECT_DISPATCH(call, construct) { \
  TestProblemBuilder builder; \
  EXPECT_CALL(builder, ReportUnhandledConstruct(construct)); \
  builder.call; \
}

TEST(ProblemBuilderTest, ReportUnhandledConstruct) {
  EXPECT_DISPATCH(AddVar(0, 0, mp::var::INTEGER), "variable");
  EXPECT_DISPATCH(AddObj(mp::obj::MIN, TestExpr(), 0), "objective");
  EXPECT_DISPATCH(AddCon(TestExpr(), 0, 0, 0), "algebraic constraint");
  EXPECT_DISPATCH(AddCon(TestExpr()), "logical constraint");
  EXPECT_DISPATCH(BeginCommonExpr(0), "common expression");
  EXPECT_DISPATCH(SetComplement(0, 0, 0), "complementarity constraint");

  // Initial (dual) values are not reported as unhandled.
  TestProblemBuilder builder;
  EXPECT_CALL(builder, ReportUnhandledConstruct("")).Times(0);
  builder.SetInitialValue(0, 0);
  builder.SetInitialDualValue(0, 0);

  EXPECT_DISPATCH(GetColumnSizeHandler(), "Jacobian column size");
  EXPECT_DISPATCH(SetFunction(0, "foo", 0, mp::func::NUMERIC), "function");
  EXPECT_DISPATCH(AddSuffix(0, 0, "foo"), "suffix");
  EXPECT_DISPATCH(MakeNumericConstant(0),
                  "numeric constant in nonlinear expression");
  EXPECT_DISPATCH(MakeVariable(0), "variable in nonlinear expression");
  EXPECT_DISPATCH(MakeUnary(mp::expr::ABS, TestExpr()), "abs");
  EXPECT_DISPATCH(MakeBinary(mp::expr::ADD, TestExpr(), TestExpr()), "+");
  EXPECT_DISPATCH(MakeIf(TestExpr(), TestExpr(), TestExpr()), "if expression");
  EXPECT_DISPATCH(BeginPLTerm(0), "piecewise-linear term");
  EXPECT_DISPATCH(EndPLTerm(TestProblemBuilder::PLTermHandler(), TestExpr()),
                  "piecewise-linear term");
  EXPECT_DISPATCH(BeginCall(0, 0), "function call");
  EXPECT_DISPATCH(EndCall(TestProblemBuilder::CallArgHandler()),
                  "function call");
  EXPECT_DISPATCH(BeginVarArg(mp::expr::MIN, 0), "min");
  EXPECT_DISPATCH(EndVarArg(TestProblemBuilder::NumericArgHandler()),
                  "vararg expression");
  EXPECT_DISPATCH(BeginSum(0), "sum");
  EXPECT_DISPATCH(EndSum(TestProblemBuilder::NumericArgHandler()), "sum");
  EXPECT_DISPATCH(BeginCount(0), "count expression");
  EXPECT_DISPATCH(EndCount(TestProblemBuilder::LogicalArgHandler()),
                  "count expression");
  EXPECT_DISPATCH(BeginNumberOf(0, TestExpr()), "numberof expression");
  EXPECT_DISPATCH(EndNumberOf(TestProblemBuilder::NumericArgHandler()),
                  "numberof expression");
  EXPECT_DISPATCH(MakeLogicalConstant(true), "logical constant");
  EXPECT_DISPATCH(MakeNot(TestExpr()), "logical not");
  EXPECT_DISPATCH(MakeBinaryLogical(mp::expr::OR, TestExpr(), TestExpr()),
                  "||");
  EXPECT_DISPATCH(MakeRelational(mp::expr::LT, TestExpr(), TestExpr()), "<");
  EXPECT_DISPATCH(MakeLogicalCount(mp::expr::ATLEAST, TestExpr(), TestExpr()),
                  "atleast");
  EXPECT_DISPATCH(MakeImplication(TestExpr(), TestExpr(), TestExpr()),
                  "implication expression");
  EXPECT_DISPATCH(BeginIteratedLogical(mp::expr::EXISTS, 0),
                  "exists");
  EXPECT_DISPATCH(EndIteratedLogical(TestProblemBuilder::LogicalArgHandler()),
                  "iterated logical expression");
  EXPECT_DISPATCH(BeginPairwise(mp::expr::ALLDIFF, 0), "alldiff expression");
  EXPECT_DISPATCH(EndPairwise(TestProblemBuilder::PairwiseArgHandler()),
                  "alldiff expression");
  EXPECT_DISPATCH(MakeStringLiteral("test"), "string literal");
}

struct TestProblemBuilder2 :
    mp::ProblemBuilder<TestProblemBuilder2, TestExpr> {};

TEST(ProblemBuilderTest, Throw) {
  TestProblemBuilder2 builder;
  EXPECT_THROW_MSG(builder.MakeVariable(0), mp::Error,
                   "unsupported: variable in nonlinear expression");
}

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
#include "mock-problem-builder.h"
#include "mp/nl-reader.h"

using mp::NLHeader;
namespace expr = mp::expr;

struct TestProblemBuilder : mp::ProblemBuilder<TestProblemBuilder, TestExpr> {
  MOCK_METHOD1(ReportUnhandledConstruct, void (const std::string &name));
};

// Test that ProblemBuilder can be used with NLProblemBuilder.
TEST(ProblemBuilderTest, UseWithNLProblemBuilder) {
  TestProblemBuilder builder;
  mp::internal::NLProblemBuilder<TestProblemBuilder> handler(builder);
  EXPECT_CALL(builder, ReportUnhandledConstruct(::testing::_)).Times(3);
  handler.OnNumber(0);
  handler.BeginCommonExpr(0, 0);
  handler.EndCommonExpr(0, TestExpr(), 0);
}

// Check that handling problem info doesn't throw an exception.
TEST(ProblemBuilderTest, SetInfo) {
  TestProblemBuilder builder;
  builder.SetInfo(mp::NLProblemInfo());
}

#define EXPECT_DISPATCH(call, construct) { \
  TestProblemBuilder builder; \
  EXPECT_CALL(builder, ReportUnhandledConstruct(construct)); \
  builder.call; \
}

TEST(ProblemBuilderTest, ReportUnhandledConstruct) {
  EXPECT_DISPATCH(AddVar(0, 0, mp::var::INTEGER), "variable");
  EXPECT_DISPATCH(AddObj(mp::obj::MIN, TestExpr(), 0), "objective");
  EXPECT_DISPATCH(AddCon(0, 0, TestExpr(), 0), "algebraic constraint");
  EXPECT_DISPATCH(AddCon(TestExpr()), "logical constraint");
  EXPECT_DISPATCH(AddCommonExpr(TestExpr()), "common expression");
  EXPECT_DISPATCH(SetComplementarity(0, 0, mp::ComplInfo(0)),
                  "complementarity constraint");

  // Initial (dual) values are not reported as unhandled.
  TestProblemBuilder builder;
  EXPECT_CALL(builder, ReportUnhandledConstruct("")).Times(0);
  builder.var(0).set_value(0);
  builder.algebraic_con(0).set_dual(0);

  // Suffixes are not reported as unhandled.
  builder.AddIntSuffix("foo", 0, 0);
  builder.AddDblSuffix("foo", 0, 0);

  EXPECT_DISPATCH(AddFunction("foo", 0, mp::func::NUMERIC), "function");
  EXPECT_DISPATCH(MakeNumericConstant(0),
                  "numeric constant in nonlinear expression");
  EXPECT_DISPATCH(MakeVariable(0), "variable in nonlinear expression");
  EXPECT_DISPATCH(MakeUnary(expr::ABS, TestExpr()), "abs");
  EXPECT_DISPATCH(MakeBinary(expr::ADD, TestExpr(), TestExpr()), "+");
  EXPECT_DISPATCH(MakeIf(TestExpr(), TestExpr(), TestExpr()), "if expression");
  EXPECT_DISPATCH(BeginPLTerm(0), "piecewise-linear term");
  EXPECT_DISPATCH(EndPLTerm(TestProblemBuilder::PLTermBuilder(), TestExpr()),
                  "piecewise-linear term");
  EXPECT_DISPATCH(BeginCall(TestProblemBuilder::Function(), 0),
                  "function call");
  EXPECT_DISPATCH(EndCall(TestProblemBuilder::CallExprBuilder()),
                  "function call");
  EXPECT_DISPATCH(BeginIterated(expr::MIN, 0), "min");
  EXPECT_DISPATCH(EndIterated(TestProblemBuilder::NumericExprBuilder()),
                  "vararg expression");
  EXPECT_DISPATCH(BeginSum(0), "sum");
  EXPECT_DISPATCH(EndSum(TestProblemBuilder::NumericExprBuilder()), "sum");
  EXPECT_DISPATCH(BeginCount(0), "count expression");
  EXPECT_DISPATCH(EndCount(TestProblemBuilder::IteratedLogicalExprBuilder()),
                  "count expression");
  EXPECT_DISPATCH(BeginNumberOf(0, TestExpr()), "numberof expression");
  EXPECT_DISPATCH(EndNumberOf(TestProblemBuilder::NumberOfExprBuilder()),
                  "numberof expression");
  EXPECT_DISPATCH(MakeLogicalConstant(true), "logical constant");
  EXPECT_DISPATCH(MakeNot(TestExpr()), "logical not");
  EXPECT_DISPATCH(MakeBinaryLogical(expr::OR, TestExpr(), TestExpr()),
                  "||");
  EXPECT_DISPATCH(MakeRelational(expr::LT, TestExpr(), TestExpr()), "<");
  EXPECT_DISPATCH(MakeLogicalCount(expr::ATLEAST, TestExpr(), TestExpr()),
                  "atleast");
  EXPECT_DISPATCH(MakeImplication(TestExpr(), TestExpr(), TestExpr()),
                  "implication expression");
  EXPECT_DISPATCH(BeginIteratedLogical(expr::EXISTS, 0),
                  "exists");
  EXPECT_DISPATCH(EndIteratedLogical(
                    TestProblemBuilder::IteratedLogicalExprBuilder()),
                  "iterated logical expression");
  EXPECT_DISPATCH(BeginPairwise(expr::ALLDIFF, 0), "alldiff expression");
  EXPECT_DISPATCH(EndPairwise(TestProblemBuilder::PairwiseExprBuilder()),
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

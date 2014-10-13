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

#include "gtest-extra.h"
#include "mp/problem-builder.h"
#include "mp/nl.h"

class TestExpr {};

struct TestProblemBuilder : mp::ProblemBuilder<TestProblemBuilder, TestExpr> {
  std::string name;
  void ReportUnhandledConstruct(fmt::StringRef name) {
    this->name = name.c_str();
  }
};

// Test that ProblemBuilder can be used with ProblemBuilderToNLAdapter.
TEST(ProblemBuilderTest, UseWithProblemBuilderToNLAdapter) {
  TestProblemBuilder builder;
  mp::ProblemBuilderToNLAdapter<TestProblemBuilder> handler(builder);
  handler.OnNumericConstant(0);
  EXPECT_EQ("numeric constant in nonlinear expression", builder.name);
}

// Check that handling problem info doesn't throw an exception.
TEST(ProblemBuilderTest, SetInfo) {
  TestProblemBuilder builder;
  builder.SetInfo(mp::ProblemInfo());
}

#define EXPECT_DISPATCH(call, construct) { \
  TestProblemBuilder builder; \
  builder.call; \
  EXPECT_EQ(construct, builder.name); \
}

TEST(ProblemBuilderTest, ReportUnhandledConstruct) {
  EXPECT_DISPATCH(SetObj(0, mp::obj::MIN, TestExpr()), "objective");
  EXPECT_DISPATCH(SetCon(0, TestExpr()), "nonlinear constraint");
  EXPECT_DISPATCH(SetLogicalCon(0, TestExpr()), "logical constraint");
  EXPECT_DISPATCH(SetCommonExpr(0, TestExpr(), 0),
                  "nonlinear defined variable");
  EXPECT_DISPATCH(SetComplement(0, 0, 0), "complementarity constraint");
  EXPECT_DISPATCH(GetLinearObjBuilder(0, 0), "linear objective");
  EXPECT_DISPATCH(GetLinearConBuilder(0, 0), "linear constraint");
  EXPECT_DISPATCH(GetLinearVarBuilder(0, 0), "linear defined variable");
  EXPECT_DISPATCH(SetVarBounds(0, 0, 0), "variable bound");
  EXPECT_DISPATCH(SetConBounds(0, 0, 0), "constraint bound");
  EXPECT_DISPATCH(SetInitialValue(0, 0), "initial value");
  EXPECT_DISPATCH(SetInitialDualValue(0, 0), "initial dual value");
  EXPECT_DISPATCH(GetColumnSizeHandler(), "Jacobian column size");
  EXPECT_DISPATCH(SetFunction(0, "foo", 0, mp::func::NUMERIC), "function");
  EXPECT_DISPATCH(AddSuffix(0, 0, "foo"), "suffix");
  EXPECT_DISPATCH(MakeNumericConstant(0),
               "numeric constant in nonlinear expression");
  EXPECT_DISPATCH(MakeVariable(0), "variable in nonlinear expression");
  EXPECT_DISPATCH(MakeUnary(mp::expr::ABS, TestExpr()), "unary expression");
  EXPECT_DISPATCH(MakeBinary(mp::expr::ADD, TestExpr(), TestExpr()),
               "binary expression");
  EXPECT_DISPATCH(MakeIf(TestExpr(), TestExpr(), TestExpr()), "if expression");
  EXPECT_DISPATCH(BeginPLTerm(0), "piecewise-linear term");
  EXPECT_DISPATCH(EndPLTerm(TestProblemBuilder::PLTermHandler(), TestExpr()),
               "piecewise-linear term");
  EXPECT_DISPATCH(BeginCall(0, 0), "function call");
  EXPECT_DISPATCH(EndCall(TestProblemBuilder::CallArgHandler()),
                  "function call");
  EXPECT_DISPATCH(BeginVarArg(mp::expr::MIN, 0), "vararg expression");
  EXPECT_DISPATCH(EndVarArg(TestProblemBuilder::NumericArgHandler()),
               "vararg expression");
  EXPECT_DISPATCH(BeginSum(0), "sum");
  EXPECT_DISPATCH(EndSum(TestProblemBuilder::NumericArgHandler()), "sum");
  EXPECT_DISPATCH(BeginCount(0), "count expression");
  EXPECT_DISPATCH(EndCount(TestProblemBuilder::LogicalArgHandler()),
               "count expression");
  EXPECT_DISPATCH(BeginNumberOf(0), "numberof expression");
  EXPECT_DISPATCH(EndNumberOf(TestProblemBuilder::NumericArgHandler()),
               "numberof expression");
  EXPECT_DISPATCH(MakeLogicalConstant(true), "logical constant");
  EXPECT_DISPATCH(MakeNot(TestExpr()), "logical not");
  EXPECT_DISPATCH(MakeBinaryLogical(mp::expr::OR, TestExpr(), TestExpr()),
               "binary logical expression");
  EXPECT_DISPATCH(MakeRelational(mp::expr::LT, TestExpr(), TestExpr()),
               "relational expression");
  EXPECT_DISPATCH(MakeLogicalCount(mp::expr::ATLEAST, TestExpr(), TestExpr()),
               "logical count expression");
  EXPECT_DISPATCH(MakeImplication(TestExpr(), TestExpr(), TestExpr()),
               "implication expression");
  EXPECT_DISPATCH(BeginIteratedLogical(mp::expr::EXISTS, 0),
               "iterated logical expression");
  EXPECT_DISPATCH(EndIteratedLogical(TestProblemBuilder::LogicalArgHandler()),
               "iterated logical expression");
  EXPECT_DISPATCH(BeginAllDiff(0), "alldiff expression");
  EXPECT_DISPATCH(EndAllDiff(TestProblemBuilder::AllDiffArgHandler()),
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

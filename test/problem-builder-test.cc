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
#include "mp/nl.h"

using mp::NLHeader;
namespace expr = mp::expr;

using testing::Field;
using testing::Return;
using testing::StrictMock;

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
  EXPECT_DISPATCH(AddCon(0, 0, TestExpr(), 0), "algebraic constraint");
  EXPECT_DISPATCH(AddCon(TestExpr()), "logical constraint");
  EXPECT_DISPATCH(BeginCommonExpr(0), "common expression");
  EXPECT_DISPATCH(SetComplement(0, 0, 0), "complementarity constraint");

  // Initial (dual) values are not reported as unhandled.
  TestProblemBuilder builder;
  EXPECT_CALL(builder, ReportUnhandledConstruct("")).Times(0);
  builder.SetInitialValue(0, 0);
  builder.SetInitialDualValue(0, 0);

  EXPECT_DISPATCH(GetColumnSizeHandler(), "Jacobian column size");
  EXPECT_DISPATCH(AddFunction("foo", 0, mp::func::NUMERIC), "function");
  EXPECT_DISPATCH(AddIntSuffix("foo", 0, 0), "integer suffix");
  EXPECT_DISPATCH(AddDblSuffix("foo", 0, 0), "double suffix");
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
  EXPECT_DISPATCH(BeginVarArg(expr::MIN, 0), "min");
  EXPECT_DISPATCH(EndVarArg(TestProblemBuilder::NumericExprBuilder()),
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

// Checks if ProblemBuilderToNLAdapter forwards arguments passed to
// adapter_func to ProblemBuilder's builder_func.
#define EXPECT_FORWARD(adapter_func, builder_func, args) { \
  StrictMock<MockProblemBuilder> builder; \
  EXPECT_CALL(builder, builder_func args); \
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder); \
  adapter.adapter_func args; \
}

// Version of EXPECT_FORWARD for methods with a return value.
#define EXPECT_FORWARD_RET(adapter_func, builder_func, args, result) { \
  StrictMock<MockProblemBuilder> builder; \
  EXPECT_CALL(builder, builder_func args).WillOnce(Return(result)); \
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder); \
  adapter.adapter_func args; \
}

TEST(NLProblemBuilderTest, Forward) {
  EXPECT_FORWARD(OnComplement, SetComplement, (66, 77, 88));

  EXPECT_FORWARD(OnInitialValue, SetInitialValue, (33, 4.4));
  EXPECT_FORWARD(OnInitialDualValue, SetInitialDualValue, (55, 6.6));

  EXPECT_FORWARD_RET(OnColumnSizes, GetColumnSizeHandler, (),
                     TestColumnSizeHandler(ID));

  // Use the same StringRef object in arguments, because StringRef objects
  // are compared as pointers and string literals they point to may not
  // be merged.
  fmt::StringRef str("foo");
  EXPECT_FORWARD_RET(OnIntSuffix, AddIntSuffix, (str, 99, 11),
                     TestSuffixHandler<0>(ID));
  EXPECT_FORWARD_RET(OnDblSuffix, AddDblSuffix, (str, 99, 11),
                     TestSuffixHandler<1>(ID));

  EXPECT_FORWARD_RET(OnNumericConstant, MakeNumericConstant,
                     (2.2), TestNumericExpr(ID));
  EXPECT_FORWARD_RET(OnVariableRef, MakeVariable, (33), TestReference(ID));
  EXPECT_FORWARD_RET(OnCommonExprRef, MakeCommonExpr, (33), TestReference(ID));
  EXPECT_FORWARD_RET(OnUnary, MakeUnary, (expr::ABS, TestNumericExpr(ID)),
                     TestNumericExpr(ID2));
  EXPECT_FORWARD_RET(OnBinary, MakeBinary,
                     (expr::ADD, TestNumericExpr(ID), TestNumericExpr(ID2)),
                     TestNumericExpr(ID3));
  EXPECT_FORWARD_RET(OnIf, MakeIf,
                     (TestLogicalExpr(ID), TestNumericExpr(ID2),
                      TestNumericExpr(ID3)), TestNumericExpr(ID4));

  EXPECT_FORWARD_RET(BeginPLTerm, BeginPLTerm, (44), TestPLTermBuilder(ID));
  EXPECT_FORWARD_RET(EndPLTerm, EndPLTerm,
                     (TestPLTermBuilder(ID), TestReference(ID2)),
                     TestNumericExpr(ID3));

  EXPECT_FORWARD_RET(BeginVarArg, BeginVarArg, (expr::MAX, 77),
                     TestVarArgExprBuilder(ID));
  EXPECT_FORWARD_RET(EndVarArg, EndVarArg, (TestVarArgExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginSum, BeginSum, (88), TestNumericExprBuilder(ID));
  EXPECT_FORWARD_RET(EndSum, EndSum, (TestNumericExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginCount, BeginCount, (99), TestCountExprBuilder(ID));
  EXPECT_FORWARD_RET(EndCount, EndCount, (TestCountExprBuilder(ID)),
                     TestCountExpr(ID2));

  EXPECT_FORWARD_RET(BeginNumberOf, BeginNumberOf, (11, TestNumericExpr(ID)),
                     TestNumberOfExprBuilder(ID2));
  EXPECT_FORWARD_RET(EndNumberOf, EndNumberOf, (TestNumberOfExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginSymbolicNumberOf, BeginSymbolicNumberOf,
                     (11, TestExpr(ID)), TestSymbolicNumberOfExprBuilder(ID2));
  EXPECT_FORWARD_RET(EndSymbolicNumberOf, EndSymbolicNumberOf,
                     (TestSymbolicNumberOfExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(OnLogicalConstant, MakeLogicalConstant, (true),
                     TestLogicalExpr(ID));
  EXPECT_FORWARD_RET(OnNot, MakeNot, (TestLogicalExpr(ID)),
                     TestLogicalExpr(ID2));
  EXPECT_FORWARD_RET(OnBinaryLogical, MakeBinaryLogical,
                     (expr::OR, TestLogicalExpr(ID), TestLogicalExpr(ID2)),
                     TestLogicalExpr(ID3));
  EXPECT_FORWARD_RET(OnRelational, MakeRelational,
                     (expr::LT, TestNumericExpr(ID), TestNumericExpr(ID2)),
                     TestLogicalExpr(ID3));
  EXPECT_FORWARD_RET(OnLogicalCount, MakeLogicalCount,
                     (expr::ATLEAST, TestNumericExpr(ID), TestCountExpr(ID2)),
                     TestLogicalExpr(ID3));
  EXPECT_FORWARD_RET(OnImplication, MakeImplication,
                     (TestLogicalExpr(ID), TestLogicalExpr(ID2),
                      TestLogicalExpr(ID3)), TestLogicalExpr(ID4));

  EXPECT_FORWARD_RET(BeginIteratedLogical, BeginIteratedLogical,
                     (expr::EXISTS, 22), TestIteratedLogicalExprBuilder(ID));
  EXPECT_FORWARD_RET(EndIteratedLogical, EndIteratedLogical,
                     (TestIteratedLogicalExprBuilder(ID)),
                     TestLogicalExpr(ID2));

  EXPECT_FORWARD_RET(BeginPairwise, BeginPairwise, (expr::ALLDIFF, 33),
                     TestPairwiseExprBuilder(ID));
  EXPECT_FORWARD_RET(EndPairwise, EndPairwise,
                     (TestPairwiseExprBuilder(ID)), TestLogicalExpr(ID2));

  EXPECT_FORWARD_RET(OnStringLiteral, MakeStringLiteral, (str), TestExpr(ID));

  EXPECT_FORWARD_RET(OnSymbolicIf, MakeSymbolicIf,
                     (TestLogicalExpr(ID), TestExpr(ID2), TestExpr(ID3)),
                     TestExpr(ID4));
}

TEST(NLProblemBuilderTest, OnHeader) {
  StrictMock<MockProblemBuilder> builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  NLHeader h;
  EXPECT_CALL(builder, SetInfo(testing::Ref(h)));
  adapter.OnHeader(h);
}

// Test that ProblemBuilderToNLAdapter passes 0 as the number of objectives
// if obj_index is set to SKIP_ALL_OBJS.
TEST(NLProblemBuilderTest, SkipAllObjs) {
  MockProblemBuilder builder;
  typedef mp::ProblemBuilderToNLAdapter<MockProblemBuilder> Adapter;
  Adapter adapter(builder, Adapter::SKIP_ALL_OBJS);
  auto header = NLHeader();
  header.num_objs = 10;
  EXPECT_CALL(builder, SetInfo(Field(&mp::ProblemInfo::num_objs, 0)));
  adapter.OnHeader(header);
}

// Test that ProblemBuilderToNLAdapter passes the total number of objectives
// if obj_index is set to NEED_ALL_OBJS.
TEST(NLProblemBuilderTest, NeedAllObjs) {
  MockProblemBuilder builder;
  typedef mp::ProblemBuilderToNLAdapter<MockProblemBuilder> Adapter;
  Adapter adapter(builder, Adapter::NEED_ALL_OBJS);
  auto header = NLHeader();
  header.num_objs = 10;
  EXPECT_CALL(builder, SetInfo(Field(&mp::ProblemInfo::num_objs, 10)));
  adapter.OnHeader(header);
}

// Test that ProblemBuilderToNLAdapter passes min(1, num_objs) as the
// number of objectives by default.
TEST(NLProblemBuilderTest, SingleObjective) {
  MockProblemBuilder builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  auto header = NLHeader();
  header.num_objs = 10;
  EXPECT_CALL(builder, SetInfo(Field(&mp::ProblemInfo::num_objs, 1)));
  adapter.OnHeader(header);
  header.num_objs = 0;
  EXPECT_CALL(builder, SetInfo(Field(&mp::ProblemInfo::num_objs, 0)));
  adapter.OnHeader(header);
}

TEST(NLProblemBuilderTest, OnVarBounds) {
  StrictMock<MockProblemBuilder> builder;
  EXPECT_CALL(builder, AddVar(7.7, 8.8, mp::var::INTEGER));
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  adapter.OnVarBounds(66, 7.7, 8.8);
}

TEST(NLProblemBuilderTest, OnObj) {
  StrictMock<MockProblemBuilder> builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  auto header = mp::NLHeader();
  header.num_objs = 1;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  adapter.OnHeader(header);
  auto expr = TestNumericExpr(ID);
  adapter.OnObj(0, mp::obj::MAX, expr);
  auto obj_builder = TestLinearObjBuilder(ID);
  EXPECT_CALL(builder, AddObj(mp::obj::MAX, expr, 11)).
      WillOnce(Return(obj_builder));
  EXPECT_EQ(obj_builder, adapter.OnLinearObjExpr(0, 11));
}

TEST(NLProblemBuilderTest, OnCon) {
  StrictMock<MockProblemBuilder> builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  auto header = mp::NLHeader();
  header.num_algebraic_cons = 1;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  adapter.OnHeader(header);
  auto expr = TestNumericExpr(ID);
  adapter.OnAlgebraicCon(0, expr);
  adapter.OnConBounds(0, 11, 22);
  auto con_builder = TestLinearConBuilder(ID);
  EXPECT_CALL(builder, AddCon(11, 22, expr, 33)).WillOnce(Return(con_builder));
  EXPECT_EQ(con_builder, adapter.OnLinearConExpr(0, 33));
}

TEST(NLProblemBuilderTest, OnLogicalCon) {
  StrictMock<MockProblemBuilder> builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  auto expr = TestLogicalExpr(ID);
  EXPECT_CALL(builder, AddCon(expr));
  adapter.OnLogicalCon(0, expr);
}

TEST(NLProblemBuilderTest, OnCommonExpr) {
  StrictMock<MockProblemBuilder> builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  auto header = mp::NLHeader();
  header.num_common_exprs_in_cons = 1;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  adapter.OnHeader(header);
  auto expr_builder = TestLinearExprBuilder(ID);
  EXPECT_CALL(builder, BeginCommonExpr(11)).WillOnce(Return(expr_builder));
  adapter.BeginCommonExpr(0, 11);
  auto expr = TestNumericExpr(ID);
  EXPECT_CALL(builder, EndCommonExpr(expr_builder, expr, 22));
  adapter.EndCommonExpr(expr_builder, expr, 22);
}

TEST(NLProblemBuilderTest, OnFunction) {
  StrictMock<MockProblemBuilder> builder;
  mp::ProblemBuilderToNLAdapter<MockProblemBuilder> adapter(builder);
  auto header = mp::NLHeader();
  header.num_funcs = 1;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  adapter.OnHeader(header);
  auto func = TestFunction(ID);
  fmt::StringRef name("f");
  EXPECT_CALL(builder, AddFunction(name, 11, mp::func::SYMBOLIC)).
      WillOnce(Return(func));
  adapter.OnFunction(0, name, 11, mp::func::SYMBOLIC);
  auto call_builder = TestCallExprBuilder(ID);
  EXPECT_CALL(builder, BeginCall(func, 11)).WillOnce(Return(call_builder));
  adapter.BeginCall(0, 11);
  EXPECT_CALL(builder, EndCall(call_builder)).
      WillOnce(Return(TestNumericExpr(ID)));
  adapter.EndCall(call_builder);
}

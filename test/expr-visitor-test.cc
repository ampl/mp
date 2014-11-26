/*
 Expression visitor tests

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
#include "mp/expr-visitor.h"

// IsSame<T, U>::VALUE is true iff T and U are the same type.
template <typename T, typename U>
struct IsSame {
  static const bool VALUE = false;
};

template <typename T>
struct IsSame<T, T> {
  static const bool VALUE = true;
};

// Checks if ExprTypes::T is a typedef of mp::T.
#define CHECK_EXPR_TYPE(T) \
  EXPECT_TRUE((IsSame<mp::internal::ExprTypes::T, mp::T>::VALUE))

TEST(ExprTypesTest, Typedefs) {
  using mp::internal::ExprTypes;
  CHECK_EXPR_TYPE(NumericExpr);
  CHECK_EXPR_TYPE(LogicalExpr);
  CHECK_EXPR_TYPE(NumericConstant);
  CHECK_EXPR_TYPE(Variable);
  CHECK_EXPR_TYPE(UnaryExpr);
  CHECK_EXPR_TYPE(BinaryExpr);
  CHECK_EXPR_TYPE(IfExpr);
  CHECK_EXPR_TYPE(PLTerm);
  CHECK_EXPR_TYPE(CallExpr);
  EXPECT_TRUE((IsSame<ExprTypes::VarArgExpr, mp::IteratedExpr>::VALUE));
  EXPECT_TRUE((IsSame<ExprTypes::SumExpr, mp::IteratedExpr>::VALUE));
  EXPECT_TRUE((IsSame<ExprTypes::NumberOfExpr, mp::IteratedExpr>::VALUE));
  CHECK_EXPR_TYPE(CountExpr);
  CHECK_EXPR_TYPE(LogicalConstant);
  CHECK_EXPR_TYPE(NotExpr);
  CHECK_EXPR_TYPE(BinaryLogicalExpr);
  CHECK_EXPR_TYPE(RelationalExpr);
  CHECK_EXPR_TYPE(LogicalCountExpr);
  CHECK_EXPR_TYPE(ImplicationExpr);
  CHECK_EXPR_TYPE(IteratedLogicalExpr);
  CHECK_EXPR_TYPE(PairwiseExpr);
}

// Dummy struct used as an argument to TestResult and TestLResult constructors
// to make sure they are not called accidentally.
struct Dummy {};

struct TestResult {
  TestResult(Dummy) {}
};
struct TestLResult {
  TestLResult(Dummy) {}
};

// Use different classes for Result and LResult to make sure that it works.
struct MockVisitor : mp::ExprVisitor<MockVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitNumericConstant, TestResult (NumericConstant e));
  MOCK_METHOD1(VisitVariable, TestResult (Variable e));

  // Unary expressions.
  MOCK_METHOD1(VisitMinus, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAbs, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitFloor, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitCeil, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitSqrt, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitPow2, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitExp, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitLog, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitLog10, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitSin, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitSinh, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitCos, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitCosh, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitTanh, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitTan, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAsin, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAsinh, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAcos, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAcosh, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAtan, TestResult (UnaryExpr e));
  MOCK_METHOD1(VisitAtanh, TestResult (UnaryExpr e));

  // Binary expressions.
  MOCK_METHOD1(VisitAdd, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitSub, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitLess, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitMul, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitDiv, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitIntDiv, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitMod, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitPow, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitPowConstBase, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitPowConstExp, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitAtan2, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitPrecision, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitRound, TestResult (BinaryExpr e));
  MOCK_METHOD1(VisitTrunc, TestResult (BinaryExpr e));

  MOCK_METHOD1(VisitIf, TestResult (IfExpr e));
  MOCK_METHOD1(VisitPLTerm, TestResult (PLTerm e));
  MOCK_METHOD1(VisitCall, TestResult (CallExpr e));
  MOCK_METHOD1(VisitMin, TestResult (VarArgExpr e));
  MOCK_METHOD1(VisitMax, TestResult (VarArgExpr e));
  MOCK_METHOD1(VisitSum, TestResult (SumExpr e));
  MOCK_METHOD1(VisitNumberOf, TestResult (NumberOfExpr e));
  MOCK_METHOD1(VisitCount, TestResult (CountExpr e));

  MOCK_METHOD1(VisitLogicalConstant, TestLResult (LogicalConstant e));
  MOCK_METHOD1(VisitNot, TestLResult (NotExpr e));

  MOCK_METHOD1(VisitOr, TestLResult (BinaryLogicalExpr e));
  MOCK_METHOD1(VisitAnd, TestLResult (BinaryLogicalExpr e));
  MOCK_METHOD1(VisitIff, TestLResult (BinaryLogicalExpr e));

  // Relational expressions.
  MOCK_METHOD1(VisitLT, TestLResult (RelationalExpr e));
  MOCK_METHOD1(VisitLE, TestLResult (RelationalExpr e));
  MOCK_METHOD1(VisitEQ, TestLResult (RelationalExpr e));
  MOCK_METHOD1(VisitGE, TestLResult (RelationalExpr e));
  MOCK_METHOD1(VisitGT, TestLResult (RelationalExpr e));
  MOCK_METHOD1(VisitNE, TestLResult (RelationalExpr e));

  // Logical count expressions.
  MOCK_METHOD1(VisitAtLeast, TestLResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitAtMost, TestLResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitExactly, TestLResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitNotAtLeast, TestLResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitNotAtMost, TestLResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitNotExactly, TestLResult (LogicalCountExpr e));

  MOCK_METHOD1(VisitImplication, TestLResult (ImplicationExpr e));
  MOCK_METHOD1(VisitForAll, TestLResult (IteratedLogicalExpr e));
  MOCK_METHOD1(VisitExists, TestLResult (IteratedLogicalExpr e));
  MOCK_METHOD1(VisitAllDiff, TestLResult (PairwiseExpr e));
  MOCK_METHOD1(VisitNotAllDiff, TestLResult (PairwiseExpr e));
};

class ExprVisitorTest : public ::testing::Test {
 protected:
  mp::ExprFactory factory_;
  ::testing::StrictMock<MockVisitor> visitor_;
  mp::Variable var_;

  ExprVisitorTest() {
    ::testing::DefaultValue<TestResult>::Set(TestResult(Dummy()));
    var_ = factory_.MakeVariable(0);
  }
};

TEST_F(ExprVisitorTest, VisitNumericConstant) {
  auto e = factory_.MakeNumericConstant(42);
  EXPECT_CALL(visitor_, VisitNumericConstant(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
}

TEST_F(ExprVisitorTest, VisitVariable) {
  auto e = factory_.MakeVariable(42);
  EXPECT_CALL(visitor_, VisitVariable(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
}

#define TEST_UNARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeUnary(mp::expr::KIND, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
  }

TEST_UNARY(MINUS, Minus)
TEST_UNARY(ABS, Abs)
TEST_UNARY(FLOOR, Floor)
TEST_UNARY(CEIL, Ceil)
TEST_UNARY(SQRT, Sqrt)
TEST_UNARY(POW2, Pow2)
TEST_UNARY(EXP, Exp)
TEST_UNARY(LOG, Log)
TEST_UNARY(LOG10, Log10)
TEST_UNARY(SIN, Sin)
TEST_UNARY(SINH, Sinh)
TEST_UNARY(COS, Cos)
TEST_UNARY(COSH, Cosh)
TEST_UNARY(TAN, Tan)
TEST_UNARY(TANH, Tanh)
TEST_UNARY(ASIN, Asin)
TEST_UNARY(ASINH, Asinh)
TEST_UNARY(ACOS, Acos)
TEST_UNARY(ACOSH, Acosh)
TEST_UNARY(ATAN, Atan)
TEST_UNARY(ATANH, Atanh)

#define TEST_BINARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeBinary(mp::expr::KIND, var_, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
  }

TEST_BINARY(ADD, Add)
TEST_BINARY(SUB, Sub)
TEST_BINARY(LESS, Less)
TEST_BINARY(MUL, Mul)
TEST_BINARY(DIV, Div)
TEST_BINARY(INT_DIV, IntDiv)
TEST_BINARY(MOD, Mod)
TEST_BINARY(POW, Pow)
TEST_BINARY(POW_CONST_BASE, PowConstBase)
TEST_BINARY(POW_CONST_EXP, PowConstExp)
TEST_BINARY(ATAN2, Atan2)
TEST_BINARY(PRECISION, Precision)
TEST_BINARY(ROUND, Round)
TEST_BINARY(TRUNC, Trunc)

TEST_F(ExprVisitorTest, VisitIf) {
  auto e = factory_.MakeIf(factory_.MakeLogicalConstant(true), var_, var_);
  EXPECT_CALL(visitor_, VisitIf(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
}

TEST_F(ExprVisitorTest, VisitPLTerm) {
  auto builder = factory_.BeginPLTerm(1);
  builder.AddSlope(-1);
  builder.AddBreakpoint(0);
  builder.AddSlope(1);
  auto e = factory_.EndPLTerm(builder, var_);
  EXPECT_CALL(visitor_, VisitPLTerm(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
}

TEST_F(ExprVisitorTest, VisitCall) {
  auto builder = factory_.BeginCall(factory_.AddFunction("f"), 0);
  auto e = factory_.EndCall(builder);
  EXPECT_CALL(visitor_, VisitCall(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
}

#define TEST_ITERATED(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto builder = factory_.BeginIterated(mp::expr::KIND, 1); \
    builder.AddArg(var_); \
    auto e = factory_.EndIterated(builder); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
  }

TEST_ITERATED(MIN, Min)
TEST_ITERATED(MAX, Max)
TEST_ITERATED(SUM, Sum)
TEST_ITERATED(NUMBEROF, NumberOf)

TEST_F(ExprVisitorTest, VisitCount) {
  auto builder = factory_.BeginCount(1);
  builder.AddArg(factory_.MakeLogicalConstant(true));
  auto e = factory_.EndCount(builder);
  EXPECT_CALL(visitor_, VisitCount(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
}

// TODO: test logical expressions

/*
struct TestVisitor : ExprVisitor<TestVisitor, TestResult, TestLResult> {
  TestResult VisitUnhandledNumericExpr(NumericExpr e) {
    TestResult result = {e};
    return result;
  }

  TestLResult VisitUnhandledLogicalExpr(LogicalExpr e) {
    TestLResult result = {e};
    return result;
  }
};

TEST_F(ExprTest, ExprVisitorForwardsUnhandled) {
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  for (int i = 0; i < size; ++i) {
    TestVisitor v;
    const OpInfo &info = OP_INFO[i];
    if (info.kind == ex::UNKNOWN || info.kind == ex::STRING) continue;
    expr raw = RawExpr(i - 1);
    Expr e(::MakeExpr(&raw));
    Expr result;
    if (NumericExpr ne = Cast<NumericExpr>(e))
      result = v.Visit(ne).expr;
    else
      result = v.Visit(Cast<LogicalExpr>(e)).expr;
    EXPECT_EQ(e, result);
  }
}

struct NullVisitor : ExprVisitor<NullVisitor, int, int> {};

TEST_F(ExprTest, ExprVisitorUnhandledThrows) {
  EXPECT_THROW(NullVisitor().Visit(MakeConst(0)), UnsupportedExprError);
  EXPECT_THROW(NullVisitor().Visit(l0), UnsupportedExprError);
}

TEST_F(ExprTest, ExprVisitorInvalidExpr) {
  expr raw = RawExpr(opcode(ex::CONSTANT));
  NumericExpr ne(::MakeExpr<NumericExpr>(&raw));
  LogicalExpr le(::MakeExpr<LogicalExpr>(&raw));
  raw.op = reinterpret_cast<efunc*>(-1);
#ifndef NDEBUG
  EXPECT_DEBUG_DEATH(NullVisitor().Visit(ne), "Assertion");
  EXPECT_DEBUG_DEATH(NullVisitor().Visit(le), "Assertion");
#endif
}
*/

// TODO

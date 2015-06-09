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
#include "mock-allocator.h"
#include "test-assert.h"
#include "mp/expr-visitor.h"

using ::testing::StrictMock;

namespace expr = mp::expr;

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
  EXPECT_TRUE((IsSame<ExprTypes::Variable, mp::Reference>::VALUE));
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

// Dummy struct used as an argument to the TestResult constructor
// to make sure it is not called accidentally.
struct Dummy {};

struct TestResult {
  TestResult(Dummy) {}
};

struct MockVisitor : mp::ExprVisitor<MockVisitor, TestResult> {
  MOCK_METHOD1(VisitNumericConstant, TestResult (NumericConstant e));
  MOCK_METHOD1(VisitVariable, TestResult (Variable e));
  MOCK_METHOD1(VisitCommonExpr, TestResult (CommonExpr e));

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

  MOCK_METHOD1(VisitLogicalConstant, TestResult (LogicalConstant e));
  MOCK_METHOD1(VisitNot, TestResult (NotExpr e));

  MOCK_METHOD1(VisitOr, TestResult (BinaryLogicalExpr e));
  MOCK_METHOD1(VisitAnd, TestResult (BinaryLogicalExpr e));
  MOCK_METHOD1(VisitIff, TestResult (BinaryLogicalExpr e));

  // Relational expressions.
  MOCK_METHOD1(VisitLT, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitLE, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitEQ, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitGE, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitGT, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitNE, TestResult (RelationalExpr e));

  // Logical count expressions.
  MOCK_METHOD1(VisitAtLeast, TestResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitAtMost, TestResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitExactly, TestResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitNotAtLeast, TestResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitNotAtMost, TestResult (LogicalCountExpr e));
  MOCK_METHOD1(VisitNotExactly, TestResult (LogicalCountExpr e));

  MOCK_METHOD1(VisitImplication, TestResult (ImplicationExpr e));
  MOCK_METHOD1(VisitExists, TestResult (IteratedLogicalExpr e));
  MOCK_METHOD1(VisitForAll, TestResult (IteratedLogicalExpr e));
  MOCK_METHOD1(VisitAllDiff, TestResult (PairwiseExpr e));
  MOCK_METHOD1(VisitNotAllDiff, TestResult (PairwiseExpr e));

  MOCK_METHOD1(VisitStringLiteral, TestResult (StringLiteral e));
  MOCK_METHOD1(VisitSymbolicIf, TestResult (SymbolicIfExpr e));
};

class ExprVisitorTest : public ::testing::Test {
 protected:
  mp::ExprFactory factory_;
  StrictMock<MockVisitor> visitor_;
  mp::Reference var_;
  mp::LogicalConstant false_;

  ExprVisitorTest() {
    ::testing::DefaultValue<TestResult>::Set(TestResult(Dummy()));
    ::testing::DefaultValue<TestResult>::Set(TestResult(Dummy()));
    var_ = factory_.MakeVariable(0);
    false_ = factory_.MakeLogicalConstant(false);
  }

  // Makes a test count expression.
  mp::CountExpr MakeCount() {
    auto builder = factory_.BeginCount(1);
    builder.AddArg(false_);
    return factory_.EndCount(builder);
  }

  // Makes a test iterated expression.
  mp::IteratedExpr MakeIterated(expr::Kind kind) {
    auto builder = factory_.BeginIterated(kind, 1);
    builder.AddArg(var_);
    return factory_.EndIterated(builder);
  }
};

// An expression visitor that doesn't handle anything.
struct NullVisitor : mp::ExprVisitor<NullVisitor, TestResult> {};

// A visitor for testing VisitNumeric and VisitLogical.
struct MockUnhandledVisitor :
    mp::ExprVisitor<MockUnhandledVisitor, TestResult> {
  MOCK_METHOD1(VisitNumeric, TestResult (NumericExpr e));
  MOCK_METHOD1(VisitLogical, TestResult (LogicalExpr e));
};

#define TEST_UNHANDLED(expr) \
  EXPECT_THROW_MSG(StrictMock<NullVisitor>().Visit(expr), \
                   mp::UnsupportedError, \
                   fmt::format("unsupported: {}", str(expr.kind())))

// Tests that UnsupportedError is thrown if expr is unhandled.
#define TEST_UNHANDLED_NUMERIC(expr) { \
    TEST_UNHANDLED(expr); \
    StrictMock<MockUnhandledVisitor> visitor; \
    EXPECT_CALL(visitor, VisitNumeric(expr)); \
    visitor.Visit(expr); \
  }

#define TEST_UNHANDLED_LOGICAL(expr) { \
    TEST_UNHANDLED(expr); \
    StrictMock<MockUnhandledVisitor> visitor; \
    EXPECT_CALL(visitor, VisitLogical(expr)); \
    visitor.Visit(expr); \
  }

TEST_F(ExprVisitorTest, VisitNumericConstant) {
  auto e = factory_.MakeNumericConstant(42);
  EXPECT_CALL(visitor_, VisitNumericConstant(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_NUMERIC(base);
}

TEST_F(ExprVisitorTest, VisitVariable) {
  auto e = factory_.MakeVariable(42);
  EXPECT_CALL(visitor_, VisitVariable(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_NUMERIC(base);
}

TEST_F(ExprVisitorTest, VisitCommonExpr) {
  auto e = factory_.MakeCommonExpr(0);
  EXPECT_CALL(visitor_, VisitCommonExpr(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_NUMERIC(base);
}

// A visitor for testing VisitUnary.
struct MockUnaryVisitor : mp::ExprVisitor<MockUnaryVisitor, TestResult> {
  MOCK_METHOD1(VisitUnary, TestResult (UnaryExpr e));
};

#define TEST_UNARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeUnary(expr::KIND, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockUnaryVisitor> unary_visitor; \
    EXPECT_CALL(unary_visitor, VisitUnary(e)); \
    unary_visitor.Visit(base); \
    TEST_UNHANDLED_NUMERIC(base); \
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

// A visitor for testing VisitBinary.
struct MockBinaryVisitor : mp::ExprVisitor<MockBinaryVisitor, TestResult> {
  MOCK_METHOD1(VisitBinary, TestResult (BinaryExpr e));
};

#define TEST_BINARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeBinary(expr::KIND, var_, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockBinaryVisitor> binary_visitor; \
    EXPECT_CALL(binary_visitor, VisitBinary(e)); \
    binary_visitor.Visit(base); \
    TEST_UNHANDLED_NUMERIC(base); \
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

// A visitor for testing VisitBinaryFunc.
struct MockBinaryFuncVisitor :
    mp::ExprVisitor<MockBinaryFuncVisitor, TestResult> {
  MOCK_METHOD1(VisitBinaryFunc, TestResult (BinaryExpr e));
};

#define TEST_BINARY_FUNC(KIND, name) \
  TEST_BINARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name##Func) { \
    auto e = factory_.MakeBinary(expr::KIND, var_, var_); \
    mp::NumericExpr base = e; \
    StrictMock<MockBinaryFuncVisitor> binary_visitor; \
    EXPECT_CALL(binary_visitor, VisitBinaryFunc(e)); \
    binary_visitor.Visit(base); \
  }

TEST_BINARY_FUNC(ATAN2, Atan2)
TEST_BINARY_FUNC(PRECISION, Precision)
TEST_BINARY_FUNC(ROUND, Round)
TEST_BINARY_FUNC(TRUNC, Trunc)

TEST_F(ExprVisitorTest, VisitIf) {
  auto e = factory_.MakeIf(false_, var_, var_);
  EXPECT_CALL(visitor_, VisitIf(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_NUMERIC(base);
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
  TEST_UNHANDLED_NUMERIC(base);
}

TEST_F(ExprVisitorTest, VisitCall) {
  auto builder = factory_.BeginCall(factory_.AddFunction("f", 0), 0);
  auto e = factory_.EndCall(builder);
  EXPECT_CALL(visitor_, VisitCall(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_NUMERIC(base);
}

// A visitor for testing VisitVarArg.
struct MockVarArgVisitor : mp::ExprVisitor<MockVarArgVisitor, TestResult> {
  MOCK_METHOD1(VisitVarArg, TestResult (VarArgExpr e));
};

#define TEST_ITERATED(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = MakeIterated(expr::KIND); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
    TEST_UNHANDLED_NUMERIC(base); \
  }

#define TEST_VARARG(KIND, name) \
  TEST_ITERATED(KIND, name) \
  TEST_F(ExprVisitorTest, VisitVarArg##name) { \
    auto e = MakeIterated(expr::KIND); \
    StrictMock<MockVarArgVisitor> vararg_visitor; \
    EXPECT_CALL(vararg_visitor, VisitVarArg(e)); \
    mp::NumericExpr base = e; \
    vararg_visitor.Visit(base); \
  }

TEST_VARARG(MIN, Min)
TEST_VARARG(MAX, Max)
TEST_ITERATED(SUM, Sum)
TEST_ITERATED(NUMBEROF, NumberOf)

TEST_F(ExprVisitorTest, VisitCount) {
  auto e = MakeCount();
  EXPECT_CALL(visitor_, VisitCount(e));
  mp::NumericExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_NUMERIC(base);
}

TEST_F(ExprVisitorTest, VisitLogicalConstant) {
  auto e = factory_.MakeLogicalConstant(true);
  EXPECT_CALL(visitor_, VisitLogicalConstant(e));
  mp::LogicalExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_LOGICAL(base);
}

TEST_F(ExprVisitorTest, VisitNot) {
  auto e = factory_.MakeNot(false_);
  EXPECT_CALL(visitor_, VisitNot(e));
  mp::LogicalExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_LOGICAL(base);
}

// A visitor for testing VisitBinaryLogical.
struct MockBinaryLogicalVisitor :
    mp::ExprVisitor<MockBinaryLogicalVisitor, TestResult> {
  MOCK_METHOD1(VisitBinaryLogical, TestResult (BinaryLogicalExpr e));
};

#define TEST_BINARY_LOGICAL(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeBinaryLogical(expr::KIND, false_, false_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockBinaryLogicalVisitor> binary_visitor; \
    EXPECT_CALL(binary_visitor, VisitBinaryLogical(e)); \
    binary_visitor.Visit(base); \
    TEST_UNHANDLED_LOGICAL(base); \
  }

TEST_BINARY_LOGICAL(OR, Or)
TEST_BINARY_LOGICAL(AND, And)
TEST_BINARY_LOGICAL(IFF, Iff)

// A visitor for testing VisitRelational.
struct MockRelationalVisitor :
    mp::ExprVisitor<MockRelationalVisitor, TestResult> {
  MOCK_METHOD1(VisitRelational, TestResult (RelationalExpr e));
};

#define TEST_RELATIONAL(name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeRelational(expr::name, var_, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockRelationalVisitor> rel_visitor; \
    EXPECT_CALL(rel_visitor, VisitRelational(e)); \
    rel_visitor.Visit(base); \
    TEST_UNHANDLED_LOGICAL(base); \
  }

TEST_RELATIONAL(LT)
TEST_RELATIONAL(LE)
TEST_RELATIONAL(EQ)
TEST_RELATIONAL(GE)
TEST_RELATIONAL(GT)
TEST_RELATIONAL(NE)

// A visitor for testing VisitLogicalCount.
struct MockLogicalCountVisitor :
    mp::ExprVisitor<MockLogicalCountVisitor, TestResult> {
  MOCK_METHOD1(VisitLogicalCount, TestResult (LogicalCountExpr e));
};

#define TEST_LOGICAL_COUNT(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeLogicalCount(expr::KIND, var_, MakeCount()); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockLogicalCountVisitor> count_visitor; \
    EXPECT_CALL(count_visitor, VisitLogicalCount(e)); \
    count_visitor.Visit(base); \
    TEST_UNHANDLED_LOGICAL(base); \
  }

TEST_LOGICAL_COUNT(ATLEAST, AtLeast)
TEST_LOGICAL_COUNT(ATMOST, AtMost)
TEST_LOGICAL_COUNT(EXACTLY, Exactly)
TEST_LOGICAL_COUNT(NOT_ATLEAST, NotAtLeast)
TEST_LOGICAL_COUNT(NOT_ATMOST, NotAtMost)
TEST_LOGICAL_COUNT(NOT_EXACTLY, NotExactly)

TEST_F(ExprVisitorTest, VisitImplication) {
  auto e = factory_.MakeImplication(false_, false_, false_);
  EXPECT_CALL(visitor_, VisitImplication(e));
  mp::LogicalExpr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED_LOGICAL(base);
}

// A visitor for testing VisitIteratedLogical.
struct MockIteratedLogicalVisitor :
    mp::ExprVisitor<MockIteratedLogicalVisitor, TestResult> {
  MOCK_METHOD1(VisitIteratedLogical, TestResult (IteratedLogicalExpr e));
};

#define TEST_ITERATED_LOGICAL(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto builder = factory_.BeginIteratedLogical(expr::KIND, 1); \
    builder.AddArg(false_); \
    auto e = factory_.EndIteratedLogical(builder); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockIteratedLogicalVisitor> iter_visitor; \
    EXPECT_CALL(iter_visitor, VisitIteratedLogical(e)); \
    iter_visitor.Visit(base); \
    TEST_UNHANDLED_LOGICAL(base); \
  }

TEST_ITERATED_LOGICAL(EXISTS, Exists)
TEST_ITERATED_LOGICAL(FORALL, ForAll)

#define TEST_PAIRWISE(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto builder = factory_.BeginPairwise(expr::KIND, 1); \
    builder.AddArg(var_); \
    auto e = factory_.EndPairwise(builder); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    TEST_UNHANDLED_LOGICAL(base); \
  }

TEST_PAIRWISE(ALLDIFF, AllDiff)
TEST_PAIRWISE(NOT_ALLDIFF, NotAllDiff)

TEST_F(ExprVisitorTest, VisitStringLiteral) {
  auto e = factory_.MakeStringLiteral("foo");
  EXPECT_CALL(visitor_, VisitStringLiteral(e));
  mp::Expr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED(base);
}

TEST_F(ExprVisitorTest, VisitSymbolicIf) {
  auto e = factory_.MakeSymbolicIf(false_, var_, var_);
  EXPECT_CALL(visitor_, VisitSymbolicIf(e));
  mp::Expr base = e;
  visitor_.Visit(base);
  TEST_UNHANDLED(base);
}

TEST_F(ExprVisitorTest, InvalidExpr) {
  using ::testing::_;
  typedef AllocatorRef< ::testing::NiceMock<MockAllocator> > Allocator;
  ::testing::NiceMock<MockAllocator> alloc;
  char buffer[100];
  mp::BasicExprFactory<Allocator> f((Allocator(&alloc)));

  // Corrupt a numeric expression and check if visiting it results in
  // assertion failure.
  EXPECT_CALL(alloc, allocate(_)).WillOnce(::testing::Return(buffer));
  auto e1 = f.MakeNumericConstant(0);
  std::fill(buffer, buffer + sizeof(buffer), 0);
  EXPECT_ASSERT(visitor_.Visit(e1), "invalid expression");

  // Corrupt a logical expression and check if visiting it results in
  // assertion failure.
  EXPECT_CALL(alloc, allocate(_)).WillOnce(::testing::Return(buffer));
  auto e2 = f.MakeLogicalConstant(false);
  std::fill(buffer, buffer + sizeof(buffer), 0);
  EXPECT_ASSERT(visitor_.Visit(e2), "invalid expression");
}

struct MockConverter : mp::ExprConverter<MockConverter, TestResult> {
  MOCK_METHOD1(VisitLT, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitLE, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitEQ, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitGE, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitGT, TestResult (RelationalExpr e));
  MOCK_METHOD1(VisitNE, TestResult (RelationalExpr e));
};

MATCHER_P2(IsRelational, kind, expr, "") {
  return arg.kind() == kind &&
      arg.lhs() == expr.lhs() && arg.rhs() == expr.rhs();
}

TEST(ExprConverterTest, ConvertLogicalCountToRelational) {
  mp::ExprFactory f;
  auto lhs = f.MakeNumericConstant(42);
  auto rhs = f.EndCount(f.BeginCount(0));
  MockConverter converter;
  auto e = f.MakeLogicalCount(expr::ATLEAST, lhs, rhs);
  EXPECT_CALL(converter, VisitLE(IsRelational(expr::LE, e)));
  converter.Visit(e);
  e = f.MakeLogicalCount(expr::ATMOST, lhs, rhs);
  EXPECT_CALL(converter, VisitGE(IsRelational(expr::GE, e)));
  converter.Visit(e);
  e = f.MakeLogicalCount(expr::EXACTLY, lhs, rhs);
  EXPECT_CALL(converter, VisitEQ(IsRelational(expr::EQ, e)));
  converter.Visit(e);
  e = f.MakeLogicalCount(expr::NOT_ATLEAST, lhs, rhs);
  EXPECT_CALL(converter, VisitGT(IsRelational(expr::GT, e)));
  converter.Visit(e);
  e = f.MakeLogicalCount(expr::NOT_ATMOST, lhs, rhs);
  EXPECT_CALL(converter, VisitLT(IsRelational(expr::LT, e)));
  converter.Visit(e);
  e = f.MakeLogicalCount(expr::NOT_EXACTLY, lhs, rhs);
  EXPECT_CALL(converter, VisitNE(IsRelational(expr::NE, e)));
  converter.Visit(e);
}

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
  MOCK_METHOD1(VisitExists, TestLResult (IteratedLogicalExpr e));
  MOCK_METHOD1(VisitForAll, TestLResult (IteratedLogicalExpr e));
  MOCK_METHOD1(VisitAllDiff, TestLResult (PairwiseExpr e));
  MOCK_METHOD1(VisitNotAllDiff, TestLResult (PairwiseExpr e));
};

class ExprVisitorTest : public ::testing::Test {
 protected:
  mp::ExprFactory factory_;
  StrictMock<MockVisitor> visitor_;
  mp::Variable var_;
  mp::LogicalConstant false_;

  ExprVisitorTest() {
    ::testing::DefaultValue<TestResult>::Set(TestResult(Dummy()));
    ::testing::DefaultValue<TestLResult>::Set(TestLResult(Dummy()));
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
  mp::IteratedExpr MakeIterated(mp::expr::Kind kind) {
    auto builder = factory_.BeginIterated(kind, 1);
    builder.AddArg(var_);
    return factory_.EndIterated(builder);
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

// A visitor for testing VisitUnary.
struct MockUnaryVisitor :
    mp::ExprVisitor<MockUnaryVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitUnary, TestResult (UnaryExpr e));
};

#define TEST_UNARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeUnary(mp::expr::KIND, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockUnaryVisitor> unary_visitor; \
    EXPECT_CALL(unary_visitor, VisitUnary(e)); \
    unary_visitor.Visit(base); \
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
struct MockBinaryVisitor :
    mp::ExprVisitor<MockBinaryVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitBinary, TestResult (BinaryExpr e));
};

#define TEST_BINARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeBinary(mp::expr::KIND, var_, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockBinaryVisitor> binary_visitor; \
    EXPECT_CALL(binary_visitor, VisitBinary(e)); \
    binary_visitor.Visit(base); \
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
    mp::ExprVisitor<MockBinaryFuncVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitBinaryFunc, TestResult (BinaryExpr e));
};

#define TEST_BINARY_FUNC(KIND, name) \
  TEST_BINARY(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name##Func) { \
    auto e = factory_.MakeBinary(mp::expr::KIND, var_, var_); \
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

// A visitor for testing VisitVarArg.
struct MockVarArgVisitor :
    mp::ExprVisitor<MockVarArgVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitVarArg, TestResult (VarArgExpr e));
};

#define TEST_ITERATED(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = MakeIterated(mp::expr::KIND); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::NumericExpr base = e; \
    visitor_.Visit(base); \
  }

#define TEST_VARARG(KIND, name) \
  TEST_ITERATED(KIND, name) \
  TEST_F(ExprVisitorTest, VisitVarArg##name) { \
    auto e = MakeIterated(mp::expr::KIND); \
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
}

TEST_F(ExprVisitorTest, VisitLogicalConstant) {
  auto e = factory_.MakeLogicalConstant(true);
  EXPECT_CALL(visitor_, VisitLogicalConstant(e));
  mp::LogicalExpr base = e;
  visitor_.Visit(base);
}

TEST_F(ExprVisitorTest, VisitNot) {
  auto e = factory_.MakeNot(false_);
  EXPECT_CALL(visitor_, VisitNot(e));
  mp::LogicalExpr base = e;
  visitor_.Visit(base);
}

// A visitor for testing VisitBinaryLogical.
struct MockBinaryLogicalVisitor :
    mp::ExprVisitor<MockBinaryLogicalVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitBinaryLogical, TestLResult (BinaryLogicalExpr e));
};

#define TEST_BINARY_LOGICAL(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeBinaryLogical(mp::expr::KIND, false_, false_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockBinaryLogicalVisitor> binary_visitor; \
    EXPECT_CALL(binary_visitor, VisitBinaryLogical(e)); \
    binary_visitor.Visit(base); \
  }

TEST_BINARY_LOGICAL(OR, Or)
TEST_BINARY_LOGICAL(AND, And)
TEST_BINARY_LOGICAL(IFF, Iff)

// A visitor for testing VisitRelational.
struct MockRelationalVisitor :
    mp::ExprVisitor<MockRelationalVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitRelational, TestLResult (RelationalExpr e));
};

#define TEST_RELATIONAL(name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeRelational(mp::expr::name, var_, var_); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockRelationalVisitor> rel_visitor; \
    EXPECT_CALL(rel_visitor, VisitRelational(e)); \
    rel_visitor.Visit(base); \
  }

TEST_RELATIONAL(LT)
TEST_RELATIONAL(LE)
TEST_RELATIONAL(EQ)
TEST_RELATIONAL(GE)
TEST_RELATIONAL(GT)
TEST_RELATIONAL(NE)

// A visitor for testing VisitLogicalCount.
struct MockLogicalCountVisitor :
    mp::ExprVisitor<MockLogicalCountVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitLogicalCount, TestLResult (LogicalCountExpr e));
};

#define TEST_LOGICAL_COUNT(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto e = factory_.MakeLogicalCount(mp::expr::KIND, var_, MakeCount()); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockLogicalCountVisitor> count_visitor; \
    EXPECT_CALL(count_visitor, VisitLogicalCount(e)); \
    count_visitor.Visit(base); \
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
}

// A visitor for testing VisitIteratedLogical.
struct MockIteratedLogicalVisitor :
    mp::ExprVisitor<MockIteratedLogicalVisitor, TestResult, TestLResult> {
  MOCK_METHOD1(VisitIteratedLogical, TestLResult (IteratedLogicalExpr e));
};

#define TEST_ITERATED_LOGICAL(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto builder = factory_.BeginIteratedLogical(mp::expr::KIND, 1); \
    builder.AddArg(false_); \
    auto e = factory_.EndIteratedLogical(builder); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
    StrictMock<MockIteratedLogicalVisitor> iter_visitor; \
    EXPECT_CALL(iter_visitor, VisitIteratedLogical(e)); \
    iter_visitor.Visit(base); \
  }

TEST_ITERATED_LOGICAL(EXISTS, Exists)
TEST_ITERATED_LOGICAL(FORALL, ForAll)

#define TEST_PAIRWISE(KIND, name) \
  TEST_F(ExprVisitorTest, Visit##name) { \
    auto builder = factory_.BeginPairwise(mp::expr::KIND, 1); \
    builder.AddArg(var_); \
    auto e = factory_.EndPairwise(builder); \
    EXPECT_CALL(visitor_, Visit##name(e)); \
    mp::LogicalExpr base = e; \
    visitor_.Visit(base); \
  }

TEST_PAIRWISE(ALLDIFF, AllDiff)
TEST_PAIRWISE(NOT_ALLDIFF, NotAllDiff)

TEST_F(ExprVisitorTest, InvalidExpr) {
  using ::testing::_;
  typedef AllocatorRef< ::testing::NiceMock<MockAllocator<char> > > Allocator;
  ::testing::NiceMock< MockAllocator<char> > alloc;
  char buffer[100];
  mp::BasicExprFactory<Allocator> f((Allocator(&alloc)));

  // Corrupt a numeric expression and check if visiting it results in
  // assertion failure.
  EXPECT_CALL(alloc, allocate(_)).WillOnce(::testing::Return(buffer));
  auto e1 = f.MakeNumericConstant(0);
  std::fill(buffer, buffer + sizeof(buffer), 0);
  EXPECT_ASSERT(visitor_.Visit(e1), "invalid numeric expression");

  // Corrupt a logical expression and check if visiting it results in
  // assertion failure.
  EXPECT_CALL(alloc, allocate(_)).WillOnce(::testing::Return(buffer));
  auto e2 = f.MakeLogicalConstant(false);
  std::fill(buffer, buffer + sizeof(buffer), 0);
  EXPECT_ASSERT(visitor_.Visit(e2), "invalid logical expression");
}

// TODO: test unhandled expressions

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
*/

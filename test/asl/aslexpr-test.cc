/*
 Tests of the the C++ interface to AMPL expression trees.

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

#include <algorithm>

#include "../gtest-extra.h"
#include "../util.h"

#include "mp/nl-reader.h"
#include "asl/aslbuilder.h"
#include "asl/aslexpr-visitor.h"
#include "asl/aslproblem.h"

namespace asl = mp::asl;
using asl::Function;
using asl::Cast;
using asl::Expr;
using asl::NumericExpr;
using asl::LogicalExpr;
using asl::UnaryExpr;
using asl::BinaryExpr;
using asl::VarArgExpr;
using asl::SumExpr;
using asl::CountExpr;
using asl::IfExpr;
using asl::PiecewiseLinearExpr;
using asl::NumericConstant;
using asl::Reference;
using asl::NumberOfExpr;
using asl::CallExpr;
using asl::LogicalConstant;
using asl::RelationalExpr;
using asl::NotExpr;
using asl::LogicalCountExpr;
using asl::BinaryLogicalExpr;
using asl::ImplicationExpr;
using asl::IteratedLogicalExpr;
using asl::PairwiseExpr;
using asl::StringLiteral;
using asl::ExprVisitor;
using asl::LinearTerm;
using asl::LinearExpr;

using mp::Error;
using mp::MakeArrayRef;
namespace ex = mp::expr;
namespace func = mp::func;

namespace {

expr RawExpr(int opcode) {
  expr e = expr();
  e.op = reinterpret_cast<efunc*>(opcode);
  return e;
}

class TestExpr : public Expr {
 public:
  static void TestProxy();
  static void TestExprIterator();

  template <typename ExprType>
  static ExprType MakeExpr(expr *e) { return Expr::Create<ExprType>(e); }
};

struct TestGrad {
  TestGrad *next;
  double coef;
  int varno;
};

double TestFunc(arglist *) { return 0; }
}  // namespace

namespace mp {
#ifndef NDEBUG
namespace internal {
template <>
bool Is<TestExpr>(ex::Kind kind) {
  return kind >= expr::FIRST_UNARY && kind <= expr::LAST_BINARY;
}
}
#endif

namespace asl {

template <>
class LinearExpr< asl::LinearTerm<TestGrad> > {
 private:
  TestGrad grad_;

 public:
  LinearExpr(const TestGrad &g) : grad_(g) {}
  LinearTerm<TestGrad> get() { return LinearTerm<TestGrad>(&grad_); }
};
}  // namespace asl
}  // namespace mp

template <typename ExprType>
ExprType MakeExpr(expr *e) { return TestExpr::MakeExpr<ExprType>(e); }

Expr MakeExpr(expr *e) { return TestExpr::MakeExpr<Expr>(e); }

void TestExpr::TestProxy() {
  expr e = RawExpr(nl_opcode(ex::DIV));
  Proxy<NumericExpr> p(&e);
  EXPECT_EQ(ex::DIV, p->kind());
}

void TestExpr::TestExprIterator() {
  using mp::asl::internal::ExprIterator;
  {
    ExprIterator<NumericExpr> i;
    EXPECT_EQ(ExprIterator<NumericExpr>(), i);
  }
  expr exprs[] = {
      RawExpr(nl_opcode(ex::DIV)),
      RawExpr(nl_opcode(ex::ADD)),
      RawExpr(nl_opcode(ex::ATAN)),
  };
  expr *const ptrs[] = {exprs, exprs + 1, exprs + 2};
  ExprIterator<NumericExpr> i(ptrs);
  EXPECT_EQ(ExprIterator<NumericExpr>(ptrs), i);
  EXPECT_NE(ExprIterator<NumericExpr>(), i);
  EXPECT_EQ(ex::DIV, (*i).kind());
  EXPECT_EQ(ex::DIV, i->kind());

  ExprIterator<NumericExpr> i2(++i);
  EXPECT_EQ(i2, i);
  EXPECT_NE(ExprIterator<NumericExpr>(ptrs), i);
  EXPECT_EQ(ExprIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(ex::ADD, i->kind());

  ExprIterator<NumericExpr> i3(i++);
  EXPECT_NE(i3, i);
  EXPECT_NE(ExprIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(ExprIterator<NumericExpr>(ptrs + 2), i);
  EXPECT_EQ(ex::ADD, i3->kind());
  EXPECT_EQ(ex::ATAN, i->kind());

  int index = 0;
  for (ExprIterator<NumericExpr>
      i(ptrs), e(ptrs + 3); i != e; ++i, ++index) {
    int code = static_cast<int>(reinterpret_cast<size_t>(ptrs[index]->op));
    EXPECT_EQ(code, nl_opcode(i->kind()));
  }
  EXPECT_EQ(3, index);
  std::vector<NumericExpr> vec;
  std::copy(ExprIterator<NumericExpr>(ptrs),
      ExprIterator<NumericExpr>(ptrs + 3), std::back_inserter(vec));
  EXPECT_EQ(ex::ADD, vec[1].kind());
}

class ExprTest : public ::testing::Test {
 protected:
  asl::internal::ASLBuilder builder;
  NumericExpr n1, n2;
  LogicalExpr l0, l1;

  enum {NUM_VARS = 50};

  NumericConstant MakeConst(double value) {
    return builder.MakeNumericConstant(value);
  }

  Reference MakeVariable(int index) { return builder.MakeVariable(index); }

  UnaryExpr MakeUnary(ex::Kind kind, NumericExpr arg) {
    return builder.MakeUnary(kind, arg);
  }

  BinaryExpr MakeBinary(ex::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder.MakeBinary(kind, lhs, rhs);
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr then_expr, NumericExpr else_expr) {
    return builder.MakeIf(condition, then_expr, else_expr);
  }

  BinaryLogicalExpr MakeBinaryLogical(
      ex::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    return builder.MakeBinaryLogical(kind, lhs, rhs);
  }

  RelationalExpr MakeRelational(ex::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder.MakeRelational(kind, lhs, rhs);
  }

  LogicalCountExpr MakeLogicalCount(
      ex::Kind kind, NumericExpr lhs, CountExpr rhs) {
    return builder.MakeLogicalCount(kind, lhs, rhs);
  }

  ImplicationExpr MakeImplication(
      LogicalExpr condition, LogicalExpr then_expr, LogicalExpr else_expr) {
    return builder.MakeImplication(condition, then_expr, else_expr);
  }

public:
  ExprTest() {
    mp::ProblemInfo info = mp::ProblemInfo();
    info.num_vars = NUM_VARS;
    info.num_objs = 1;
    info.num_funcs = 2;
    int flags = asl::internal::ASL_STANDARD_OPCODES | ASL_allow_missing_funcs;
    builder.SetInfo(info);
    builder.set_flags(flags);
    n1 = builder.MakeNumericConstant(1);
    n2 = builder.MakeNumericConstant(2);
    l0 = builder.MakeLogicalConstant(false);
    l1 = builder.MakeLogicalConstant(true);
  }
};

TEST_F(ExprTest, NumericKinds) {
  const ex::Kind kinds[] = {
      ex::NUMBER,
      ex::VARIABLE,
      ex::FIRST_UNARY,
      ex::FIRST_BINARY,
      ex::IF,
      ex::FIRST_VARARG,
      ex::SUM,
      ex::COUNT,
      ex::PLTERM,
      ex::NUMBEROF
  };
  int i = 0, n = sizeof(kinds) / sizeof(*kinds);
  EXPECT_GT(n, 0);
  EXPECT_EQ(ex::FIRST_NUMERIC, ex::FIRST_EXPR);
  for (; i < n; ++i) {
    ex::Kind kind = kinds[i];
    EXPECT_GE(kind, ex::FIRST_NUMERIC);
    EXPECT_LE(kind, ex::LAST_NUMERIC);
    for (int j = i + 1; j < n; ++j)
      EXPECT_NE(kind, kinds[j]);  // Check if all different.
    EXPECT_LT(kind, ex::FIRST_LOGICAL);
  }
  EXPECT_EQ(i, n);
}

TEST_F(ExprTest, LogicalKinds) {
  const ex::Kind kinds[] = {
      ex::BOOL,
      ex::NOT,
      ex::FIRST_RELATIONAL,
      ex::FIRST_BINARY_LOGICAL,
      ex::IMPLICATION,
      ex::FIRST_ITERATED_LOGICAL,
      ex::ALLDIFF
  };
  int i = 0, n = sizeof(kinds) / sizeof(*kinds);
  EXPECT_GT(n, 0);
  EXPECT_LT(ex::FIRST_LOGICAL, ex::LAST_EXPR);
  for (; i < n; ++i) {
    ex::Kind kind = kinds[i];
    EXPECT_GE(kind, ex::FIRST_LOGICAL);
    EXPECT_LE(kind, ex::LAST_LOGICAL);
    for (int j = i + 1; j < n; ++j)
      EXPECT_NE(kind, kinds[j]);  // Check if all different.
    EXPECT_GT(kind, ex::LAST_NUMERIC);
  }
  EXPECT_EQ(i, n);
}

TEST_F(ExprTest, Proxy) {
  TestExpr::TestProxy();
}

TEST_F(ExprTest, ExprIterator) {
  TestExpr::TestExprIterator();
}

TEST_F(ExprTest, ExprCtor) {
  {
    Expr e;
    EXPECT_FALSE(e);
  }
  {
    expr raw = RawExpr(nl_opcode(ex::SUB));
    Expr e(::MakeExpr(&raw));
    EXPECT_EQ(ex::SUB, e.kind());
  }
  {
    expr raw = RawExpr(nl_opcode(ex::OR));
    Expr e(::MakeExpr(&raw));
    EXPECT_EQ(ex::OR, e.kind());
  }
}

TEST_F(ExprTest, SafeBool) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = RawExpr(42);
  Expr e2(::MakeExpr(&raw2));
  EXPECT_TRUE(e2);
}

struct OpInfo {
  int code;
  ex::Kind kind;
};

const OpInfo OP_INFO[] = {
  {-1,                 ex::UNKNOWN},
  {ex::ADD,            ex::FIRST_BINARY},
  {ex::SUB,            ex::FIRST_BINARY},
  {ex::MUL,            ex::FIRST_BINARY},
  {ex::DIV,            ex::FIRST_BINARY},
  {ex::MOD,            ex::FIRST_BINARY},
  {ex::POW,            ex::FIRST_BINARY},
  {ex::LESS,           ex::FIRST_BINARY},
  { 7,                 ex::UNKNOWN},
  { 8,                 ex::UNKNOWN},
  { 9,                 ex::UNKNOWN},
  {10,                 ex::UNKNOWN},
  {ex::MIN,            ex::FIRST_VARARG},
  {ex::MAX,            ex::FIRST_VARARG},
  {ex::FLOOR,          ex::FIRST_UNARY},
  {ex::CEIL,           ex::FIRST_UNARY},
  {ex::ABS,            ex::FIRST_UNARY},
  {ex::MINUS,          ex::FIRST_UNARY},
  {17,                 ex::UNKNOWN},
  {18,                 ex::UNKNOWN},
  {19,                 ex::UNKNOWN},
  {ex::OR,             ex::FIRST_BINARY_LOGICAL},
  {ex::AND,            ex::FIRST_BINARY_LOGICAL},
  {ex::LT,             ex::FIRST_RELATIONAL},
  {ex::LE,             ex::FIRST_RELATIONAL},
  {ex::EQ,             ex::FIRST_RELATIONAL},
  {25,                 ex::UNKNOWN},
  {26,                 ex::UNKNOWN},
  {27,                 ex::UNKNOWN},
  {ex::GE,             ex::FIRST_RELATIONAL},
  {ex::GT,             ex::FIRST_RELATIONAL},
  {ex::NE,             ex::FIRST_RELATIONAL},
  {31,                 ex::UNKNOWN},
  {32,                 ex::UNKNOWN},
  {33,                 ex::UNKNOWN},
  {ex::NOT,            ex::NOT},
  {ex::IF,             ex::IF},
  {36,                 ex::UNKNOWN},
  {ex::TANH,           ex::FIRST_UNARY},
  {ex::TAN,            ex::FIRST_UNARY},
  {ex::SQRT,           ex::FIRST_UNARY},
  {ex::SINH,           ex::FIRST_UNARY},
  {ex::SIN,            ex::FIRST_UNARY},
  {ex::LOG10,          ex::FIRST_UNARY},
  {ex::LOG,            ex::FIRST_UNARY},
  {ex::EXP,            ex::FIRST_UNARY},
  {ex::COSH,           ex::FIRST_UNARY},
  {ex::COS,            ex::FIRST_UNARY},
  {ex::ATANH,          ex::FIRST_UNARY},
  {ex::ATAN2,          ex::FIRST_BINARY},
  {ex::ATAN,           ex::FIRST_UNARY},
  {ex::ASINH,          ex::FIRST_UNARY},
  {ex::ASIN,           ex::FIRST_UNARY},
  {ex::ACOSH,          ex::FIRST_UNARY},
  {ex::ACOS,           ex::FIRST_UNARY},
  {ex::SUM,            ex::SUM},
  {ex::TRUNC_DIV,      ex::FIRST_BINARY},
  {ex::PRECISION,      ex::FIRST_BINARY},
  {ex::ROUND,          ex::FIRST_BINARY},
  {ex::TRUNC,          ex::FIRST_BINARY},
  {ex::COUNT,          ex::COUNT},
  {ex::NUMBEROF,       ex::NUMBEROF},
  {ex::NUMBEROF_SYM,   ex::UNKNOWN},
  {ex::ATLEAST,        ex::FIRST_LOGICAL_COUNT},
  {ex::ATMOST,         ex::FIRST_LOGICAL_COUNT},
  {ex::PLTERM,         ex::PLTERM},
  {ex::IFSYM,          ex::UNKNOWN},
  {ex::EXACTLY,        ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_ATLEAST,    ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_ATMOST,     ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_EXACTLY,    ex::FIRST_LOGICAL_COUNT},
  {ex::AND,            ex::FIRST_ITERATED_LOGICAL},
  {ex::OR,             ex::FIRST_ITERATED_LOGICAL},
  {ex::IMPLICATION,    ex::IMPLICATION},
  {ex::IFF,            ex::FIRST_BINARY_LOGICAL},
  {ex::ALLDIFF,        ex::FIRST_PAIRWISE},
  {ex::NOT_ALLDIFF,    ex::FIRST_PAIRWISE},
  {ex::POW_CONST_EXP,  ex::FIRST_BINARY},
  {ex::POW2,           ex::FIRST_UNARY},
  {ex::POW_CONST_BASE, ex::FIRST_BINARY},
  {ex::CALL,           ex::CALL},
  {ex::NUMBER,         ex::NUMBER},
  {ex::STRING,         ex::STRING},
  {ex::VARIABLE,       ex::VARIABLE},
  {ex::LAST_EXPR + 1,  ex::UNKNOWN},
  {777,                ex::UNKNOWN}
};

template <typename ExprType>
void TestAssertInCreate(int opcode) {
  expr e = RawExpr(opcode);
  EXPECT_DEBUG_DEATH(
        MakeExpr<ExprType>(&e), "Assertion") << opcode;  // NOLINT(*)
}

template <typename ExprType>
std::size_t CheckExpr(ex::Kind start, ex::Kind end = ex::UNKNOWN,
    ex::Kind bad_kind = ex::PLTERM) {
  if (end == ex::UNKNOWN)
    end = start;
  {
    // Check default ctor.
    ExprType e;
    EXPECT_FALSE(e);
  }
  TestAssertInCreate<ExprType>(nl_opcode(bad_kind));
  std::size_t expr_count = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  for (int i = 0; i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int opcode = i - 1;
    expr raw = RawExpr(opcode);
    bool is_this_kind = (info.kind >= start && info.kind <= end) ||
        (info.kind == ex::NUMBER && start == ex::FIRST_LOGICAL);
    if (info.kind != ex::UNKNOWN) {
      Expr e(::MakeExpr(&raw));
      EXPECT_EQ(is_this_kind, Cast<ExprType>(e) != 0);
      bool cast_result = Cast<ExprType>(e);
      EXPECT_EQ(is_this_kind, cast_result);
    }
    if (!is_this_kind) continue;
    ExprType e(MakeExpr<ExprType>(&raw));
    EXPECT_EQ(opcode, nl_opcode(e.kind()));
    ++expr_count;
  }
  EXPECT_GT(expr_count, 0u);
  return expr_count;
}

TEST_F(ExprTest, Expr) {
  EXPECT_EQ(67u, CheckExpr<Expr>(ex::FIRST_EXPR, ex::LAST_EXPR, ex::UNKNOWN));
  TestAssertInCreate<Expr>(7);
  TestAssertInCreate<Expr>(sizeof(OP_INFO) / sizeof(*OP_INFO) - 1);
  TestAssertInCreate<Expr>(777);
}

// Test if Expr::Create() uses internal::Is() to check whether expression
// is of specified type. This allows testing Is() functions instead of doing
// expensive death tests. Is() is specialized for TestExpr to accept unary
// and binary but not other expression kinds.
TEST_F(ExprTest, CreateUsesIs) {
  expr raw1 = RawExpr(nl_opcode(ex::ADD));  // binary
  ::MakeExpr<TestExpr>(&raw1);
  expr raw2 = RawExpr(nl_opcode(ex::MINUS));  // unary
  ::MakeExpr<TestExpr>(&raw2);
  TestAssertInCreate<TestExpr>(nl_opcode(ex::PLTERM));  // neither
}

TEST_F(ExprTest, EqualityOperator) {
  expr raw1 = expr(), raw2 = expr();
  EXPECT_TRUE(::MakeExpr(&raw1) != Expr());
  EXPECT_TRUE(::MakeExpr(&raw1) != ::MakeExpr(&raw2));
  EXPECT_TRUE(Expr() == Expr());
  EXPECT_TRUE(::MakeExpr(&raw1) == ::MakeExpr(&raw1));
  EXPECT_TRUE(::MakeExpr(&raw2) == ::MakeExpr(&raw2));
}

TEST_F(ExprTest, NumericExpr) {
  EXPECT_EQ(45u,
      CheckExpr<NumericExpr>(ex::FIRST_NUMERIC, ex::LAST_NUMERIC, ex::NOT));
}

TEST_F(ExprTest, LogicalExpr) {
  EXPECT_EQ(22u, CheckExpr<LogicalExpr>(ex::FIRST_LOGICAL, ex::LAST_LOGICAL));
}

TEST_F(ExprTest, NumericConstant) {
  EXPECT_EQ(1u, CheckExpr<NumericConstant>(ex::NUMBER));
  asl::NumericConstant expr = builder.MakeNumericConstant(42);
  EXPECT_EQ(ex::NUMBER, expr.kind());
  EXPECT_EQ(42, expr.value());
}

TEST_F(ExprTest, Variable) {
  EXPECT_EQ(1u, CheckExpr<Reference>(ex::VARIABLE));
  Reference var = builder.MakeVariable(0);
  EXPECT_EQ(ex::VARIABLE, var.kind());
  EXPECT_EQ(0, var.index());
  var = builder.MakeVariable(9);
  EXPECT_EQ(9, var.index());
  EXPECT_DEBUG_DEATH(builder.MakeVariable(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(
        builder.MakeVariable(NUM_VARS);, "Assertion");  // NOLINT(*)
}

TEST_F(ExprTest, UnaryExpr) {
  const ex::Kind kinds[] = {
      ex::FLOOR, ex::CEIL, ex::ABS, ex::MINUS, ex::TANH, ex::TAN, ex::SQRT,
      ex::SINH, ex::SIN, ex::LOG10, ex::LOG, ex::EXP, ex::COSH, ex::COS,
      ex::ATANH, ex::ATAN, ex::ASINH, ex::ASIN, ex::ACOSH, ex::ACOS, ex::POW2
  };
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds, CheckExpr<UnaryExpr>(ex::FIRST_UNARY));
  for (std::size_t i = 0; i < num_kinds; ++i) {
    asl::UnaryExpr expr = builder.MakeUnary(kinds[i], n1);
    EXPECT_EQ(kinds[i], expr.kind());
    EXPECT_EQ(n1, expr.arg());
  }
  EXPECT_THROW_MSG(builder.MakeUnary(ex::ADD, n1), Error,
    fmt::format("invalid unary expression kind {}", ex::ADD));
}

#define EXPECT_BINARY(expr, expected_opcode, expected_lhs, expected_rhs) \
  EXPECT_EQ(expected_opcode, expr.kind()); \
  EXPECT_EQ(expected_lhs, expr.lhs()); \
  EXPECT_EQ(expected_rhs, expr.rhs())

TEST_F(ExprTest, BinaryExpr) {
  const ex::Kind kinds[] = {
    ex::ADD, ex::SUB, ex::MUL, ex::DIV, ex::MOD, ex::POW, ex::LESS,
    ex::ATAN2, ex::TRUNC_DIV, ex::PRECISION, ex::ROUND, ex::TRUNC,
    ex::POW_CONST_BASE, ex::POW_CONST_EXP
  };
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds, CheckExpr<BinaryExpr>(ex::FIRST_BINARY));
  for (std::size_t i = 0; i < num_kinds; ++i) {
    asl::BinaryExpr expr = builder.MakeBinary(kinds[i], n1, n2);
    EXPECT_BINARY(expr, kinds[i], n1, n2);
  }
  EXPECT_THROW_MSG(builder.MakeBinary(ex::MINUS, n1, n2), Error,
    fmt::format("invalid binary expression kind {}", ex::MINUS));
}

TEST_F(ExprTest, IfExpr) {
  EXPECT_EQ(1u, CheckExpr<IfExpr>(ex::IF));
  asl::IfExpr expr = builder.MakeIf(l1, n1, n2);
  EXPECT_EQ(ex::IF, expr.kind());
  EXPECT_EQ(l1, expr.condition());
  EXPECT_EQ(n1, expr.then_expr());
  EXPECT_EQ(n2, expr.else_expr());
}

TEST_F(ExprTest, PiecewiseLinearExpr) {
  EXPECT_EQ(1u,
            CheckExpr<PiecewiseLinearExpr>(ex::PLTERM, ex::PLTERM, ex::ADD));
  enum { NUM_BREAKPOINTS = 2 };
  double breakpoints[NUM_BREAKPOINTS] = { 11, 22 };
  double slopes[NUM_BREAKPOINTS + 1] = {33, 44, 55};
  Reference var = builder.MakeVariable(2);
  asl::PiecewiseLinearExpr expr = builder.MakePiecewiseLinear(
      NUM_BREAKPOINTS, breakpoints, slopes, var);
  EXPECT_EQ(ex::PLTERM, expr.kind());
  EXPECT_EQ(NUM_BREAKPOINTS, expr.num_breakpoints());
  EXPECT_EQ(NUM_BREAKPOINTS + 1, expr.num_slopes());
  for (int i = 0; i < NUM_BREAKPOINTS; ++i) {
    EXPECT_EQ(breakpoints[i], expr.breakpoint(i));
    EXPECT_EQ(slopes[i], expr.slope(i));
  }
  EXPECT_EQ(slopes[NUM_BREAKPOINTS], expr.slope(NUM_BREAKPOINTS));
  EXPECT_EQ(var, expr.arg());
#ifndef NDEBUG
  EXPECT_DEBUG_DEATH(
      builder.MakePiecewiseLinear(-1, breakpoints, slopes, var);,
      "Assertion");  // NOLINT(*)
#endif
}

TEST_F(ExprTest, CallExpr) {
  EXPECT_EQ(1u, CheckExpr<CallExpr>(ex::CALL));
  enum {NUM_ARGS = 3};
  Function f = builder.RegisterFunction(
        "foo", TestFunc, NUM_ARGS, func::SYMBOLIC);
  const Expr args[NUM_ARGS] = {n1, n2, builder.MakeStringLiteral("abc")};
  CallExpr expr = builder.MakeCall(f, args);
  EXPECT_EQ(ex::CALL, expr.kind());
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  EXPECT_EQ(f, expr.function());
  int arg_index = 0;
  for (CallExpr::iterator
      i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
    EXPECT_EQ(args[arg_index], *i);
    EXPECT_EQ(args[arg_index], expr.arg(arg_index));
  }
  EXPECT_EQ(NUM_ARGS, arg_index);
}

TEST_F(ExprTest, VarArgExpr) {
  const ex::Kind kinds[] = {ex::MIN, ex::MAX};
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds, CheckExpr<VarArgExpr>(ex::FIRST_VARARG));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  for (size_t i = 0; i < num_kinds; ++i) {
    VarArgExpr expr = builder.MakeVarArg(kinds[i], args);
    EXPECT_EQ(kinds[i], expr.kind());
    int arg_index = 0;
    VarArgExpr::iterator j = expr.begin(), end = expr.end();
    for (; j != end; ++j, ++arg_index)
      EXPECT_EQ(args[arg_index], *j);
    EXPECT_FALSE(*j);
    EXPECT_EQ(NUM_ARGS, arg_index);
    j = expr.begin();
    VarArgExpr::iterator j2 = j++;
    EXPECT_EQ(args[0], *j2);
    EXPECT_EQ(args[1], *j);
  }
  EXPECT_THROW_MSG(builder.MakeVarArg(ex::MINUS, args), Error,
      fmt::format("invalid vararg expression kind {}", ex::MINUS));
}

TEST_F(ExprTest, SumExpr) {
  EXPECT_EQ(1u, CheckExpr<SumExpr>(ex::SUM));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  asl::SumExpr expr = builder.MakeSum(args);
  EXPECT_EQ(ex::SUM, expr.kind());
  int arg_index = 0;
  for (asl::SumExpr::iterator
      i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
    EXPECT_EQ(args[arg_index], *i);
  }
}

TEST_F(ExprTest, CountExpr) {
  EXPECT_EQ(1u, CheckExpr<CountExpr>(ex::COUNT));
  enum {NUM_ARGS = 2};
  LogicalExpr args[NUM_ARGS] = {l1, l0};
  asl::CountExpr expr = builder.MakeCount(args);
  EXPECT_EQ(ex::COUNT, expr.kind());
  int arg_index = 0;
  for (asl::CountExpr::iterator
      i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
    EXPECT_EQ(args[arg_index], *i);
  }
}

TEST_F(ExprTest, NumberOfExpr) {
  EXPECT_EQ(1u, CheckExpr<NumberOfExpr>(ex::NUMBEROF));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  asl::NumberOfExpr expr = builder.MakeNumberOf(args);
  EXPECT_EQ(ex::NUMBEROF, expr.kind());
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  for (int i = 0; i < NUM_ARGS; ++i)
    EXPECT_EQ(args[i], expr.arg(i));
#ifndef NDEBUG
  EXPECT_DEBUG_DEATH(
      builder.MakeNumberOf(MakeArrayRef(args, 0));, "Assertion");  // NOLINT(*)
#endif
}

TEST_F(ExprTest, LogicalConstant) {
  EXPECT_EQ(1u, CheckExpr<LogicalConstant>(ex::BOOL));
  asl::LogicalConstant expr = builder.MakeLogicalConstant(true);
  EXPECT_EQ(ex::NUMBER, expr.kind());
  EXPECT_TRUE(expr.value());
  EXPECT_FALSE(builder.MakeLogicalConstant(false).value());
}

TEST_F(ExprTest, NotExpr) {
  EXPECT_EQ(1u, CheckExpr<NotExpr>(ex::NOT));
  asl::NotExpr expr = builder.MakeNot(l1);
  EXPECT_EQ(ex::NOT, expr.kind());
  EXPECT_EQ(l1, expr.arg());
}

TEST_F(ExprTest, BinaryLogicalExpr) {
  const ex::Kind kinds[] = {ex::OR, ex::AND, ex::IFF};
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds, CheckExpr<BinaryLogicalExpr>(ex::FIRST_BINARY_LOGICAL));
  for (size_t i = 0; i < num_kinds; ++i) {
    asl::BinaryLogicalExpr expr =
        builder.MakeBinaryLogical(kinds[i], l1, l0);
    EXPECT_BINARY(expr, kinds[i], l1, l0);
  }
  EXPECT_THROW_MSG(builder.MakeBinaryLogical(ex::MINUS, l1, l0), Error,
    fmt::format("invalid binary logical expression kind {}", ex::MINUS));
}

TEST_F(ExprTest, RelationalExpr) {
  const ex::Kind kinds[] = {ex::LT, ex::LE, ex::EQ, ex::GE, ex::GT, ex::NE};
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds, CheckExpr<RelationalExpr>(ex::FIRST_RELATIONAL));
  for (size_t i = 0; i < num_kinds; ++i) {
    asl::RelationalExpr expr = builder.MakeRelational(kinds[i], n1, n2);
    EXPECT_BINARY(expr, kinds[i], n1, n2);
  }
  EXPECT_THROW_MSG(builder.MakeRelational(ex::MINUS, n1, n2), Error,
    fmt::format("invalid relational expression kind {}", ex::MINUS));
}

TEST_F(ExprTest, LogicalCountExpr) {
  const ex::Kind kinds[] = {
    ex::ATLEAST, ex::ATMOST, ex::EXACTLY,
    ex::NOT_ATLEAST, ex::NOT_ATMOST, ex::NOT_EXACTLY
  };
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds, CheckExpr<LogicalCountExpr>(ex::FIRST_LOGICAL_COUNT));
  asl::CountExpr count = builder.MakeCount(MakeArrayRef(&l1, 1));
  for (size_t i = 0; i < num_kinds; ++i) {
    asl::LogicalCountExpr expr =
        builder.MakeLogicalCount(kinds[i], n1, count);
    EXPECT_BINARY(expr, kinds[i], n1, count);
  }
  EXPECT_THROW_MSG(builder.MakeLogicalCount(ex::MINUS, n1, count), Error,
    fmt::format("invalid logical count expression kind {}", ex::MINUS));
}

TEST_F(ExprTest, ImplicationExpr) {
  EXPECT_EQ(1u, CheckExpr<ImplicationExpr>(ex::IMPLICATION));
  LogicalExpr condition = builder.MakeLogicalConstant(true);
  asl::ImplicationExpr expr = builder.MakeImplication(condition, l0, l1);
  EXPECT_EQ(ex::IMPLICATION, expr.kind());
  EXPECT_EQ(condition, expr.condition());
  EXPECT_EQ(l0, expr.then_expr());
  EXPECT_EQ(l1, expr.else_expr());
}

TEST_F(ExprTest, IteratedLogicalExpr) {
  const ex::Kind kinds[] = {ex::EXISTS, ex::FORALL};
  std::size_t num_kinds = sizeof(kinds) / sizeof(*kinds);
  EXPECT_EQ(num_kinds,
    CheckExpr<IteratedLogicalExpr>(ex::FIRST_ITERATED_LOGICAL));
  enum {NUM_ARGS = 3};
  LogicalExpr args[NUM_ARGS] = {l0, l1, builder.MakeLogicalConstant(false)};
  for (size_t i = 0; i < num_kinds; ++i) {
    IteratedLogicalExpr expr = builder.MakeIteratedLogical(kinds[i], args);
    EXPECT_EQ(kinds[i], expr.kind());
    int arg_index = 0;
    for (IteratedLogicalExpr::iterator
        i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
      EXPECT_EQ(args[arg_index], *i);
    }
  }
  EXPECT_THROW_MSG(builder.MakeIteratedLogical(ex::MINUS, args), Error,
        fmt::format("invalid iterated logical expression kind {}", ex::MINUS));
}

TEST_F(ExprTest, PairwiseExpr) {
  EXPECT_EQ(2u, CheckExpr<PairwiseExpr>(ex::FIRST_PAIRWISE));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  PairwiseExpr expr = builder.MakeAllDiff(args);
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  int index = 0;
  for (PairwiseExpr::iterator
       i = expr.begin(), end = expr.end(); i != end; ++i, ++index) {
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index], expr.arg(index));
  }
  EXPECT_EQ(NUM_ARGS, index);
}

TEST_F(ExprTest, StringLiteral) {
  EXPECT_EQ(1u, CheckExpr<StringLiteral>(ex::STRING));
  StringLiteral e = builder.MakeStringLiteral("abc");
  EXPECT_STREQ("abc", e.value());
}

struct TestLResult {
  asl::LogicalExpr expr;
};

TestLResult MakeResult(asl::LogicalExpr e) {
  TestLResult result = {e};
  return result;
}

TEST_F(ExprTest, LinearTerm) {
  TestGrad g = {0, 11, 22};
  LinearExpr< LinearTerm<TestGrad> > expr(g);
  LinearTerm<TestGrad> term = expr.get();
  EXPECT_EQ(11, term.coef());
  EXPECT_EQ(22, term.var_index());
}

struct TestTerm {
  typedef TestGrad Grad;
  Grad *grad_;
  TestTerm(Grad *g) : grad_(g) {}
};

TEST_F(ExprTest, LinearExprIterator) {
  auto g2 = TestGrad();
  auto g1 = TestGrad();
  g1.next = &g2;
  LinearExpr<TestTerm>::iterator i(&g1), i2(i), i3(&g2), end(0);
  EXPECT_TRUE(i == i2);
  EXPECT_TRUE(i != i3);
  EXPECT_TRUE(i != end);
  EXPECT_EQ(&g1, i->grad_);
  EXPECT_EQ(&g1, (*i).grad_);
  EXPECT_TRUE(++i == i3);
  EXPECT_TRUE(i == i3);
  EXPECT_TRUE(i2++ != i3);
  EXPECT_TRUE(i2 == i3);
  EXPECT_EQ(&g2, i->grad_);
}

TEST_F(ExprTest, IsZero) {
  EXPECT_TRUE(IsZero(MakeConst(0)));
  EXPECT_FALSE(IsZero(MakeConst(1)));
  EXPECT_FALSE(IsZero(MakeConst(-1)));
  EXPECT_FALSE(IsZero(MakeConst(10)));
  EXPECT_FALSE(IsZero(MakeVariable(0)));
}

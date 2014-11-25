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

#include "mp/nl.h"
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
using asl::Variable;
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
using mp::UnsupportedExprError;
using mp::InvalidNumericExprError;
using mp::InvalidLogicalExprError;
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
  static void TestArrayIterator();

  template <typename ExprT>
  static ExprT MakeExpr(expr *e) { return Expr::Create<ExprT>(e); }
};

struct TestGrad {
  TestGrad *next;
  double coef;
  int varno;
};

double TestFunc(arglist *) { return 0; }
}  // namespace

namespace mp {
namespace asl {
#ifndef NDEBUG
namespace internal {
template <>
bool Is<TestExpr>(ex::Kind kind) {
  return kind >= expr::FIRST_UNARY && kind <= expr::LAST_BINARY;
}
}
#endif

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

template <typename ExprT>
ExprT MakeExpr(expr *e) { return TestExpr::MakeExpr<ExprT>(e); }

Expr MakeExpr(expr *e) { return TestExpr::MakeExpr<Expr>(e); }

void TestExpr::TestProxy() {
  expr e = RawExpr(opcode(ex::DIV));
  Proxy<NumericExpr> p(&e);
  EXPECT_EQ(ex::DIV, p->kind());
}

void TestExpr::TestArrayIterator() {
  {
    ArrayIterator<NumericExpr> i;
    EXPECT_EQ(ArrayIterator<NumericExpr>(), i);
  }
  expr exprs[] = {
      RawExpr(opcode(ex::DIV)),
      RawExpr(opcode(ex::ADD)),
      RawExpr(opcode(ex::ATAN)),
  };
  expr *const ptrs[] = {exprs, exprs + 1, exprs + 2};
  ArrayIterator<NumericExpr> i(ptrs);
  EXPECT_EQ(ArrayIterator<NumericExpr>(ptrs), i);
  EXPECT_NE(ArrayIterator<NumericExpr>(), i);
  EXPECT_EQ(ex::DIV, (*i).kind());
  EXPECT_EQ(ex::DIV, i->kind());

  ArrayIterator<NumericExpr> i2(++i);
  EXPECT_EQ(i2, i);
  EXPECT_NE(ArrayIterator<NumericExpr>(ptrs), i);
  EXPECT_EQ(ArrayIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(ex::ADD, i->kind());

  ArrayIterator<NumericExpr> i3(i++);
  EXPECT_NE(i3, i);
  EXPECT_NE(ArrayIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(ArrayIterator<NumericExpr>(ptrs + 2), i);
  EXPECT_EQ(ex::ADD, i3->kind());
  EXPECT_EQ(ex::ATAN, i->kind());

  int index = 0;
  for (ArrayIterator<NumericExpr>
      i(ptrs), e(ptrs + 3); i != e; ++i, ++index) {
    int code = static_cast<int>(reinterpret_cast<size_t>(ptrs[index]->op));
    EXPECT_EQ(code, opcode(i->kind()));
  }
  EXPECT_EQ(3, index);
  std::vector<NumericExpr> vec;
  std::copy(ArrayIterator<NumericExpr>(ptrs),
      ArrayIterator<NumericExpr>(ptrs + 3), std::back_inserter(vec));
  EXPECT_EQ(ex::ADD, vec[1].kind());
}

class ExprTest : public ::testing::Test {
 protected:
  asl::internal::ASLBuilder builder;
  NumericExpr n1, n2;
  LogicalConstant l0, l1;

  enum {NUM_VARS = 50};

  NumericConstant MakeConst(double value) {
    return builder.MakeNumericConstant(value);
  }

  Variable MakeVariable(int index) { return builder.MakeVariable(index); }

  UnaryExpr MakeUnary(ex::Kind kind, NumericExpr arg) {
    return builder.MakeUnary(kind, arg);
  }

  BinaryExpr MakeBinary(ex::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder.MakeBinary(kind, lhs, rhs);
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return builder.MakeIf(condition, true_expr, false_expr);
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
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return builder.MakeImplication(condition, true_expr, false_expr);
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
      ex::FIRST_UNARY,
      ex::FIRST_BINARY,
      ex::FIRST_VARARG,
      ex::SUM,
      ex::COUNT,
      ex::IF,
      ex::PLTERM,
      ex::VARIABLE,
      ex::NUMBEROF,
      ex::CONSTANT
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
    if (kind != ex::CONSTANT)
      EXPECT_LT(kind, ex::FIRST_LOGICAL);
  }
  EXPECT_EQ(i, n);
}

TEST_F(ExprTest, LogicalKinds) {
  const ex::Kind kinds[] = {
      ex::CONSTANT,
      ex::FIRST_RELATIONAL,
      ex::NOT,
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
    if (kind != ex::CONSTANT)
      EXPECT_GT(kind, ex::LAST_NUMERIC);
  }
  EXPECT_EQ(i, n);
}

TEST_F(ExprTest, Proxy) {
  TestExpr::TestProxy();
}

TEST_F(ExprTest, ArrayIterator) {
  TestExpr::TestArrayIterator();
}

TEST_F(ExprTest, ExprCtor) {
  {
    Expr e;
    EXPECT_FALSE(e);
  }
  {
    expr raw = RawExpr(opcode(ex::SUB));
    Expr e(::MakeExpr(&raw));
    EXPECT_EQ(ex::SUB, e.kind());
  }
  {
    expr raw = RawExpr(opcode(ex::OR));
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
  const char *str;
  ex::Kind kind;
};

const OpInfo OP_INFO[] = {
  {-1,                 "unknown",               ex::UNKNOWN},
  {ex::ADD,            "+",                     ex::FIRST_BINARY},
  {ex::SUB,            "-",                     ex::FIRST_BINARY},
  {ex::MUL,            "*",                     ex::FIRST_BINARY},
  {ex::DIV,            "/",                     ex::FIRST_BINARY},
  {ex::MOD,            "mod",                   ex::FIRST_BINARY},
  {ex::POW,            "^",                     ex::FIRST_BINARY},
  {ex::LESS,           "less",                  ex::FIRST_BINARY},
  { 7,                 "unknown",               ex::UNKNOWN},
  { 8,                 "unknown",               ex::UNKNOWN},
  { 9,                 "unknown",               ex::UNKNOWN},
  {10,                 "unknown",               ex::UNKNOWN},
  {ex::MIN,            "min",                   ex::FIRST_VARARG},
  {ex::MAX,            "max",                   ex::FIRST_VARARG},
  {ex::FLOOR,          "floor",                 ex::FIRST_UNARY},
  {ex::CEIL,           "ceil",                  ex::FIRST_UNARY},
  {ex::ABS,            "abs",                   ex::FIRST_UNARY},
  {ex::MINUS,          "unary -",               ex::FIRST_UNARY},
  {17,                 "unknown",               ex::UNKNOWN},
  {18,                 "unknown",               ex::UNKNOWN},
  {19,                 "unknown",               ex::UNKNOWN},
  {ex::OR,             "||",                    ex::FIRST_BINARY_LOGICAL},
  {ex::AND,            "&&",                    ex::FIRST_BINARY_LOGICAL},
  {ex::LT,             "<",                     ex::FIRST_RELATIONAL},
  {ex::LE,             "<=",                    ex::FIRST_RELATIONAL},
  {ex::EQ,             "=",                     ex::FIRST_RELATIONAL},
  {25,                 "unknown",               ex::UNKNOWN},
  {26,                 "unknown",               ex::UNKNOWN},
  {27,                 "unknown",               ex::UNKNOWN},
  {ex::GE,             ">=",                    ex::FIRST_RELATIONAL},
  {ex::GT,             ">",                     ex::FIRST_RELATIONAL},
  {ex::NE,             "!=",                    ex::FIRST_RELATIONAL},
  {31,                 "unknown",               ex::UNKNOWN},
  {32,                 "unknown",               ex::UNKNOWN},
  {33,                 "unknown",               ex::UNKNOWN},
  {ex::NOT,            "!",                     ex::NOT},
  {ex::IF,             "if",                    ex::IF},
  {36,                 "unknown",               ex::UNKNOWN},
  {ex::TANH,           "tanh",                  ex::FIRST_UNARY},
  {ex::TAN,            "tan",                   ex::FIRST_UNARY},
  {ex::SQRT,           "sqrt",                  ex::FIRST_UNARY},
  {ex::SINH,           "sinh",                  ex::FIRST_UNARY},
  {ex::SIN,            "sin",                   ex::FIRST_UNARY},
  {ex::LOG10,          "log10",                 ex::FIRST_UNARY},
  {ex::LOG,            "log",                   ex::FIRST_UNARY},
  {ex::EXP,            "exp",                   ex::FIRST_UNARY},
  {ex::COSH,           "cosh",                  ex::FIRST_UNARY},
  {ex::COS,            "cos",                   ex::FIRST_UNARY},
  {ex::ATANH,          "atanh",                 ex::FIRST_UNARY},
  {ex::ATAN2,          "atan2",                 ex::FIRST_BINARY},
  {ex::ATAN,           "atan",                  ex::FIRST_UNARY},
  {ex::ASINH,          "asinh",                 ex::FIRST_UNARY},
  {ex::ASIN,           "asin",                  ex::FIRST_UNARY},
  {ex::ACOSH,          "acosh",                 ex::FIRST_UNARY},
  {ex::ACOS,           "acos",                  ex::FIRST_UNARY},
  {ex::SUM,            "sum",                   ex::SUM},
  {ex::INT_DIV,        "div",                   ex::FIRST_BINARY},
  {ex::PRECISION,      "precision",             ex::FIRST_BINARY},
  {ex::ROUND,          "round",                 ex::FIRST_BINARY},
  {ex::TRUNC,          "trunc",                 ex::FIRST_BINARY},
  {ex::COUNT,          "count",                 ex::COUNT},
  {ex::NUMBEROF,       "numberof",              ex::NUMBEROF},
  {ex::NUMBEROF_SYM,   "string numberof",       ex::UNKNOWN},
  {ex::ATLEAST,        "atleast",               ex::FIRST_LOGICAL_COUNT},
  {ex::ATMOST,         "atmost",                ex::FIRST_LOGICAL_COUNT},
  {ex::PLTERM,         "piecewise-linear term", ex::PLTERM},
  {ex::IFSYM,          "string if-then-else",   ex::UNKNOWN},
  {ex::EXACTLY,        "exactly",               ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_ATLEAST,    "!atleast",              ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_ATMOST,     "!atmost",               ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_EXACTLY,    "!exactly",              ex::FIRST_LOGICAL_COUNT},
  {ex::AND,            "forall",                ex::FIRST_ITERATED_LOGICAL},
  {ex::OR,             "exists",                ex::FIRST_ITERATED_LOGICAL},
  {ex::IMPLICATION,    "==>",                   ex::IMPLICATION},
  {ex::IFF,            "<==>",                  ex::FIRST_BINARY_LOGICAL},
  {ex::ALLDIFF,        "alldiff",               ex::FIRST_PAIRWISE},
  {ex::NOT_ALLDIFF,    "!alldiff",              ex::FIRST_PAIRWISE},
  {ex::POW_CONST_EXP,  "^",                     ex::FIRST_BINARY},
  {ex::POW2,           "^2",                    ex::FIRST_UNARY},
  {ex::POW_CONST_BASE, "^",                     ex::FIRST_BINARY},
  {ex::CALL,           "function call",         ex::CALL},
  {ex::CONSTANT,       "constant",              ex::CONSTANT},
  {ex::STRING,         "string",                ex::STRING},
  {ex::VARIABLE,       "variable",              ex::VARIABLE},
  {ex::LAST_EXPR + 1,  "unknown",               ex::UNKNOWN},
  {777,                "unknown",               ex::UNKNOWN}
};

template <typename ExprT>
void TestAssertInCreate(int opcode) {
  expr e = RawExpr(opcode);
  EXPECT_DEBUG_DEATH(MakeExpr<ExprT>(&e), "Assertion") << opcode;  // NOLINT(*)
}

template <typename ExprT>
std::size_t CheckExpr(ex::Kind start, ex::Kind end = ex::UNKNOWN,
    ex::Kind bad_kind = ex::PLTERM) {
  if (end == ex::UNKNOWN)
    end = start;
  {
    // Check default ctor.
    ExprT e;
    EXPECT_FALSE(e);
  }
  TestAssertInCreate<ExprT>(opcode(bad_kind));
  std::size_t expr_count = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  for (int i = 0; i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int opcode = i - 1;
    const char *opstr = info.str;
    expr raw = RawExpr(opcode);
    bool is_this_kind = info.kind >= start && info.kind <= end;
    if (info.kind != ex::UNKNOWN) {
      Expr e(::MakeExpr(&raw));
      EXPECT_EQ(is_this_kind, asl::internal::Is<ExprT>(e));
      bool cast_result = Cast<ExprT>(e);
      EXPECT_EQ(is_this_kind, cast_result);
    }
    if (!is_this_kind) continue;
    ExprT e(MakeExpr<ExprT>(&raw));
    EXPECT_EQ(opcode, ex::opcode(e.kind()));
    EXPECT_STREQ(opstr, e.opstr());
    ++expr_count;
  }
  EXPECT_GT(expr_count, 0u);
  return expr_count;
}

TEST_F(ExprTest, Expr) {
  EXPECT_EQ(67, CheckExpr<Expr>(ex::FIRST_EXPR, ex::LAST_EXPR, ex::UNKNOWN));
  TestAssertInCreate<Expr>(7);
  TestAssertInCreate<Expr>(sizeof(OP_INFO) / sizeof(*OP_INFO) - 1);
  TestAssertInCreate<Expr>(777);
}

// Test if Expr::Create() uses internal::Is() to check whether expression
// is of specified type. This allows testing Is() functions instead of doing
// expensive death tests. Is() is specialized for TestExpr to accept unary
// and binary but not other expression kinds.
TEST_F(ExprTest, CreateUsesIs) {
  expr raw1 = RawExpr(opcode(ex::ADD));  // binary
  ::MakeExpr<TestExpr>(&raw1);
  expr raw2 = RawExpr(opcode(ex::MINUS));  // unary
  ::MakeExpr<TestExpr>(&raw2);
  TestAssertInCreate<TestExpr>(opcode(ex::PLTERM));  // neither
}

TEST_F(ExprTest, EqualityOperator) {
  expr raw1 = expr(), raw2 = expr();
  EXPECT_TRUE(::MakeExpr(&raw1) != Expr());
  EXPECT_TRUE(::MakeExpr(&raw1) != ::MakeExpr(&raw2));
  EXPECT_TRUE(Expr() == Expr());
  EXPECT_TRUE(::MakeExpr(&raw1) == ::MakeExpr(&raw1));
  EXPECT_TRUE(::MakeExpr(&raw2) == ::MakeExpr(&raw2));
}

TEST_F(ExprTest, EqualNumericConstant) {
  EXPECT_TRUE(Equal(MakeConst(0.42), MakeConst(0.42)));
  EXPECT_FALSE(Equal(MakeConst(0.42), MakeConst(42)));
}

TEST_F(ExprTest, EqualVariable) {
  EXPECT_TRUE(Equal(MakeVariable(0), MakeVariable(0)));
  EXPECT_FALSE(Equal(MakeVariable(0), MakeVariable(1)));
  EXPECT_FALSE(Equal(MakeVariable(0), MakeConst(0)));
}

TEST_F(ExprTest, EqualUnaryExpr) {
  NumericExpr e = MakeUnary(ex::MINUS, MakeVariable(0));
  EXPECT_TRUE(Equal(e, MakeUnary(ex::MINUS, MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeVariable(0)));
  EXPECT_FALSE(Equal(e, MakeUnary(ex::FLOOR, MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeUnary(ex::MINUS, MakeVariable(1))));
}

TEST_F(ExprTest, EqualBinaryExpr) {
  NumericExpr e = MakeBinary(ex::ADD, MakeVariable(0), MakeConst(42));
  EXPECT_TRUE(Equal(e, MakeBinary(ex::ADD, MakeVariable(0), MakeConst(42))));
  EXPECT_FALSE(Equal(e, MakeBinary(ex::SUB, MakeVariable(0), MakeConst(42))));
  EXPECT_FALSE(Equal(e, MakeBinary(ex::ADD, MakeConst(42), MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeBinary(ex::ADD, MakeVariable(0), MakeConst(0))));
  EXPECT_FALSE(Equal(MakeConst(42), e));
}

TEST_F(ExprTest, EqualIfExpr) {
  NumericExpr e =
      MakeIf(builder.MakeLogicalConstant(0), MakeVariable(1), MakeConst(42));
  EXPECT_TRUE(Equal(e,
      MakeIf(builder.MakeLogicalConstant(0), MakeVariable(1), MakeConst(42))));
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  EXPECT_FALSE(Equal(e, builder.MakeSum(args)));
  EXPECT_FALSE(Equal(e,
    MakeIf(builder.MakeLogicalConstant(0), MakeVariable(1), MakeConst(0))));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualPiecewiseLinear) {
  double breaks[] = {5, 10};
  double slopes[] = {-1, 0, 1};
  Variable x = MakeVariable(0), y = MakeVariable(1);
  NumericExpr e = builder.MakePiecewiseLinear(2, breaks, slopes, x);
  EXPECT_TRUE(Equal(e, builder.MakePiecewiseLinear(2, breaks, slopes, x)));
  EXPECT_FALSE(Equal(e, builder.MakePiecewiseLinear(1, breaks, slopes, x)));
  EXPECT_FALSE(Equal(e, builder.MakePiecewiseLinear(2, breaks, slopes, y)));
  double breaks2[] = {5, 11};
  EXPECT_FALSE(Equal(e, builder.MakePiecewiseLinear(2, breaks2, slopes, x)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualVarArgExpr) {
  NumericExpr args1[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args.
  NumericExpr args2[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  NumericExpr e = builder.MakeVarArg(ex::MIN, args1);
  EXPECT_TRUE(Equal(e, builder.MakeVarArg(ex::MIN, args2)));
  NumericExpr args3[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_FALSE(Equal(e, builder.MakeVarArg(ex::MIN, args3)));
  EXPECT_FALSE(Equal(builder.MakeVarArg(ex::MIN, args3),
                     builder.MakeVarArg(ex::MIN, args1)));
  EXPECT_FALSE(Equal(e, builder.MakeVarArg(ex::MAX, args1)));
  NumericExpr args4[] = {MakeVariable(0), MakeVariable(1), MakeConst(0)};
  EXPECT_FALSE(Equal(e, builder.MakeVarArg(ex::MIN, args4)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualSumExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  NumericExpr e = builder.MakeSum(args);
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args.
  NumericExpr args2[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  EXPECT_TRUE(Equal(e, builder.MakeSum(args2)));
  NumericExpr args3[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_FALSE(Equal(e, builder.MakeSum(args3)));
  EXPECT_FALSE(Equal(builder.MakeSum(args3), e));
  LogicalExpr args4[] = {l0, l1, l1};
  EXPECT_FALSE(Equal(e, builder.MakeCount(args4)));
  NumericExpr args5[] = {MakeVariable(0), MakeVariable(1), MakeConst(0)};
  EXPECT_FALSE(Equal(e, builder.MakeSum(args5)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualCountExpr) {
  LogicalExpr args[] = {l0, l1, l1};
  NumericExpr e = builder.MakeCount(args);
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args.
  LogicalExpr args2[] = {
    builder.MakeLogicalConstant(false),
    builder.MakeLogicalConstant(true),
    builder.MakeLogicalConstant(true)
  };
  EXPECT_TRUE(Equal(e, builder.MakeCount(args2)));
  LogicalExpr args3[] = {l0, l1};
  EXPECT_FALSE(Equal(e, builder.MakeCount(args3)));
  EXPECT_FALSE(Equal(builder.MakeCount(args3), e));
  NumericExpr args4[] = {MakeConst(0), MakeConst(1), MakeConst(1)};
  EXPECT_FALSE(Equal(e, builder.MakeSum(args4)));
  LogicalExpr args5[] = {l0, l1, l0};
  EXPECT_FALSE(Equal(e, builder.MakeCount(args5)));
}

TEST_F(ExprTest, NumericExpr) {
  EXPECT_EQ(45,
      CheckExpr<NumericExpr>(ex::FIRST_NUMERIC, ex::LAST_NUMERIC, ex::NOT));
}

TEST_F(ExprTest, LogicalExpr) {
  EXPECT_EQ(22, CheckExpr<LogicalExpr>(ex::FIRST_LOGICAL, ex::LAST_LOGICAL));
}

TEST_F(ExprTest, NumericConstant) {
  EXPECT_EQ(1, CheckExpr<NumericConstant>(ex::CONSTANT));
  asl::NumericConstant expr = builder.MakeNumericConstant(42);
  EXPECT_EQ(ex::CONSTANT, expr.kind());
  EXPECT_EQ(42, expr.value());
}

TEST_F(ExprTest, Variable) {
  EXPECT_EQ(1, CheckExpr<Variable>(ex::VARIABLE));
  asl::Variable var = builder.MakeVariable(0);
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
    ex::ATAN2, ex::INT_DIV, ex::PRECISION, ex::ROUND, ex::TRUNC,
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
  EXPECT_EQ(1, CheckExpr<IfExpr>(ex::IF));
  asl::IfExpr expr = builder.MakeIf(l1, n1, n2);
  EXPECT_EQ(ex::IF, expr.kind());
  EXPECT_EQ(l1, expr.condition());
  EXPECT_EQ(n1, expr.true_expr());
  EXPECT_EQ(n2, expr.false_expr());
}

TEST_F(ExprTest, PiecewiseLinearExpr) {
  EXPECT_EQ(1, CheckExpr<PiecewiseLinearExpr>(ex::PLTERM, ex::PLTERM, ex::ADD));
  enum { NUM_BREAKPOINTS = 2 };
  double breakpoints[NUM_BREAKPOINTS] = { 11, 22 };
  double slopes[NUM_BREAKPOINTS + 1] = {33, 44, 55};
  asl::Variable var = builder.MakeVariable(2);
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
  EXPECT_EQ(2, expr.var_index());
#ifndef NDEBUG
  EXPECT_DEBUG_DEATH(
      builder.MakePiecewiseLinear(-1, breakpoints, slopes, var);,
      "Assertion");  // NOLINT(*)
#endif
}

TEST_F(ExprTest, CallExpr) {
  EXPECT_EQ(1, CheckExpr<CallExpr>(ex::CALL));
  enum {NUM_ARGS = 3};
  Function f = builder.AddFunction("foo", TestFunc, NUM_ARGS, func::SYMBOLIC);
  const Expr args[NUM_ARGS] = {n1, n2, builder.MakeStringLiteral("abc")};
  CallExpr expr = builder.MakeCall(f, args);
  EXPECT_EQ(ex::CALL, expr.kind());
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  EXPECT_EQ(f, expr.function());
  int arg_index = 0;
  for (CallExpr::iterator
      i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
    EXPECT_EQ(args[arg_index], *i);
    EXPECT_EQ(args[arg_index], expr[arg_index]);
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
  EXPECT_EQ(1, CheckExpr<SumExpr>(ex::SUM));
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
  EXPECT_EQ(1, CheckExpr<CountExpr>(ex::COUNT));
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
  EXPECT_EQ(1, CheckExpr<NumberOfExpr>(ex::NUMBEROF));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  asl::NumberOfExpr expr = builder.MakeNumberOf(args);
  EXPECT_EQ(ex::NUMBEROF, expr.kind());
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  for (int i = 0; i < NUM_ARGS; ++i)
    EXPECT_EQ(args[i], expr[i]);
#ifndef NDEBUG
  EXPECT_DEBUG_DEATH(
      builder.MakeNumberOf(MakeArrayRef(args, 0));, "Assertion");  // NOLINT(*)
#endif
}

TEST_F(ExprTest, LogicalConstant) {
  EXPECT_EQ(1, CheckExpr<LogicalConstant>(ex::CONSTANT));
  asl::LogicalConstant expr = builder.MakeLogicalConstant(true);
  EXPECT_EQ(ex::CONSTANT, expr.kind());
  EXPECT_TRUE(expr.value());
  EXPECT_FALSE(builder.MakeLogicalConstant(false).value());
}

TEST_F(ExprTest, NotExpr) {
  EXPECT_EQ(1, CheckExpr<NotExpr>(ex::NOT));
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
  EXPECT_EQ(1, CheckExpr<ImplicationExpr>(ex::IMPLICATION));
  LogicalExpr condition = builder.MakeLogicalConstant(true);
  asl::ImplicationExpr expr = builder.MakeImplication(condition, l0, l1);
  EXPECT_EQ(ex::IMPLICATION, expr.kind());
  EXPECT_EQ(condition, expr.condition());
  EXPECT_EQ(l0, expr.true_expr());
  EXPECT_EQ(l1, expr.false_expr());
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
  EXPECT_EQ(2, CheckExpr<PairwiseExpr>(ex::FIRST_PAIRWISE));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  PairwiseExpr expr = builder.MakeAllDiff(args);
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  int index = 0;
  for (PairwiseExpr::iterator
       i = expr.begin(), end = expr.end(); i != end; ++i, ++index) {
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index], expr[index]);
  }
  EXPECT_EQ(NUM_ARGS, index);
}

TEST_F(ExprTest, StringLiteral) {
  EXPECT_EQ(1, CheckExpr<StringLiteral>(ex::STRING));
  StringLiteral e = builder.MakeStringLiteral("abc");
  EXPECT_STREQ("abc", e.value());
}

struct TestResult {
  NumericExpr expr;
};

TestResult MakeResult(NumericExpr e) {
  TestResult result = {e};
  return result;
}

struct TestLResult {
  LogicalExpr expr;
};

TestLResult MakeResult(LogicalExpr e) {
  TestLResult result = {e};
  return result;
}

// Use different classes for Result and LResult to check that it is possible.
struct FullTestVisitor : ExprVisitor<FullTestVisitor, TestResult, TestLResult> {
  TestResult VisitAdd(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitSub(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitMul(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitDiv(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitMod(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitPow(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitLess(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitMin(VarArgExpr e) { return MakeResult(e); }
  TestResult VisitMax(VarArgExpr e) { return MakeResult(e); }
  TestResult VisitFloor(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitCeil(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAbs(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitMinus(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitIf(IfExpr e) { return MakeResult(e); }
  TestResult VisitTanh(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitTan(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitSqrt(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitSinh(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitSin(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitLog10(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitLog(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitExp(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitCosh(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitCos(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAtanh(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAtan2(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitAtan(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAsinh(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAsin(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAcosh(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAcos(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitSum(SumExpr e) { return MakeResult(e); }
  TestResult VisitIntDiv(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitPrecision(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitRound(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitTrunc(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitCount(CountExpr e) { return MakeResult(e); }
  TestResult VisitNumberOf(NumberOfExpr e) { return MakeResult(e); }
  TestResult VisitCall(CallExpr e) { return MakeResult(e); }
  TestResult VisitPiecewiseLinear(PiecewiseLinearExpr e) {
    return MakeResult(e);
  }
  TestResult VisitPowConstExp(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitPow2(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitPowConstBase(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitNumericConstant(NumericConstant e) { return MakeResult(e); }
  TestResult VisitVariable(Variable e) { return MakeResult(e); }
  TestLResult VisitOr(BinaryLogicalExpr e) { return MakeResult(e); }
  TestLResult VisitAnd(BinaryLogicalExpr e) { return MakeResult(e); }
  TestLResult VisitLT(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitLE(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitEQ(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGE(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGT(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitNE(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitNot(NotExpr e) { return MakeResult(e); }
  TestLResult VisitAtLeast(LogicalCountExpr e) { return MakeResult(e); }
  TestLResult VisitAtMost(LogicalCountExpr e) { return MakeResult(e); }
  TestLResult VisitExactly(LogicalCountExpr e) { return MakeResult(e); }
  TestLResult VisitNotAtLeast(LogicalCountExpr e) { return MakeResult(e); }
  TestLResult VisitNotAtMost(LogicalCountExpr e) { return MakeResult(e); }
  TestLResult VisitNotExactly(LogicalCountExpr e) { return MakeResult(e); }
  TestLResult VisitForAll(IteratedLogicalExpr e) { return MakeResult(e); }
  TestLResult VisitExists(IteratedLogicalExpr e) { return MakeResult(e); }
  TestLResult VisitImplication(ImplicationExpr e) { return MakeResult(e); }
  TestLResult VisitIff(BinaryLogicalExpr e) { return MakeResult(e); }
  TestLResult VisitAllDiff(PairwiseExpr e) { return MakeResult(e); }
  TestLResult VisitNotAllDiff(PairwiseExpr e) { return MakeResult(e); }
  TestLResult VisitLogicalConstant(LogicalConstant e) { return MakeResult(e); }
};

TEST_F(ExprTest, ExprVisitorHandlesAll) {
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  for (int i = 0; i < size; ++i) {
    FullTestVisitor v;
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

struct TestConverter : asl::ExprConverter<TestConverter, void, TestLResult> {
  TestLResult VisitLT(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitLE(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitEQ(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGE(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGT(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitNE(RelationalExpr e) { return MakeResult(e); }
};

void CheckConversion(ex::Kind from_kind, ex::Kind to_kind) {
  expr lhs = RawExpr(opcode(ex::CONSTANT)), rhs = lhs;
  expr raw = RawExpr(opcode(from_kind));
  raw.L.e = &lhs;
  raw.R.e = &rhs;
  TestConverter converter;
  RelationalExpr expr =
      Cast<RelationalExpr>(converter.Visit(::MakeExpr<LogicalExpr>(&raw)).expr);
  EXPECT_EQ(to_kind, expr.kind());
  EXPECT_EQ(::MakeExpr(&lhs), expr.lhs());
  EXPECT_EQ(::MakeExpr(&rhs), expr.rhs());
}

TEST_F(ExprTest, ConvertLogicalCountToRelational) {
  CheckConversion(ex::ATLEAST, ex::LE);
  CheckConversion(ex::ATMOST, ex::GE);
  CheckConversion(ex::EXACTLY, ex::EQ);
  CheckConversion(ex::NOT_ATLEAST, ex::GT);
  CheckConversion(ex::NOT_ATMOST, ex::LT);
  CheckConversion(ex::NOT_EXACTLY, ex::NE);
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

struct Var {
  int index;
};

struct CreateVar {
 private:
  int index_;

 public:
  CreateVar() : index_(0) {}

  Var operator()() {
    Var v = {++index_};
    return v;
  }
};

TEST_F(ExprTest, NumberOfMap) {
  asl::NumberOfMap<Var, CreateVar> map((CreateVar()));
  EXPECT_TRUE(map.begin() == map.end());
  NumericExpr args1[] = {MakeConst(11), MakeVariable(0)};
  NumberOfExpr e1 = builder.MakeNumberOf(args1);
  NumericExpr args2[] = {MakeConst(22), MakeVariable(1)};
  NumberOfExpr e2 = builder.MakeNumberOf(args2);
  map.Add(11, e1);
  map.Add(22, e2);
  NumericExpr args3[] = {MakeConst(33), MakeVariable(0)};
  map.Add(33, builder.MakeNumberOf(args3));
  asl::NumberOfMap<Var, CreateVar>::iterator i = map.begin();
  EXPECT_EQ(e1, i->expr);
  EXPECT_EQ(2u, i->values.size());
  EXPECT_EQ(1, i->values.find(11)->second.index);
  EXPECT_EQ(3, i->values.find(33)->second.index);
  ++i;
  EXPECT_EQ(e2, i->expr);
  EXPECT_EQ(1u, i->values.size());
  EXPECT_EQ(2, i->values.find(22)->second.index);
  ++i;
  EXPECT_TRUE(i == map.end());
}

TEST_F(ExprTest, IsZero) {
  EXPECT_TRUE(IsZero(MakeConst(0)));
  EXPECT_FALSE(IsZero(MakeConst(1)));
  EXPECT_FALSE(IsZero(MakeConst(-1)));
  EXPECT_FALSE(IsZero(MakeConst(10)));
  EXPECT_FALSE(IsZero(MakeVariable(0)));
}

// Checks if WriteExpr produces the expected output for expr.
static ::testing::AssertionResult CheckWrite(
    const char *, const char *expr_str,
    const std::string &expected_output, NumericExpr expr) {
  fmt::MemoryWriter w;
  WriteExpr(w, asl::LinearObjExpr(), expr);
  std::string actual_output = w.str();
  if (expected_output == actual_output)
    return ::testing::AssertionSuccess();
  return ::testing::AssertionFailure()
      << "Output of: Writer(w, LinearObjExpr(), " << expr_str << ")\n"
      << "   Actual: " << actual_output << "\n"
      << " Expected: " << expected_output << "\n";
}

#define CHECK_WRITE(expected_output, expr) \
  EXPECT_PRED_FORMAT2(CheckWrite, expected_output, expr)

TEST_F(ExprTest, WriteNumericConstant) {
  CHECK_WRITE("0", MakeConst(0));
  CHECK_WRITE("42", MakeConst(42));
  CHECK_WRITE("12.34", MakeConst(12.34));
}

TEST_F(ExprTest, WriteVariable) {
  CHECK_WRITE("x1", MakeVariable(0));
  CHECK_WRITE("x3", MakeVariable(2));
}

TEST_F(ExprTest, WriteUnaryExpr) {
  NumericExpr x1 = MakeVariable(0);
  CHECK_WRITE("-x1", MakeUnary(ex::MINUS, x1));
  CHECK_WRITE("x1 ^ 2", MakeUnary(ex::POW2, x1));
  int count = 0;
  for (int i = 0, size = sizeof(OP_INFO) / sizeof(*OP_INFO); i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    ex::Kind code = static_cast<ex::Kind>(info.code);
    if (info.kind != ex::FIRST_UNARY || code == ex::MINUS || code == ex::POW2)
      continue;
    CHECK_WRITE(fmt::format("{}(x1)", info.str), MakeUnary(code, x1));
    ++count;
  }
  EXPECT_EQ(19, count);
}

TEST_F(ExprTest, WriteBinaryExpr) {
  Variable x1 = MakeVariable(0);
  NumericConstant n42 = MakeConst(42);
  CHECK_WRITE("x1 + 42", MakeBinary(ex::ADD, x1, n42));
  CHECK_WRITE("x1 - 42", MakeBinary(ex::SUB, x1, n42));
  CHECK_WRITE("x1 * 42", MakeBinary(ex::MUL, x1, n42));
  CHECK_WRITE("x1 / 42", MakeBinary(ex::DIV, x1, n42));
  CHECK_WRITE("x1 mod 42", MakeBinary(ex::MOD, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW_CONST_BASE, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW_CONST_EXP, x1, n42));
  CHECK_WRITE("x1 less 42", MakeBinary(ex::LESS, x1, n42));
  CHECK_WRITE("x1 div 42", MakeBinary(ex::INT_DIV, x1, n42));
}

TEST_F(ExprTest, WriteBinaryFunc) {
  auto x1 = MakeVariable(0);
  auto n42 = MakeConst(42);
  CHECK_WRITE("atan2(x1, 42)", MakeBinary(ex::ATAN2, x1, n42));
  CHECK_WRITE("precision(x1, 42)", MakeBinary(ex::PRECISION, x1, n42));
  CHECK_WRITE("round(x1, 42)", MakeBinary(ex::ROUND, x1, n42));
  CHECK_WRITE("trunc(x1, 42)", MakeBinary(ex::TRUNC, x1, n42));
}

TEST_F(ExprTest, WriteIfExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if x1 = 0 then 1",
      MakeIf(MakeRelational(ex::EQ, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 0 else 1",
      MakeIf(MakeRelational(ex::EQ, MakeVariable(0), n0), n0, n1));
}

TEST_F(ExprTest, WritePiecewiseLinearExpr) {
  double breaks[] = {5, 10};
  double slopes[] = {-1, 0, 1};
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43",
              builder.MakePiecewiseLinear(2, breaks, slopes, MakeVariable(42)));
}

TEST_F(ExprTest, WriteCallExpr) {
  asl::Function f = builder.AddFunction("foo", TestFunc, -1);
  Expr args[] = {
      MakeConst(3),
      MakeBinary(ex::ADD, MakeVariable(0), MakeConst(5)),
      MakeConst(7),
      MakeVariable(1)
  };
  CHECK_WRITE("foo()", builder.MakeCall(f, MakeArrayRef(args, 0)));
  CHECK_WRITE("foo(3)", builder.MakeCall(f, MakeArrayRef(args, 1)));
  CHECK_WRITE("foo(3, x1 + 5, 7)", builder.MakeCall(f, MakeArrayRef(args, 3)));
  CHECK_WRITE("foo(3, x1 + 5, 7, x2)", builder.MakeCall(f, args));
}

TEST_F(ExprTest, WriteVarArgExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  CHECK_WRITE("min(x1, x2, 42)", builder.MakeVarArg(ex::MIN, args));
  CHECK_WRITE("max(x1, x2, 42)", builder.MakeVarArg(ex::MAX, args));
}

TEST_F(ExprTest, WriteSumExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  CHECK_WRITE("/* sum */ (x1 + x2 + 42)", builder.MakeSum(args));
  NumericExpr args2[] = {
    MakeBinary(ex::ADD, MakeVariable(0), MakeVariable(1)), MakeConst(42)
  };
  CHECK_WRITE("/* sum */ ((x1 + x2) + 42)", builder.MakeSum(args2));
}

TEST_F(ExprTest, WriteCountExpr) {
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  CHECK_WRITE("count(x1 = 0, 1, 0)", builder.MakeCount(args));
}

TEST_F(ExprTest, WriteNumberOfExpr) {
  NumericExpr args[] = {MakeConst(42), MakeConst(43), MakeConst(44)};
  CHECK_WRITE("numberof 42 in (43, 44)", builder.MakeNumberOf(args));
}

TEST_F(ExprTest, WriteNotExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if !(x1 = 0) then 1",
      MakeIf(builder.MakeNot(
               MakeRelational(ex::EQ, MakeVariable(0), n0)), n1, n0));
}

TEST_F(ExprTest, WriteBinaryLogicalExpr) {
  auto e1 = MakeRelational(ex::GT, MakeVariable(0), MakeConst(0));
  auto e2 = MakeRelational(ex::LT, MakeVariable(0), MakeConst(10));
  CHECK_WRITE("if x1 > 0 || x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(ex::OR, e1, e2), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 > 0 && x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(ex::AND, e1, e2), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 > 0 <==> x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(ex::IFF, e1, e2), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteRelationalExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if x1 < 0 then 1",
      MakeIf(MakeRelational(ex::LT, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 <= 0 then 1",
      MakeIf(MakeRelational(ex::LE, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 1",
      MakeIf(MakeRelational(ex::EQ, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 >= 0 then 1",
      MakeIf(MakeRelational(ex::GE, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 > 0 then 1",
      MakeIf(MakeRelational(ex::GT, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 != 0 then 1",
      MakeIf(MakeRelational(ex::NE, MakeVariable(0), n0), n1, n0));
}

TEST_F(ExprTest, WriteLogicalCountExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1), value = MakeConst(42);
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  asl::CountExpr count = builder.MakeCount(args);
  CHECK_WRITE("if atleast 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::ATLEAST, value, count), n1, n0));
  CHECK_WRITE("if atmost 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::ATMOST, value, count), n1, n0));
  CHECK_WRITE("if exactly 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::EXACTLY, value, count), n1, n0));
  CHECK_WRITE("if !atleast 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATLEAST, value, count), n1, n0));
  CHECK_WRITE("if !atmost 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATMOST, value, count), n1, n0));
  CHECK_WRITE("if !exactly 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_EXACTLY, value, count), n1, n0));
}

TEST_F(ExprTest, WriteImplicationExpr) {
  auto e1 = MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0));
  CHECK_WRITE("if x1 = 0 ==> 1 then 1",
      MakeIf(MakeImplication(e1, l1, l0), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 = 0 ==> 0 else 1 then 1",
      MakeIf(MakeImplication(e1, l0, l1), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteIteratedLogicalExpr) {
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  CHECK_WRITE("if /* forall */ (x1 = 0 && 1 && 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ex::FORALL, args),
             MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if /* exists */ (x1 = 0 || 1 || 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ex::EXISTS, args),
             MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WritePairwiseExpr) {
  NumericExpr args[] = {MakeConst(42), MakeConst(43), MakeConst(44)};
  CHECK_WRITE("if alldiff(42, 43, 44) then 1",
      MakeIf(builder.MakeAllDiff(args), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteStringLiteral) {
  Expr args[] = {builder.MakeStringLiteral("abc")};
  asl::Function f = builder.AddFunction("f", TestFunc, 1, func::SYMBOLIC);
  CHECK_WRITE("f('abc')", builder.MakeCall(f, args));
  args[0] = builder.MakeStringLiteral("ab'c");
  CHECK_WRITE("f('ab''c')", builder.MakeCall(f, args));
  args[0] = builder.MakeStringLiteral("ab\nc");
  CHECK_WRITE("f('ab\\\nc')", builder.MakeCall(f, args));
}

TEST_F(ExprTest, UnaryExprPrecedence) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("--x1", MakeUnary(ex::MINUS, MakeUnary(ex::MINUS, x1)));
  CHECK_WRITE("-(x1 ^ x1)", MakeUnary(ex::MINUS, MakeBinary(ex::POW, x1, x1)));
}

TEST_F(ExprTest, UnaryFuncPrecedence) {
  auto x1 = MakeVariable(0);
  int count = 0;
  for (int i = 0, size = sizeof(OP_INFO) / sizeof(*OP_INFO); i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    ex::Kind code = static_cast<ex::Kind>(info.code);
    if (info.kind != ex::FIRST_UNARY || code == ex::MINUS || code == ex::POW2)
      continue;
    CHECK_WRITE(fmt::format("{0}({0}(x1))", info.str),
        MakeUnary(code, MakeUnary(code, x1)));
    CHECK_WRITE(fmt::format("{0}(x1 + x1)", info.str),
        MakeUnary(code, MakeBinary(ex::ADD, x1, x1)));
    ++count;
  }
  EXPECT_EQ(19, count);
}

TEST_F(ExprTest, Pow2Precedence) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("(x1 ^ 2) ^ 2", MakeUnary(ex::POW2, MakeUnary(ex::POW2, x1)));
  CHECK_WRITE("(x1 * x1) ^ 2",
              MakeUnary(ex::POW2, MakeBinary(ex::MUL, x1, x1)));
}

TEST_F(ExprTest, AdditiveExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 + x2 + x3",
      MakeBinary(ex::ADD, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("x1 + x2 - x3",
      MakeBinary(ex::SUB, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("x1 + x2 less x3",
      MakeBinary(ex::LESS, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("x1 + (x2 + x3)",
      MakeBinary(ex::ADD, x1, MakeBinary(ex::ADD, x2, x3)));
  CHECK_WRITE("(x1 + x2) * x3",
      MakeBinary(ex::MUL, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("if 1 then x1 else x2 + x3",
      MakeIf(l1, x1, MakeBinary(ex::ADD, x2, x3)));
}

TEST_F(ExprTest, MultiplicativeExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 * x2 * x3",
      MakeBinary(ex::MUL, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * x2 / x3",
      MakeBinary(ex::DIV, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * x2 div x3",
      MakeBinary(ex::INT_DIV, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * x2 mod x3",
      MakeBinary(ex::MOD, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * (x2 * x3)",
      MakeBinary(ex::MUL, x1, MakeBinary(ex::MUL, x2, x3)));
  CHECK_WRITE("(x1 * x2) ^ x3",
      MakeBinary(ex::POW, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("(x1 + x2) * x3",
      MakeBinary(ex::MUL, MakeBinary(ex::ADD, x1, x2), x3));
}

TEST_F(ExprTest, ExponentiationExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 ^ x2 ^ x3",
      MakeBinary(ex::POW, x1, MakeBinary(ex::POW, x2, x3)));
  CHECK_WRITE("x1 ^ x2 ^ 3",
      MakeBinary(ex::POW, x1, MakeBinary(ex::POW_CONST_BASE, x2,
                                         MakeConst(3))));
  CHECK_WRITE("x1 ^ 3 ^ x2",
      MakeBinary(ex::POW, x1, MakeBinary(ex::POW_CONST_EXP, MakeConst(3), x2)));
  CHECK_WRITE("(x1 ^ 2) ^ 3",
      MakeBinary(ex::POW_CONST_BASE, MakeBinary(ex::POW_CONST_BASE, x1,
                                                MakeConst(2)), MakeConst(3)));
  CHECK_WRITE("-x1 ^ -x2",
      MakeBinary(ex::POW, MakeUnary(ex::MINUS, x1), MakeUnary(ex::MINUS, x2)));
  CHECK_WRITE("x1 ^ (x2 * x3)",
      MakeBinary(ex::POW, x1, MakeBinary(ex::MUL, x2, x3)));
}

TEST_F(ExprTest, BinaryFuncPrecedence) {
  auto x1 = MakeVariable(0);
  auto e = MakeBinary(ex::ADD, x1, x1);
  CHECK_WRITE("atan2(atan2(x1, x1), x1 + x1)",
      MakeBinary(ex::ATAN2, MakeBinary(ex::ATAN2, x1, x1), e));
  CHECK_WRITE("precision(precision(x1, x1), x1 + x1)",
      MakeBinary(ex::PRECISION, MakeBinary(ex::PRECISION, x1, x1), e));
  CHECK_WRITE("round(round(x1, x1), x1 + x1)",
      MakeBinary(ex::ROUND, MakeBinary(ex::ROUND, x1, x1), e));
  CHECK_WRITE("trunc(trunc(x1, x1), x1 + x1)",
      MakeBinary(ex::TRUNC, MakeBinary(ex::TRUNC, x1, x1), e));
}

TEST_F(ExprTest, IfExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), n2 = MakeConst(2);
  LogicalExpr e = MakeBinaryLogical(ex::OR, l0, l1);
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1",
      MakeIf(e, MakeIf(e, n1, n0), n0));
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1 else 2",
      MakeIf(e, MakeIf(e, n1, n2), n0));
  CHECK_WRITE("if 0 || 1 then (if 0 || 1 then 1) else 2",
      MakeIf(e, MakeIf(e, n1, n0), n2));
  CHECK_WRITE("if 0 || 1 then 0 else if 0 || 1 then 1 else 2",
      MakeIf(e, n0, MakeIf(e, n1, n2)));
  CHECK_WRITE("if !(0 || 1) then x1 + 1",
      MakeIf(builder.MakeNot(e), MakeBinary(ex::ADD, MakeVariable(0), n1), n0));
}

TEST_F(ExprTest, PiecewiseLinearExprPrecedence) {
  double breakpoints[] = {5, 10};
  double slopes[] = {-1, 0, 1};
  NumericExpr e = builder.MakePiecewiseLinear(
        2, breakpoints, slopes, MakeVariable(42));
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43 ^ 2",
      MakeBinary(ex::POW, e, MakeConst(2)));
}

TEST_F(ExprTest, CallExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1);
  auto f = builder.AddFunction("foo", TestFunc, -1);
  Expr args[] = {
      builder.MakeCall(f, MakeArrayRef<Expr>(0, 0)),
      MakeBinary(ex::ADD, x1, MakeConst(5)),
      MakeConst(7),
      MakeUnary(ex::FLOOR, x2)
  };
  CHECK_WRITE("foo(foo(), x1 + 5, 7, floor(x2))", builder.MakeCall(f, args));
}

TEST_F(ExprTest, VarArgExprPrecedence) {
  NumericExpr x1 = MakeVariable(0), x2 = MakeVariable(1);
  NumericExpr e = MakeBinary(ex::ADD, x1, x2);
  NumericExpr args[] = {e, e};
  CHECK_WRITE("min(x1 + x2, x1 + x2)", builder.MakeVarArg(ex::MIN, args));
  CHECK_WRITE("max(x1 + x2, x1 + x2)", builder.MakeVarArg(ex::MAX, args));
  NumericExpr args2[] = {x1}, args3[] = {x2};
  NumericExpr args4[] = {
    builder.MakeVarArg(ex::MIN, args2), builder.MakeVarArg(ex::MIN, args3)
  };
  CHECK_WRITE("min(min(x1), min(x2))", builder.MakeVarArg(ex::MIN, args4));
  NumericExpr args5[] = {
    builder.MakeVarArg(ex::MAX, args2), builder.MakeVarArg(ex::MAX, args3)
  };
  CHECK_WRITE("max(max(x1), max(x2))", builder.MakeVarArg(ex::MAX, args5));
}

TEST_F(ExprTest, SumExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  NumericExpr args1[] = {x2, x3};
  NumericExpr args2[] = {x1, builder.MakeSum(args1)};
  CHECK_WRITE("/* sum */ (x1 + /* sum */ (x2 + x3))", builder.MakeSum(args2));
  NumericExpr args3[] = {x1, MakeBinary(ex::MUL, x2, x3)};
  CHECK_WRITE("/* sum */ (x1 + x2 * x3)", builder.MakeSum(args3));
  NumericExpr args4[] = {MakeBinary(ex::ADD, x1, x2), x3};
  CHECK_WRITE("/* sum */ ((x1 + x2) + x3)", builder.MakeSum(args4));
}

TEST_F(ExprTest, CountExprPrecedence) {
  LogicalExpr e = MakeBinaryLogical(ex::OR, l0, l1);
  LogicalExpr args[] = {e, e};
  CHECK_WRITE("count(0 || 1, 0 || 1)", builder.MakeCount(args));
}

TEST_F(ExprTest, NumberOfExprPrecedence) {
  NumericExpr x1 = MakeVariable(0), x2 = MakeVariable(1);
  NumericExpr args[] = {MakeConst(42), x1, x2};
  NumericExpr e = builder.MakeNumberOf(args);
  NumericExpr args2[] = {e, e, e};
  CHECK_WRITE("numberof numberof 42 in (x1, x2) in ("
      "numberof 42 in (x1, x2), numberof 42 in (x1, x2))",
              builder.MakeNumberOf(args2));
  NumericExpr e2 = MakeBinary(ex::ADD, x1, x2);
  NumericExpr args3[] = {e2, e2, e2};
  CHECK_WRITE("numberof x1 + x2 in (x1 + x2, x1 + x2)",
              builder.MakeNumberOf(args3));
}

TEST_F(ExprTest, NotExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), x = MakeVariable(0);
  CHECK_WRITE("if !!(x1 = 0) then 1",
      MakeIf(builder.MakeNot(
               builder.MakeNot(MakeRelational(ex::EQ, x, n0))), n1, n0));
}

TEST_F(ExprTest, LogicalOrExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 || 1 || 1 then 1",
      MakeIf(MakeBinaryLogical(ex::OR, MakeBinaryLogical(ex::OR, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 || (1 || 1) then 1",
      MakeIf(MakeBinaryLogical(ex::OR, l0, MakeBinaryLogical(ex::OR, l1, l1)),
          n1, n0));
  CHECK_WRITE("if 0 || 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::OR, l0, MakeBinaryLogical(ex::AND, l1, l1)),
          n1, n0));
}

TEST_F(ExprTest, LogicalAndExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 && 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::AND, MakeBinaryLogical(ex::AND, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 && (1 && 1) then 1",
      MakeIf(MakeBinaryLogical(ex::AND, l0, MakeBinaryLogical(ex::AND, l1, l1)),
          n1, n0));
  CHECK_WRITE("if 0 <= 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::AND,
                               MakeRelational(ex::LE, n0, n1), l1), n1, n0));
}

TEST_F(ExprTest, IffExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 <==> 1 <==> 1 then 1",
      MakeIf(MakeBinaryLogical(ex::IFF, MakeBinaryLogical(ex::IFF, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 <==> (1 <==> 1) then 1",
      MakeIf(MakeBinaryLogical(ex::IFF, l0, MakeBinaryLogical(ex::IFF, l1, l1)),
          n1, n0));
  CHECK_WRITE("if (0 <==> 1) && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::AND, MakeBinaryLogical(ex::IFF, l0, l1), l1),
          n1, n0));
}

TEST_F(ExprTest, RelationalExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  auto e1 = MakeBinary(ex::ADD, MakeVariable(0), n1);
  auto e2 = MakeBinary(ex::ADD, MakeVariable(1), n1);
  CHECK_WRITE("if x1 + 1 < x2 + 1 then 1",
      MakeIf(MakeRelational(ex::LT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 <= x2 + 1 then 1",
      MakeIf(MakeRelational(ex::LE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 = x2 + 1 then 1",
      MakeIf(MakeRelational(ex::EQ, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 >= x2 + 1 then 1",
      MakeIf(MakeRelational(ex::GE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 > x2 + 1 then 1",
      MakeIf(MakeRelational(ex::GT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 != x2 + 1 then 1",
      MakeIf(MakeRelational(ex::NE, e1, e2), n1, n0));
}

TEST_F(ExprTest, LogicalCountExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), lhs = MakeConst(42);
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1
  };
  LogicalExpr count1 =
      builder.MakeLogicalCount(ex::ATLEAST, lhs, builder.MakeCount(args));
  LogicalExpr args2[] = {count1, l1};
  CountExpr count2 = builder.MakeCount(args2);
  CHECK_WRITE("if atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::ATLEAST, lhs, count2), n1, n0));
  CHECK_WRITE("if atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::ATMOST, lhs, count2), n1, n0));
  CHECK_WRITE("if exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::EXACTLY, lhs, count2), n1, n0));
  CHECK_WRITE("if !atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATLEAST, lhs, count2), n1, n0));
  CHECK_WRITE("if !atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATMOST, lhs, count2), n1, n0));
  CHECK_WRITE("if !exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_EXACTLY, lhs, count2), n1, n0));

  args2[0] = l0;
  CountExpr count = builder.MakeCount(args2);
  CHECK_WRITE("if atleast 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::ATLEAST, lhs, count), l0), n1, n0));
  CHECK_WRITE("if atmost 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::ATMOST, lhs, count), l0), n1, n0));
  CHECK_WRITE("if exactly 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::EXACTLY, lhs, count), l0), n1, n0));
  CHECK_WRITE("if !atleast 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::NOT_ATLEAST, lhs, count), l0),
             n1, n0));
  CHECK_WRITE("if !atmost 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::NOT_ATMOST, lhs, count), l0),
             n1, n0));
  CHECK_WRITE("if !exactly 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::NOT_EXACTLY, lhs, count), l0),
             n1, n0));
}

TEST_F(ExprTest, IteratedLogicalExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  LogicalExpr args[] = {MakeBinaryLogical(ex::AND, l0, l0), l0};
  CHECK_WRITE("if /* forall */ ((0 && 0) && 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ex::FORALL, args), n1, n0));
  args[0] = MakeBinaryLogical(ex::OR, l0, l0);
  CHECK_WRITE("if /* exists */ ((0 || 0) || 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ex::EXISTS, args), n1, n0));
  args[0] = l0;
  LogicalExpr args2[] = {builder.MakeIteratedLogical(ex::FORALL, args), l0};
  CHECK_WRITE("if /* forall */ (/* forall */ (0 && 0) && 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ex::FORALL, args2), n1, n0));
  args2[0] = builder.MakeIteratedLogical(ex::EXISTS, args);
  CHECK_WRITE("if /* exists */ (/* exists */ (0 || 0) || 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ex::EXISTS, args2), n1, n0));
}

TEST_F(ExprTest, ImplicationExprPrecedence) {
  NumericConstant n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 ==> 1 ==> 0 then 1",
      MakeIf(MakeImplication(MakeImplication(l0, l1, l0), l0, l0), n1, n0));
  CHECK_WRITE("if 0 ==> 1 ==> 0 else 1 then 1",
      MakeIf(MakeImplication(MakeImplication(l0, l1, l0), l0, l1), n1, n0));
  CHECK_WRITE("if 0 ==> 1 else 0 ==> 1 then 1",
      MakeIf(MakeImplication(l0, l1, MakeImplication(l0, l1, l0)), n1, n0));
  CHECK_WRITE("if 0 ==> (1 ==> 0) else 1 then 1",
      MakeIf(MakeImplication(l0, MakeImplication(l1, l0, l0), l1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 ==> 0 else 1) then 1",
      MakeIf(MakeImplication(l0, MakeImplication(l1, l0, l1), l0), n1, n0));
  CHECK_WRITE("if 0 ==> 1 || 0 else 1 then 1",
              MakeIf(MakeImplication(
                       l0, MakeBinaryLogical(ex::OR, l1, l0), l1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 <==> 0) else 1 then 1",
              MakeIf(MakeImplication(
                       l0, MakeBinaryLogical(ex::IFF, l1, l0), l1), n1, n0));
}

TEST_F(ExprTest, PairwiseExprPrecedence) {
  NumericConstant n0 = MakeConst(0), n1 = MakeConst(1);
  NumericExpr args[] = {
    MakeBinary(ex::ADD, n0, n1), MakeBinary(ex::ADD, n0, n1)
  };
  CHECK_WRITE("if alldiff(0 + 1, 0 + 1) then 1",
      MakeIf(builder.MakeAllDiff(args), n1, n0));
}

#ifdef MP_USE_UNORDERED_MAP

using asl::internal::HashCombine;

TEST_F(ExprTest, HashNumericConstant) {
  size_t hash = HashCombine<int>(0, ex::CONSTANT);
  hash = HashCombine(hash, 42.0);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(builder.MakeNumericConstant(42)));
}

TEST_F(ExprTest, HashVariable) {
  size_t hash = HashCombine<int>(0, ex::VARIABLE);
  hash = HashCombine(hash, 42);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(MakeVariable(42)));
}

template <typename Base>
void CheckHash(asl::BasicUnaryExpr<Base> e) {
  size_t hash = HashCombine<int>(0, e.kind());
  hash = HashCombine<Base>(hash, e.arg());
  EXPECT_EQ(hash, std::hash<Base>()(e));
}

TEST_F(ExprTest, HashUnaryExpr) {
  CheckHash(builder.MakeUnary(ex::MINUS, builder.MakeVariable(0)));
}

template <typename Expr, typename Base, typename Arg>
void CheckHashBinary(Expr e) {
  size_t hash = HashCombine<int>(0, e.kind());
  hash = HashCombine<Arg>(hash, e.lhs());
  hash = HashCombine<Arg>(hash, e.rhs());
  EXPECT_EQ(hash, std::hash<Base>()(e));
}

template <typename Base, typename Arg>
void CheckHash(asl::BasicBinaryExpr<Base, Arg> e) {
  CheckHashBinary<asl::BasicBinaryExpr<Base, Arg>, Base, Arg>(e);
}

TEST_F(ExprTest, HashBinaryExpr) {
  CheckHash(builder.MakeBinary(ex::ADD, builder.MakeVariable(9), n2));
}

template <typename Base>
void CheckHash(asl::BasicIfExpr<Base> e) {
  size_t hash = HashCombine<int>(0, e.kind());
  hash = HashCombine<LogicalExpr>(hash, e.condition());
  hash = HashCombine<Base>(hash, e.true_expr());
  hash = HashCombine<Base>(hash, e.false_expr());
  EXPECT_EQ(hash, std::hash<Base>()(e));
}

TEST_F(ExprTest, HashIfExpr) {
  CheckHash(builder.MakeIf(l0, builder.MakeVariable(2), n1));
}

TEST_F(ExprTest, HashPiecewiseLinearExpr) {
  size_t hash = HashCombine<int>(0, ex::PLTERM);
  enum {NUM_BREAKPOINTS = 2};
  double breakpoints[NUM_BREAKPOINTS] = {5, 10};
  double slopes[NUM_BREAKPOINTS + 1] = {-1, 0, 1};
  for (size_t i = 0; i < NUM_BREAKPOINTS; ++i) {
    hash = HashCombine(hash, slopes[i]);
    hash = HashCombine(hash, breakpoints[i]);
  }
  hash = HashCombine(hash, slopes[NUM_BREAKPOINTS]);
  hash = HashCombine(hash, 9);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(
    builder.MakePiecewiseLinear(NUM_BREAKPOINTS, breakpoints, slopes,
                                builder.MakeVariable(9))));
}

TEST_F(ExprTest, HashCallExpr) {
  Variable var = MakeVariable(9);
  Expr args[2] = {n1, var};
  Function f = builder.AddFunction("foo", TestFunc, 2, func::SYMBOLIC);
  size_t hash = HashCombine<int>(0, ex::CALL);
  hash = HashCombine(hash, f.name());
  hash = HashCombine<NumericExpr>(hash, n1);
  hash = HashCombine<NumericExpr>(hash, var);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(builder.MakeCall(f, args)));
}

template <typename Expr, typename Arg, typename Base>
size_t CheckHash(Expr e) {
  size_t hash = HashCombine<int>(0, e.kind());
  for (typename Expr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    hash = HashCombine<Arg>(hash, *i);
  EXPECT_EQ(hash, std::hash<Base>()(e));
  return hash;
}

template <typename BaseT, typename ArgT, int ID>
void CheckHash(asl::BasicIteratedExpr<BaseT, ArgT, ID> e) {
  CheckHash<asl::BasicIteratedExpr<BaseT, ArgT, ID>, ArgT, BaseT>(e);
}

TEST_F(ExprTest, HashNumericVarArgExpr) {
  // Test computing hash for numeric expressions with abitrary numeric
  // arguments: VarArgExpr, SumExpr, NumberOfExpr.
  NumericExpr args[] = {
    builder.MakeVariable(0),
    builder.MakeVariable(1),
    builder.MakeNumericConstant(42)
  };
  CheckHash<VarArgExpr, NumericExpr, NumericExpr>(
        builder.MakeVarArg(ex::MIN, args));
  CheckHash(builder.MakeSum(args));
  CheckHash(builder.MakeNumberOf(args));
}

TEST_F(ExprTest, HashCountExpr) {
  LogicalExpr args[] = {l0, l1, l1};
  CheckHash(builder.MakeCount(args));
}

TEST_F(ExprTest, HashLogicalConstant) {
  size_t hash = HashCombine<int>(0, ex::CONSTANT);
  hash = HashCombine(hash, true);
  EXPECT_EQ(hash, std::hash<LogicalExpr>()(l1));
}

TEST_F(ExprTest, HashNotExpr) {
  CheckHash(builder.MakeNot(l1));
}

TEST_F(ExprTest, HashBinaryLogicalExpr) {
  CheckHash(builder.MakeBinaryLogical(ex::OR, l1, l0));
}

TEST_F(ExprTest, HashRelationalExpr) {
  CheckHash(builder.MakeRelational(ex::LT, builder.MakeVariable(6), n2));
}

TEST_F(ExprTest, HashLogicalCountExpr) {
  LogicalExpr args[] = {l0, l1};
  CheckHashBinary<LogicalCountExpr, LogicalExpr, NumericExpr>(
        builder.MakeLogicalCount(ex::ATMOST, n1, builder.MakeCount(args)));
}

TEST_F(ExprTest, HashImplicationExpr) {
  CheckHash(builder.MakeImplication(l1, l0, builder.MakeLogicalConstant(true)));
}

TEST_F(ExprTest, HashIteratedLogical) {
  LogicalExpr args[] = {l0, l1, builder.MakeLogicalConstant(false)};
  CheckHash(builder.MakeIteratedLogical(ex::EXISTS, args));
}

TEST_F(ExprTest, HashAllDiff) {
  NumericExpr args[] = {n1, n2, builder.MakeVariable(4)};
  CheckHash(builder.MakeAllDiff(args));
}

struct TestString {
  const char *str;
};

namespace std {

template <>
struct hash<TestString> {
  std::size_t operator()(TestString ts) const {
    size_t hash = asl::internal::HashCombine<int>(0, ex::STRING);
    for (const char *s = ts.str; *s; ++s)
      hash = asl::internal::HashCombine(hash, *s);
    return hash;
  }
};
}

TEST_F(ExprTest, HashStringLiteral) {
  // String literal can only occur as a function argument, so test
  // it as a part of a call expression.
  Expr args[] = {builder.MakeStringLiteral("test")};
  Function f = builder.AddFunction("foo", TestFunc, 1, func::SYMBOLIC);
  size_t hash = HashCombine<int>(0, ex::CALL);
  hash = HashCombine(hash, f.name());
  TestString ts = {"test"};
  hash = HashCombine(hash, ts);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(builder.MakeCall(f, args)));
}

TEST_F(ExprTest, HashNumberOfArgs) {
  size_t hash = HashCombine<NumericExpr>(0, MakeVariable(11));
  hash = HashCombine<NumericExpr>(hash, MakeConst(22));
  NumericExpr args[] = {MakeConst(42), MakeVariable(11), MakeConst(22)};
  EXPECT_EQ(hash, asl::internal::HashNumberOfArgs()(
              builder.MakeNumberOf(args)));
}

TEST_F(ExprTest, EqualNumberOfArgs) {
  using asl::internal::EqualNumberOfArgs;
  NumericExpr args1[] = {MakeConst(0), MakeVariable(11), MakeConst(22)};
  NumericExpr args2[] = {MakeConst(1), MakeVariable(11), MakeConst(22)};
  EXPECT_TRUE(EqualNumberOfArgs()(
                builder.MakeNumberOf(args1), builder.MakeNumberOf(args2)));
  NumericExpr args3[] = {MakeConst(1), MakeVariable(11)};
  EXPECT_FALSE(EqualNumberOfArgs()(
                 builder.MakeNumberOf(args1), builder.MakeNumberOf(args3)));
  NumericExpr args4[] = {MakeConst(1), MakeVariable(11), MakeConst(33)};
  EXPECT_FALSE(EqualNumberOfArgs()(
                 builder.MakeNumberOf(args1), builder.MakeNumberOf(args4)));
}

struct TestNumberOf {
  NumberOfExpr expr;

  TestNumberOf(NumberOfExpr e) : expr(e) {}
};

TEST_F(ExprTest, MatchNumberOfArgs) {
  using asl::internal::MatchNumberOfArgs;
  NumericExpr args1[] = {MakeConst(1), MakeVariable(11), MakeConst(22)};
  NumericExpr args2[] = {MakeConst(0), MakeVariable(11), MakeConst(22)};
  EXPECT_TRUE(MatchNumberOfArgs<TestNumberOf>(
      builder.MakeNumberOf(args1))(TestNumberOf(builder.MakeNumberOf(args2))));
  NumericExpr args3[] = {MakeConst(1), MakeVariable(11)};
  EXPECT_FALSE(MatchNumberOfArgs<TestNumberOf>(
      builder.MakeNumberOf(args3))(TestNumberOf(builder.MakeNumberOf(args2))));
  NumericExpr args4[] = {MakeConst(1), MakeVariable(11), MakeConst(33)};
  EXPECT_FALSE(MatchNumberOfArgs<TestNumberOf>(
      builder.MakeNumberOf(args4))(TestNumberOf(builder.MakeNumberOf(args2))));
}
#endif

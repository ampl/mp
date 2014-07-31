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

#include "gtest/gtest.h"
#include "tests/util.h"

#include "solvers/util/aslbuilder.h"
#include "solvers/util/nl.h"
#include "solvers/util/problem.h"

using ampl::Function;
using ampl::Cast;
using ampl::Expr;
using ampl::NumericExpr;
using ampl::LogicalExpr;
using ampl::UnaryExpr;
using ampl::BinaryExpr;
using ampl::VarArgExpr;
using ampl::SumExpr;
using ampl::CountExpr;
using ampl::IfExpr;
using ampl::PiecewiseLinearExpr;
using ampl::NumericConstant;
using ampl::Variable;
using ampl::NumberOfExpr;
using ampl::CallExpr;
using ampl::LogicalConstant;
using ampl::RelationalExpr;
using ampl::NotExpr;
using ampl::LogicalCountExpr;
using ampl::BinaryLogicalExpr;
using ampl::ImplicationExpr;
using ampl::IteratedLogicalExpr;
using ampl::AllDiffExpr;
using ampl::StringLiteral;
using ampl::ExprVisitor;
using ampl::internal::MakeArrayRef;

using ampl::LinearTerm;
using ampl::LinearExpr;

using ampl::Error;
using ampl::UnsupportedExprError;
using ampl::InvalidNumericExprError;
using ampl::InvalidLogicalExprError;

namespace {

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
}

namespace ampl {
#ifndef NDEBUG
namespace internal {
template <>
bool Is<TestExpr>(Expr e) {
  return e.kind() == Expr::BINARY || e.kind() == Expr::UNARY;
}
}
#endif

template <>
class LinearExpr< LinearTerm<TestGrad> > {
 private:
  TestGrad grad_;

 public:
  LinearExpr(const TestGrad &g) : grad_(g) {}
  LinearTerm<TestGrad> get() { return LinearTerm<TestGrad>(&grad_); }
};
}

namespace {

template <typename ExprT>
ExprT MakeExpr(expr *e) { return TestExpr::MakeExpr<ExprT>(e); }

Expr MakeExpr(expr *e) { return TestExpr::MakeExpr<Expr>(e); }

void TestExpr::TestProxy() {
  expr e = {reinterpret_cast<efunc*>(OPDIV)};
  Proxy<NumericExpr> p(&e);
  EXPECT_EQ(OPDIV, p->opcode());
}

void TestExpr::TestArrayIterator() {
  {
    ArrayIterator<NumericExpr> i;
    EXPECT_EQ(ArrayIterator<NumericExpr>(), i);
  }
  expr exprs[] = {
      {reinterpret_cast<efunc*>(OPDIV)},
      {reinterpret_cast<efunc*>(OPPLUS)},
      {reinterpret_cast<efunc*>(OP_atan)},
  };
  expr *const ptrs[] = {exprs, exprs + 1, exprs + 2};
  ArrayIterator<NumericExpr> i(ptrs);
  EXPECT_EQ(ArrayIterator<NumericExpr>(ptrs), i);
  EXPECT_NE(ArrayIterator<NumericExpr>(), i);
  EXPECT_EQ(OPDIV, (*i).opcode());
  EXPECT_EQ(OPDIV, i->opcode());

  ArrayIterator<NumericExpr> i2(++i);
  EXPECT_EQ(i2, i);
  EXPECT_NE(ArrayIterator<NumericExpr>(ptrs), i);
  EXPECT_EQ(ArrayIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(OPPLUS, i->opcode());

  ArrayIterator<NumericExpr> i3(i++);
  EXPECT_NE(i3, i);
  EXPECT_NE(ArrayIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(ArrayIterator<NumericExpr>(ptrs + 2), i);
  EXPECT_EQ(OPPLUS, i3->opcode());
  EXPECT_EQ(OP_atan, i->opcode());

  int index = 0;
  for (ArrayIterator<NumericExpr>
      i(ptrs), e(ptrs + 3); i != e; ++i, ++index) {
    int code = static_cast<int>(reinterpret_cast<size_t>(ptrs[index]->op));
    EXPECT_EQ(code, i->opcode());
  }
  EXPECT_EQ(3, index);
  std::vector<NumericExpr> vec;
  std::copy(ArrayIterator<NumericExpr>(ptrs),
      ArrayIterator<NumericExpr>(ptrs + 3), std::back_inserter(vec));
  EXPECT_EQ(OPPLUS, vec[1].opcode());
}

class ExprTest : public ::testing::Test {
 protected:
  ampl::internal::ASLBuilder builder;
  NumericExpr n1, n2;
  LogicalConstant l0, l1;

  enum {NUM_VARS = 50};

  NumericConstant MakeConst(double value) {
    return builder.MakeNumericConstant(value);
  }

  Variable MakeVariable(int index) { return builder.MakeVariable(index); }

  UnaryExpr MakeUnary(int opcode, NumericExpr arg) {
    return builder.MakeUnary(opcode, arg);
  }

  BinaryExpr MakeBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return builder.MakeBinary(opcode, lhs, rhs);
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return builder.MakeIf(condition, true_expr, false_expr);
  }

  BinaryLogicalExpr MakeBinaryLogical(
      int opcode, LogicalExpr lhs, LogicalExpr rhs) {
    return builder.MakeBinaryLogical(opcode, lhs, rhs);
  }

  RelationalExpr MakeRelational(int opcode, NumericExpr lhs, NumericExpr rhs) {
    return builder.MakeRelational(opcode, lhs, rhs);
  }

  LogicalCountExpr MakeLogicalCount(
      int opcode, NumericExpr lhs, CountExpr rhs) {
    return builder.MakeLogicalCount(opcode, lhs, rhs);
  }

  ImplicationExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return builder.MakeImplication(condition, true_expr, false_expr);
  }

public:
  ExprTest() {
    ampl::NLHeader header = {};
    header.num_vars = NUM_VARS;
    header.num_objs = 1;
    header.num_funcs = 2;
    int flags = ampl::internal::ASL_STANDARD_OPCODES | ASL_allow_missing_funcs;
    builder.BeginBuild("", header, flags);
    n1 = builder.MakeNumericConstant(1);
    n2 = builder.MakeNumericConstant(2);
    l0 = builder.MakeLogicalConstant(false);
    l1 = builder.MakeLogicalConstant(true);
  }
};

TEST_F(ExprTest, NumericKinds) {
  const Expr::Kind kinds[] = {
      Expr::UNARY,
      Expr::BINARY,
      Expr::VARARG,
      Expr::SUM,
      Expr::COUNT,
      Expr::IF,
      Expr::PLTERM,
      Expr::VARIABLE,
      Expr::NUMBEROF,
      Expr::CONSTANT
  };
  int i = 0, n = sizeof(kinds) / sizeof(*kinds);
  EXPECT_GT(n, 0);
  EXPECT_EQ(Expr::NUMERIC_START, Expr::EXPR_START);
  for (; i < n; ++i) {
    Expr::Kind kind = kinds[i];
    EXPECT_GE(kind, Expr::NUMERIC_START);
    EXPECT_LE(kind, Expr::NUMERIC_END);
    for (int j = i + 1; j < n; ++j)
      EXPECT_NE(kind, kinds[j]);  // Check if all different.
    if (kind != Expr::CONSTANT)
      EXPECT_LT(kind, Expr::LOGICAL_START);
  }
  EXPECT_EQ(i, n);
}

TEST_F(ExprTest, LogicalKinds) {
  const Expr::Kind kinds[] = {
      Expr::CONSTANT,
      Expr::RELATIONAL,
      Expr::NOT,
      Expr::BINARY_LOGICAL,
      Expr::IMPLICATION,
      Expr::ITERATED_LOGICAL,
      Expr::ALLDIFF
  };
  int i = 0, n = sizeof(kinds) / sizeof(*kinds);
  EXPECT_GT(n, 0);
  EXPECT_LT(Expr::LOGICAL_END, Expr::EXPR_END);
  for (; i < n; ++i) {
    Expr::Kind kind = kinds[i];
    EXPECT_GE(kind, Expr::LOGICAL_START);
    EXPECT_LE(kind, Expr::LOGICAL_END);
    for (int j = i + 1; j < n; ++j)
      EXPECT_NE(kind, kinds[j]);  // Check if all different.
    if (kind != Expr::CONSTANT)
      EXPECT_GT(kind, Expr::NUMERIC_END);
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
    expr raw = {reinterpret_cast<efunc*>(OPMINUS)};
    Expr e(::MakeExpr(&raw));
    EXPECT_EQ(OPMINUS, e.opcode());
  }
  {
    expr raw = {reinterpret_cast<efunc*>(OPOR)};
    Expr e(::MakeExpr(&raw));
    EXPECT_EQ(OPOR, e.opcode());
  }
}

TEST_F(ExprTest, SafeBool) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = {reinterpret_cast<efunc*>(42)};
  Expr e2(::MakeExpr(&raw2));
  EXPECT_TRUE(e2);
}

struct OpInfo {
  int code;
  const char *str;
  Expr::Kind kind;
};

const OpInfo OP_INFO[] = {
  {-1, "unknown"},
  {OPPLUS,  "+",    Expr::BINARY},
  {OPMINUS, "-",    Expr::BINARY},
  {OPMULT,  "*",    Expr::BINARY},
  {OPDIV,   "/",    Expr::BINARY},
  {OPREM,   "mod",  Expr::BINARY},
  {OPPOW,   "^",    Expr::BINARY},
  {OPLESS,  "less", Expr::BINARY},
  { 7, "unknown"},
  { 8, "unknown"},
  { 9, "unknown"},
  {10, "unknown"},
  {MINLIST,  "min",     Expr::VARARG},
  {MAXLIST,  "max",     Expr::VARARG},
  {FLOOR,    "floor",   Expr::UNARY},
  {CEIL,     "ceil",    Expr::UNARY},
  {ABS,      "abs",     Expr::UNARY},
  {OPUMINUS, "unary -", Expr::UNARY},
  {17, "unknown"},
  {18, "unknown"},
  {19, "unknown"},
  {OPOR,  "||", Expr::BINARY_LOGICAL},
  {OPAND, "&&", Expr::BINARY_LOGICAL},
  {LT,    "<",  Expr::RELATIONAL},
  {LE,    "<=", Expr::RELATIONAL},
  {EQ,    "=",  Expr::RELATIONAL},
  {25, "unknown"},
  {26, "unknown"},
  {27, "unknown"},
  {GE, ">=", Expr::RELATIONAL},
  {GT, ">",  Expr::RELATIONAL},
  {NE, "!=", Expr::RELATIONAL},
  {31, "unknown"},
  {32, "unknown"},
  {33, "unknown"},
  {OPNOT, "!", Expr::NOT},
  {OPIFnl, "if", Expr::IF},
  {36, "unknown"},
  {OP_tanh,  "tanh",  Expr::UNARY},
  {OP_tan,   "tan",   Expr::UNARY},
  {OP_sqrt,  "sqrt",  Expr::UNARY},
  {OP_sinh,  "sinh",  Expr::UNARY},
  {OP_sin,   "sin",   Expr::UNARY},
  {OP_log10, "log10", Expr::UNARY},
  {OP_log,   "log",   Expr::UNARY},
  {OP_exp,   "exp",   Expr::UNARY},
  {OP_cosh,  "cosh",  Expr::UNARY},
  {OP_cos,   "cos",   Expr::UNARY},
  {OP_atanh, "atanh", Expr::UNARY},
  {OP_atan2, "atan2", Expr::BINARY},
  {OP_atan,  "atan",  Expr::UNARY},
  {OP_asinh, "asinh", Expr::UNARY},
  {OP_asin,  "asin",  Expr::UNARY},
  {OP_acosh, "acosh", Expr::UNARY},
  {OP_acos,  "acos",  Expr::UNARY},
  {OPSUMLIST, "sum",  Expr::SUM},
  {OPintDIV,  "div",  Expr::BINARY},
  {OPprecision, "precision", Expr::BINARY},
  {OPround,     "round",     Expr::BINARY},
  {OPtrunc,     "trunc",     Expr::BINARY},
  {OPCOUNT,     "count",           Expr::COUNT},
  {OPNUMBEROF,  "numberof",        Expr::NUMBEROF},
  {OPNUMBEROFs, "string numberof", Expr::UNKNOWN},
  {OPATLEAST, "atleast", Expr::LOGICAL_COUNT},
  {OPATMOST,  "atmost",  Expr::LOGICAL_COUNT},
  {OPPLTERM, "pl term", Expr::PLTERM},
  {OPIFSYM,  "string if-then-else", Expr::UNKNOWN},
  {OPEXACTLY,    "exactly",  Expr::LOGICAL_COUNT},
  {OPNOTATLEAST, "!atleast", Expr::LOGICAL_COUNT},
  {OPNOTATMOST,  "!atmost",  Expr::LOGICAL_COUNT},
  {OPNOTEXACTLY, "!exactly", Expr::LOGICAL_COUNT},
  {ANDLIST, "forall", Expr::ITERATED_LOGICAL},
  {ORLIST,  "exists", Expr::ITERATED_LOGICAL},
  {OPIMPELSE, "==>", Expr::IMPLICATION},
  {OP_IFF, "<==>", Expr::BINARY_LOGICAL},
  {OPALLDIFF, "alldiff", Expr::ALLDIFF},
  {OP1POW, "^",  Expr::BINARY},
  {OP2POW, "^2", Expr::UNARY},
  {OPCPOW, "^",  Expr::BINARY},
  {OPFUNCALL, "function call", Expr::CALL},
  {OPNUM, "number", Expr::CONSTANT},
  {OPHOL, "string", Expr::STRING},
  {OPVARVAL, "variable", Expr::VARIABLE},
  {N_OPS, "unknown"},
  {777,   "unknown"}
};

template <typename ExprT>
void TestAssertInCreate(int opcode) {
  expr raw = {reinterpret_cast<efunc*>(opcode)};
  EXPECT_DEBUG_DEATH(MakeExpr<ExprT>(&raw);, "Assertion");  // NOLINT(*)
}

template <typename ExprT>
int CheckExpr(Expr::Kind start, Expr::Kind end = Expr::UNKNOWN,
    int bad_opcode = OPPLTERM) {
  if (end == Expr::UNKNOWN)
    end = start;
  {
    // Check default ctor.
    ExprT e;
    EXPECT_FALSE(e);
  }
  TestAssertInCreate<ExprT>(bad_opcode);
  int expr_count = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  EXPECT_EQ(N_OPS + 3, size);
  for (int i = 0; i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int opcode = info.code;
    const char *opstr = info.str;
    expr raw = {reinterpret_cast<efunc*>(opcode)};
    bool is_this_kind = info.kind >= start && info.kind <= end;
    if (info.kind != Expr::UNKNOWN) {
      Expr e(::MakeExpr(&raw));
      EXPECT_EQ(is_this_kind, ampl::internal::Is<ExprT>(e));
      bool cast_result = Cast<ExprT>(e);
      EXPECT_EQ(is_this_kind, cast_result);
    }
    if (!is_this_kind) continue;
    ExprT e(MakeExpr<ExprT>(&raw));
    EXPECT_EQ(opcode, e.opcode());
    EXPECT_STREQ(opstr, e.opstr());
    ++expr_count;
  }
  EXPECT_GT(expr_count, 0);
  return expr_count;
}

TEST_F(ExprTest, Expr) {
  EXPECT_EQ(66, CheckExpr<Expr>(Expr::EXPR_START, Expr::EXPR_END, -1));
  TestAssertInCreate<Expr>(7);
  TestAssertInCreate<Expr>(N_OPS);
  TestAssertInCreate<Expr>(777);
}

// Test if Expr::Create() uses internal::Is() to check whether expression
// is of specified type. This allows testing Is() functions instead of doing
// expensive death tests. Is() is specialized for TestExpr to accept unary
// and binary but not other expression kinds.
TEST_F(ExprTest, CreateUsesIs) {
  expr raw1 = {reinterpret_cast<efunc*>(OPPLUS)};  // binary
  ::MakeExpr<TestExpr>(&raw1);
  expr raw2 = {reinterpret_cast<efunc*>(OPUMINUS)};  // unary
  ::MakeExpr<TestExpr>(&raw2);
  TestAssertInCreate<TestExpr>(OPPLTERM);  // neither
}

TEST_F(ExprTest, EqualityOperator) {
  expr raw1 = {}, raw2 = {};
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
  NumericExpr e = MakeUnary(OPUMINUS, MakeVariable(0));
  EXPECT_TRUE(Equal(e, MakeUnary(OPUMINUS, MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeVariable(0)));
  EXPECT_FALSE(Equal(e, MakeUnary(FLOOR, MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeUnary(OPUMINUS, MakeVariable(1))));
}

TEST_F(ExprTest, EqualBinaryExpr) {
  NumericExpr e = MakeBinary(OPPLUS, MakeVariable(0), MakeConst(42));
  EXPECT_TRUE(Equal(e, MakeBinary(OPPLUS, MakeVariable(0), MakeConst(42))));
  EXPECT_FALSE(Equal(e, MakeBinary(OPMINUS, MakeVariable(0), MakeConst(42))));
  EXPECT_FALSE(Equal(e, MakeBinary(OPPLUS, MakeConst(42), MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeBinary(OPPLUS, MakeVariable(0), MakeConst(0))));
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
  NumericExpr e = builder.MakeVarArg(MINLIST, args1);
  EXPECT_TRUE(Equal(e, builder.MakeVarArg(MINLIST, args2)));
  NumericExpr args3[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_FALSE(Equal(e, builder.MakeVarArg(MINLIST, args3)));
  EXPECT_FALSE(Equal(builder.MakeVarArg(MINLIST, args3),
                     builder.MakeVarArg(MINLIST, args1)));
  EXPECT_FALSE(Equal(e, builder.MakeVarArg(MAXLIST, args1)));
  NumericExpr args4[] = {MakeVariable(0), MakeVariable(1), MakeConst(0)};
  EXPECT_FALSE(Equal(e, builder.MakeVarArg(MINLIST, args4)));
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
  EXPECT_FALSE(Equal(e, l1));
}

TEST_F(ExprTest, EqualStringLiteral) {
  StringLiteral s1 = builder.MakeStringLiteral("abc");
  StringLiteral s2 = builder.MakeStringLiteral("abc");
  EXPECT_NE(s1.value(), s2.value());
  EXPECT_TRUE(Equal(s1, s2));
  EXPECT_FALSE(Equal(s1, builder.MakeStringLiteral("def")));
}

TEST_F(ExprTest, NumericExpr) {
  EXPECT_EQ(45,
      CheckExpr<NumericExpr>(Expr::NUMERIC_START, Expr::NUMERIC_END, OPNOT));
}

TEST_F(ExprTest, LogicalExpr) {
  EXPECT_EQ(21,
      CheckExpr<LogicalExpr>(Expr::LOGICAL_START, Expr::LOGICAL_END));
}

TEST_F(ExprTest, NumericConstant) {
  EXPECT_EQ(1, CheckExpr<NumericConstant>(Expr::CONSTANT));
  ampl::NumericConstant expr = builder.MakeNumericConstant(42);
  EXPECT_EQ(OPNUM, expr.opcode());
  EXPECT_EQ(42, expr.value());
}

TEST_F(ExprTest, Variable) {
  EXPECT_EQ(1, CheckExpr<Variable>(Expr::VARIABLE));
  ampl::Variable var = builder.MakeVariable(0);
  EXPECT_EQ(OPVARVAL, var.opcode());
  EXPECT_EQ(0, var.index());
  var = builder.MakeVariable(9);
  EXPECT_EQ(9, var.index());
  EXPECT_DEBUG_DEATH(builder.MakeVariable(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(builder.MakeVariable(NUM_VARS);, "Assertion");  // NOLINT(*)
}

TEST_F(ExprTest, UnaryExpr) {
  const int opcodes[] = {
      FLOOR, CEIL, ABS, OPUMINUS, OP_tanh, OP_tan, OP_sqrt,
      OP_sinh, OP_sin, OP_log10, OP_log, OP_exp, OP_cosh, OP_cos,
      OP_atanh, OP_atan, OP_asinh, OP_asin, OP_acosh, OP_acos, OP2POW
  };
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes, CheckExpr<UnaryExpr>(Expr::UNARY));
  for (std::size_t i = 0; i < num_opcodes; ++i) {
    ampl::UnaryExpr expr = builder.MakeUnary(opcodes[i], n1);
    EXPECT_EQ(opcodes[i], expr.opcode());
    EXPECT_EQ(n1, expr.arg());
  }
  EXPECT_THROW_MSG(builder.MakeUnary(OPPLUS, n1), Error,
    fmt::format("invalid unary expression code {}", OPPLUS));
}

#define EXPECT_BINARY(expr, expected_opcode, expected_lhs, expected_rhs) \
  EXPECT_EQ(expected_opcode, expr.opcode()); \
  EXPECT_EQ(expected_lhs, expr.lhs()); \
  EXPECT_EQ(expected_rhs, expr.rhs())

TEST_F(ExprTest, BinaryExpr) {
  const int opcodes[] = {
      OPPLUS, OPMINUS, OPMULT, OPDIV, OPREM, OPPOW, OPLESS, OP_atan2,
      OPintDIV, OPprecision, OPround, OPtrunc, OP1POW, OPCPOW
  };
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes, CheckExpr<BinaryExpr>(Expr::BINARY));
  for (std::size_t i = 0; i < num_opcodes; ++i) {
    ampl::BinaryExpr expr = builder.MakeBinary(opcodes[i], n1, n2);
    EXPECT_BINARY(expr, opcodes[i], n1, n2);
  }
  EXPECT_THROW_MSG(builder.MakeBinary(OPUMINUS, n1, n2), Error,
    fmt::format("invalid binary expression code {}", OPUMINUS));
}

TEST_F(ExprTest, IfExpr) {
  EXPECT_EQ(1, CheckExpr<IfExpr>(Expr::IF));
  ampl::IfExpr expr = builder.MakeIf(l1, n1, n2);
  EXPECT_EQ(OPIFnl, expr.opcode());
  EXPECT_EQ(l1, expr.condition());
  EXPECT_EQ(n1, expr.true_expr());
  EXPECT_EQ(n2, expr.false_expr());
}

TEST_F(ExprTest, PiecewiseLinearExpr) {
  EXPECT_EQ(1,
      CheckExpr<PiecewiseLinearExpr>(Expr::PLTERM, Expr::PLTERM, OPPLUS));
  enum { NUM_BREAKPOINTS = 2 };
  double breakpoints[NUM_BREAKPOINTS] = { 11, 22 };
  double slopes[NUM_BREAKPOINTS + 1] = {33, 44, 55};
  ampl::Variable var = builder.MakeVariable(2);
  ampl::PiecewiseLinearExpr expr = builder.MakePiecewiseLinear(
      NUM_BREAKPOINTS, breakpoints, slopes, var);
  EXPECT_EQ(OPPLTERM, expr.opcode());
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
  EXPECT_EQ(1, CheckExpr<CallExpr>(Expr::CALL));
  enum {NUM_ARGS = 3};
  Function f = builder.AddFunction("foo", TestFunc, NUM_ARGS, Function::SYMBOLIC);
  const Expr args[NUM_ARGS] = {n1, n2, builder.MakeStringLiteral("abc")};
  CallExpr expr = builder.MakeCall(f, args);
  EXPECT_EQ(OPFUNCALL, expr.opcode());
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
  const int opcodes[] = {MINLIST, MAXLIST};
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes, CheckExpr<VarArgExpr>(Expr::VARARG));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  for (size_t i = 0, n = sizeof(opcodes) / sizeof(*opcodes); i < n; ++i) {
    VarArgExpr expr = builder.MakeVarArg(opcodes[i], args);
    EXPECT_EQ(opcodes[i], expr.opcode());
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
  EXPECT_THROW_MSG(builder.MakeVarArg(OPUMINUS, args), Error,
      fmt::format("invalid vararg expression code {}", OPUMINUS));
}

TEST_F(ExprTest, SumExpr) {
  EXPECT_EQ(1, CheckExpr<SumExpr>(Expr::SUM));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  ampl::SumExpr expr = builder.MakeSum(args);
  EXPECT_EQ(OPSUMLIST, expr.opcode());
  int arg_index = 0;
  for (ampl::SumExpr::iterator
      i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
    EXPECT_EQ(args[arg_index], *i);
  }
}

TEST_F(ExprTest, CountExpr) {
  EXPECT_EQ(1, CheckExpr<CountExpr>(Expr::COUNT));
  enum {NUM_ARGS = 2};
  LogicalExpr args[NUM_ARGS] = {l1, l0};
  ampl::CountExpr expr = builder.MakeCount(args);
  EXPECT_EQ(OPCOUNT, expr.opcode());
  int arg_index = 0;
  for (ampl::CountExpr::iterator
      i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
    EXPECT_EQ(args[arg_index], *i);
  }
}

TEST_F(ExprTest, NumberOfExpr) {
  EXPECT_EQ(1, CheckExpr<NumberOfExpr>(Expr::NUMBEROF));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  ampl::NumberOfExpr expr = builder.MakeNumberOf(args);
  EXPECT_EQ(OPNUMBEROF, expr.opcode());
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  for (int i = 0; i < NUM_ARGS; ++i)
    EXPECT_EQ(args[i], expr[i]);
#ifndef NDEBUG
  EXPECT_DEBUG_DEATH(
      builder.MakeNumberOf(MakeArrayRef(args, 0));, "Assertion");  // NOLINT(*)
#endif
}

TEST_F(ExprTest, LogicalConstant) {
  EXPECT_EQ(1, CheckExpr<LogicalConstant>(Expr::CONSTANT));
  ampl::LogicalConstant expr = builder.MakeLogicalConstant(true);
  EXPECT_EQ(OPNUM, expr.opcode());
  EXPECT_TRUE(expr.value());
  EXPECT_FALSE(builder.MakeLogicalConstant(false).value());
}

TEST_F(ExprTest, NotExpr) {
  EXPECT_EQ(1, CheckExpr<NotExpr>(Expr::NOT));
  ampl::NotExpr expr = builder.MakeNot(l1);
  EXPECT_EQ(OPNOT, expr.opcode());
  EXPECT_EQ(l1, expr.arg());
}

TEST_F(ExprTest, BinaryLogicalExpr) {
  const int opcodes[] = {OPOR, OPAND, OP_IFF};
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes, CheckExpr<BinaryLogicalExpr>(Expr::BINARY_LOGICAL));
  for (size_t i = 0; i < num_opcodes; ++i) {
    ampl::BinaryLogicalExpr expr =
        builder.MakeBinaryLogical(opcodes[i], l1, l0);
    EXPECT_BINARY(expr, opcodes[i], l1, l0);
  }
  EXPECT_THROW_MSG(builder.MakeBinaryLogical(OPUMINUS, l1, l0), Error,
    fmt::format("invalid binary expression code {}", OPUMINUS));
}

TEST_F(ExprTest, RelationalExpr) {
  const int opcodes[] = {LT, LE, EQ, GE, GT, NE};
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes, CheckExpr<RelationalExpr>(Expr::RELATIONAL));
  for (size_t i = 0; i < num_opcodes; ++i) {
    ampl::RelationalExpr expr = builder.MakeRelational(opcodes[i], n1, n2);
    EXPECT_BINARY(expr, opcodes[i], n1, n2);
  }
  EXPECT_THROW_MSG(builder.MakeRelational(OPUMINUS, n1, n2), Error,
    fmt::format("invalid binary expression code {}", OPUMINUS));
}

TEST_F(ExprTest, LogicalCountExpr) {
  const int opcodes[] = {
    OPATLEAST, OPATMOST, OPEXACTLY, OPNOTATLEAST, OPNOTATMOST, OPNOTEXACTLY
  };
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes, CheckExpr<LogicalCountExpr>(Expr::LOGICAL_COUNT));
  ampl::CountExpr count = builder.MakeCount(MakeArrayRef(&l1, 1));
  for (size_t i = 0; i < num_opcodes; ++i) {
    ampl::LogicalCountExpr expr =
        builder.MakeLogicalCount(opcodes[i], n1, count);
    EXPECT_BINARY(expr, opcodes[i], n1, count);
  }
  EXPECT_THROW_MSG(builder.MakeLogicalCount(OPUMINUS, n1, count), Error,
    fmt::format("invalid binary expression code {}", OPUMINUS));
}

TEST_F(ExprTest, ImplicationExpr) {
  EXPECT_EQ(1, CheckExpr<ImplicationExpr>(Expr::IMPLICATION));
  LogicalExpr condition = builder.MakeLogicalConstant(true);
  ampl::ImplicationExpr expr = builder.MakeImplication(condition, l0, l1);
  EXPECT_EQ(OPIMPELSE, expr.opcode());
  EXPECT_EQ(condition, expr.condition());
  EXPECT_EQ(l0, expr.true_expr());
  EXPECT_EQ(l1, expr.false_expr());
}

TEST_F(ExprTest, IteratedLogicalExpr) {
  const int opcodes[] = {ORLIST, ANDLIST};
  std::size_t num_opcodes = sizeof(opcodes) / sizeof(*opcodes);
  EXPECT_EQ(num_opcodes,
    CheckExpr<IteratedLogicalExpr>(Expr::ITERATED_LOGICAL));
  enum {NUM_ARGS = 3};
  LogicalExpr args[NUM_ARGS] = {l0, l1, builder.MakeLogicalConstant(false)};
  for (size_t i = 0; i < num_opcodes; ++i) {
    IteratedLogicalExpr expr = builder.MakeIteratedLogical(opcodes[i], args);
    EXPECT_EQ(opcodes[i], expr.opcode());
    int arg_index = 0;
    for (IteratedLogicalExpr::iterator
        i = expr.begin(), end = expr.end(); i != end; ++i, ++arg_index) {
      EXPECT_EQ(args[arg_index], *i);
    }
  }
  EXPECT_THROW_MSG(builder.MakeIteratedLogical(OPUMINUS, args), Error,
        fmt::format("invalid iterated logical expression code {}", OPUMINUS));
}

TEST_F(ExprTest, AllDiffExpr) {
  EXPECT_EQ(1, CheckExpr<AllDiffExpr>(Expr::ALLDIFF));
  enum {NUM_ARGS = 3};
  NumericExpr args[NUM_ARGS] = {n1, n2, builder.MakeNumericConstant(3)};
  AllDiffExpr expr = builder.MakeAllDiff(args);
  EXPECT_EQ(NUM_ARGS, expr.num_args());
  int index = 0;
  for (AllDiffExpr::iterator
       i = expr.begin(), end = expr.end(); i != end; ++i, ++index) {
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index], expr[index]);
  }
  EXPECT_EQ(NUM_ARGS, index);
}

TEST_F(ExprTest, StringLiteral) {
  EXPECT_EQ(1, CheckExpr<StringLiteral>(Expr::STRING));
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
  TestResult VisitPlus(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitMinus(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitMult(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitDiv(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitRem(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitPow(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitNumericLess(BinaryExpr e) { return MakeResult(e); }
  TestResult VisitMin(VarArgExpr e) { return MakeResult(e); }
  TestResult VisitMax(VarArgExpr e) { return MakeResult(e); }
  TestResult VisitFloor(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitCeil(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitAbs(UnaryExpr e) { return MakeResult(e); }
  TestResult VisitUnaryMinus(UnaryExpr e) { return MakeResult(e); }
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
  TestLResult VisitLess(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitLessEqual(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitEqual(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGreaterEqual(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGreater(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitNotEqual(RelationalExpr e) { return MakeResult(e); }
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
  TestLResult VisitAllDiff(AllDiffExpr e) { return MakeResult(e); }
  TestLResult VisitLogicalConstant(LogicalConstant e) { return MakeResult(e); }
};

TEST_F(ExprTest, ExprVisitorHandlesAll) {
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  for (int i = 0; i < size; ++i) {
    FullTestVisitor v;
    const OpInfo &info = OP_INFO[i];
    if (info.kind == Expr::UNKNOWN || info.kind == Expr::STRING) continue;
    expr raw = {reinterpret_cast<efunc*>(info.code)};
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
    if (info.kind == Expr::UNKNOWN || info.kind == Expr::STRING) continue;
    expr raw = {reinterpret_cast<efunc*>(info.code)};
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

TEST_F(ExprTest, ExprVisitorInvalidThrows) {
  expr raw = {reinterpret_cast<efunc*>(OPNUM)};
  NumericExpr ne(::MakeExpr<NumericExpr>(&raw));
  LogicalExpr le(::MakeExpr<LogicalExpr>(&raw));
  raw.op = reinterpret_cast<efunc*>(-1);
  EXPECT_THROW(NullVisitor().Visit(ne), InvalidNumericExprError);
  EXPECT_THROW(NullVisitor().Visit(le), InvalidLogicalExprError);
}

struct TestConverter : ampl::ExprConverter<TestConverter, void, TestLResult> {
  TestLResult VisitLess(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitLessEqual(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitEqual(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGreaterEqual(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitGreater(RelationalExpr e) { return MakeResult(e); }
  TestLResult VisitNotEqual(RelationalExpr e) { return MakeResult(e); }
};

void CheckConversion(int from_opcode, int to_opcode) {
  expr lhs = {reinterpret_cast<efunc*>(OPNUM)}, rhs = lhs;
  expr raw = {reinterpret_cast<efunc*>(from_opcode)};
  raw.L.e = &lhs;
  raw.R.e = &rhs;
  TestConverter converter;
  RelationalExpr expr =
      Cast<RelationalExpr>(converter.Visit(::MakeExpr<LogicalExpr>(&raw)).expr);
  EXPECT_EQ(to_opcode, expr.opcode());
  EXPECT_EQ(::MakeExpr(&lhs), expr.lhs());
  EXPECT_EQ(::MakeExpr(&rhs), expr.rhs());
}

TEST_F(ExprTest, ConvertLogicalCountToRelational) {
  CheckConversion(OPATLEAST, LE);
  CheckConversion(OPATMOST, GE);
  CheckConversion(OPEXACTLY, EQ);
  CheckConversion(OPNOTATLEAST, GT);
  CheckConversion(OPNOTATMOST, LT);
  CheckConversion(OPNOTEXACTLY, NE);
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
  TestGrad g2 = {0};
  TestGrad g1 = {&g2};
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
  ampl::NumberOfMap<Var, CreateVar> map((CreateVar()));
  EXPECT_TRUE(map.begin() == map.end());
  NumericExpr args1[] = {MakeConst(11), MakeVariable(0)};
  NumberOfExpr e1 = builder.MakeNumberOf(args1);
  NumericExpr args2[] = {MakeConst(22), MakeVariable(1)};
  NumberOfExpr e2 = builder.MakeNumberOf(args2);
  map.Add(11, e1);
  map.Add(22, e2);
  NumericExpr args3[] = {MakeConst(33), MakeVariable(0)};
  map.Add(33, builder.MakeNumberOf(args3));
  ampl::NumberOfMap<Var, CreateVar>::iterator i = map.begin();
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
  fmt::Writer w;
  WriteExpr(w, ampl::LinearObjExpr(), expr);
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
  CHECK_WRITE("-x1", MakeUnary(OPUMINUS, x1));
  CHECK_WRITE("x1 ^ 2", MakeUnary(OP2POW, x1));
  int count = 0;
  for (int i = 0, size = sizeof(OP_INFO) / sizeof(*OP_INFO); i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int code = info.code;
    if (info.kind != Expr::UNARY || code == OPUMINUS || code == OP2POW)
      continue;
    CHECK_WRITE(fmt::format("{}(x1)", info.str), MakeUnary(code, x1));
    ++count;
  }
  EXPECT_EQ(19, count);
}

TEST_F(ExprTest, WriteBinaryExpr) {
  Variable x1 = MakeVariable(0);
  NumericConstant n42 = MakeConst(42);
  CHECK_WRITE("x1 + 42", MakeBinary(OPPLUS, x1, n42));
  CHECK_WRITE("x1 - 42", MakeBinary(OPMINUS, x1, n42));
  CHECK_WRITE("x1 * 42", MakeBinary(OPMULT, x1, n42));
  CHECK_WRITE("x1 / 42", MakeBinary(OPDIV, x1, n42));
  CHECK_WRITE("x1 mod 42", MakeBinary(OPREM, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(OPPOW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(OPPOW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(OP1POW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(OPCPOW, x1, n42));
  CHECK_WRITE("x1 less 42", MakeBinary(OPLESS, x1, n42));
  CHECK_WRITE("x1 div 42", MakeBinary(OPintDIV, x1, n42));
}

TEST_F(ExprTest, WriteBinaryFunc) {
  auto x1 = MakeVariable(0);
  auto n42 = MakeConst(42);
  CHECK_WRITE("atan2(x1, 42)", MakeBinary(OP_atan2, x1, n42));
  CHECK_WRITE("precision(x1, 42)", MakeBinary(OPprecision, x1, n42));
  CHECK_WRITE("round(x1, 42)", MakeBinary(OPround, x1, n42));
  CHECK_WRITE("trunc(x1, 42)", MakeBinary(OPtrunc, x1, n42));
}

TEST_F(ExprTest, WriteIfExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if x1 = 0 then 1",
      MakeIf(MakeRelational(EQ, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 0 else 1",
      MakeIf(MakeRelational(EQ, MakeVariable(0), n0), n0, n1));
}

TEST_F(ExprTest, WritePiecewiseLinearExpr) {
  double breaks[] = {5, 10};
  double slopes[] = {-1, 0, 1};
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43",
              builder.MakePiecewiseLinear(2, breaks, slopes, MakeVariable(42)));
}

TEST_F(ExprTest, WriteCallExpr) {
  ampl::Function f = builder.AddFunction("foo", TestFunc, -1);
  Expr args[] = {
      MakeConst(3),
      MakeBinary(OPPLUS, MakeVariable(0), MakeConst(5)),
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
  CHECK_WRITE("min(x1, x2, 42)", builder.MakeVarArg(MINLIST, args));
  CHECK_WRITE("max(x1, x2, 42)", builder.MakeVarArg(MAXLIST, args));
}

TEST_F(ExprTest, WriteSumExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  CHECK_WRITE("/* sum */ (x1 + x2 + 42)", builder.MakeSum(args));
  NumericExpr args2[] = {
    MakeBinary(OPPLUS, MakeVariable(0), MakeVariable(1)), MakeConst(42)
  };
  CHECK_WRITE("/* sum */ ((x1 + x2) + 42)", builder.MakeSum(args2));
}

TEST_F(ExprTest, WriteCountExpr) {
  LogicalExpr args[] = {
    MakeRelational(EQ, MakeVariable(0), MakeConst(0)), l1, l0
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
      MakeIf(builder.MakeNot(MakeRelational(EQ, MakeVariable(0), n0)), n1, n0));
}

TEST_F(ExprTest, WriteBinaryLogicalExpr) {
  auto e1 = MakeRelational(GT, MakeVariable(0), MakeConst(0));
  auto e2 = MakeRelational(LT, MakeVariable(0), MakeConst(10));
  CHECK_WRITE("if x1 > 0 || x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(OPOR, e1, e2), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 > 0 && x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(OPAND, e1, e2), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 > 0 <==> x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(OP_IFF, e1, e2), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteRelationalExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if x1 < 0 then 1",
      MakeIf(MakeRelational(LT, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 <= 0 then 1",
      MakeIf(MakeRelational(LE, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 1",
      MakeIf(MakeRelational(EQ, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 >= 0 then 1",
      MakeIf(MakeRelational(GE, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 > 0 then 1",
      MakeIf(MakeRelational(GT, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 != 0 then 1",
      MakeIf(MakeRelational(NE, MakeVariable(0), n0), n1, n0));
}

TEST_F(ExprTest, WriteLogicalCountExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1), value = MakeConst(42);
  LogicalExpr args[] = {
    MakeRelational(EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  ampl::CountExpr count = builder.MakeCount(args);
  CHECK_WRITE("if atleast 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(OPATLEAST, value, count), n1, n0));
  CHECK_WRITE("if atmost 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(OPATMOST, value, count), n1, n0));
  CHECK_WRITE("if exactly 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(OPEXACTLY, value, count), n1, n0));
  CHECK_WRITE("if !atleast 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(OPNOTATLEAST, value, count), n1, n0));
  CHECK_WRITE("if !atmost 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(OPNOTATMOST, value, count), n1, n0));
  CHECK_WRITE("if !exactly 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(OPNOTEXACTLY, value, count), n1, n0));
}

TEST_F(ExprTest, WriteImplicationExpr) {
  auto e1 = MakeRelational(EQ, MakeVariable(0), MakeConst(0));
  CHECK_WRITE("if x1 = 0 ==> 1 then 1",
      MakeIf(MakeImplication(e1, l1, l0), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 = 0 ==> 0 else 1 then 1",
      MakeIf(MakeImplication(e1, l0, l1), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteIteratedLogicalExpr) {
  LogicalExpr args[] = {
    MakeRelational(EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  CHECK_WRITE("if /* forall */ (x1 = 0 && 1 && 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ANDLIST, args),
             MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if /* exists */ (x1 = 0 || 1 || 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ORLIST, args),
             MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteAllDiffExpr) {
  NumericExpr args[] = {MakeConst(42), MakeConst(43), MakeConst(44)};
  CHECK_WRITE("if alldiff(42, 43, 44) then 1",
      MakeIf(builder.MakeAllDiff(args), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprTest, WriteStringLiteral) {
  Expr args[] = {builder.MakeStringLiteral("abc")};
  ampl::Function f = builder.AddFunction("f", TestFunc, 1, Function::SYMBOLIC);
  CHECK_WRITE("f('abc')", builder.MakeCall(f, args));
  args[0] = builder.MakeStringLiteral("ab'c");
  CHECK_WRITE("f('ab''c')", builder.MakeCall(f, args));
  args[0] = builder.MakeStringLiteral("ab\nc");
  CHECK_WRITE("f('ab\\\nc')", builder.MakeCall(f, args));
}

TEST_F(ExprTest, UnaryExprPrecedence) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("--x1", MakeUnary(OPUMINUS, MakeUnary(OPUMINUS, x1)));
  CHECK_WRITE("-(x1 ^ x1)", MakeUnary(OPUMINUS, MakeBinary(OPPOW, x1, x1)));
}

TEST_F(ExprTest, UnaryFuncPrecedence) {
  auto x1 = MakeVariable(0);
  int count = 0;
  for (int i = 0, size = sizeof(OP_INFO) / sizeof(*OP_INFO); i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int code = info.code;
    if (info.kind != Expr::UNARY || code == OPUMINUS || code == OP2POW)
      continue;
    CHECK_WRITE(fmt::format("{0}({0}(x1))", info.str),
        MakeUnary(code, MakeUnary(code, x1)));
    CHECK_WRITE(fmt::format("{0}(x1 + x1)", info.str),
        MakeUnary(code, MakeBinary(OPPLUS, x1, x1)));
    ++count;
  }
  EXPECT_EQ(19, count);
}

TEST_F(ExprTest, Pow2Precedence) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("(x1 ^ 2) ^ 2", MakeUnary(OP2POW, MakeUnary(OP2POW, x1)));
  CHECK_WRITE("(x1 * x1) ^ 2", MakeUnary(OP2POW, MakeBinary(OPMULT, x1, x1)));
}

TEST_F(ExprTest, AdditiveExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 + x2 + x3",
      MakeBinary(OPPLUS, MakeBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("x1 + x2 - x3",
      MakeBinary(OPMINUS, MakeBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("x1 + x2 less x3",
      MakeBinary(OPLESS, MakeBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("x1 + (x2 + x3)",
      MakeBinary(OPPLUS, x1, MakeBinary(OPPLUS, x2, x3)));
  CHECK_WRITE("(x1 + x2) * x3",
      MakeBinary(OPMULT, MakeBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("if 1 then x1 else x2 + x3",
      MakeIf(l1, x1, MakeBinary(OPPLUS, x2, x3)));
}

TEST_F(ExprTest, MultiplicativeExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 * x2 * x3",
      MakeBinary(OPMULT, MakeBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * x2 / x3",
      MakeBinary(OPDIV, MakeBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * x2 div x3",
      MakeBinary(OPintDIV, MakeBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * x2 mod x3",
      MakeBinary(OPREM, MakeBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * (x2 * x3)",
      MakeBinary(OPMULT, x1, MakeBinary(OPMULT, x2, x3)));
  CHECK_WRITE("(x1 * x2) ^ x3",
      MakeBinary(OPPOW, MakeBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("(x1 + x2) * x3",
      MakeBinary(OPMULT, MakeBinary(OPPLUS, x1, x2), x3));
}

TEST_F(ExprTest, ExponentiationExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 ^ x2 ^ x3",
      MakeBinary(OPPOW, x1, MakeBinary(OPPOW, x2, x3)));
  CHECK_WRITE("x1 ^ x2 ^ 3",
      MakeBinary(OPPOW, x1, MakeBinary(OP1POW, x2, MakeConst(3))));
  CHECK_WRITE("x1 ^ 3 ^ x2",
      MakeBinary(OPPOW, x1, MakeBinary(OPCPOW, MakeConst(3), x2)));
  CHECK_WRITE("(x1 ^ 2) ^ 3",
      MakeBinary(OP1POW, MakeBinary(OP1POW, x1, MakeConst(2)), MakeConst(3)));
  CHECK_WRITE("-x1 ^ -x2",
      MakeBinary(OPPOW, MakeUnary(OPUMINUS, x1), MakeUnary(OPUMINUS, x2)));
  CHECK_WRITE("x1 ^ (x2 * x3)",
      MakeBinary(OPPOW, x1, MakeBinary(OPMULT, x2, x3)));
}

TEST_F(ExprTest, BinaryFuncPrecedence) {
  auto x1 = MakeVariable(0);
  auto e = MakeBinary(OPPLUS, x1, x1);
  CHECK_WRITE("atan2(atan2(x1, x1), x1 + x1)",
      MakeBinary(OP_atan2, MakeBinary(OP_atan2, x1, x1), e));
  CHECK_WRITE("precision(precision(x1, x1), x1 + x1)",
      MakeBinary(OPprecision, MakeBinary(OPprecision, x1, x1), e));
  CHECK_WRITE("round(round(x1, x1), x1 + x1)",
      MakeBinary(OPround, MakeBinary(OPround, x1, x1), e));
  CHECK_WRITE("trunc(trunc(x1, x1), x1 + x1)",
      MakeBinary(OPtrunc, MakeBinary(OPtrunc, x1, x1), e));
}

TEST_F(ExprTest, IfExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), n2 = MakeConst(2);
  LogicalExpr e = MakeBinaryLogical(OPOR, l0, l1);
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1",
      MakeIf(e, MakeIf(e, n1, n0), n0));
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1 else 2",
      MakeIf(e, MakeIf(e, n1, n2), n0));
  CHECK_WRITE("if 0 || 1 then (if 0 || 1 then 1) else 2",
      MakeIf(e, MakeIf(e, n1, n0), n2));
  CHECK_WRITE("if 0 || 1 then 0 else if 0 || 1 then 1 else 2",
      MakeIf(e, n0, MakeIf(e, n1, n2)));
  CHECK_WRITE("if !(0 || 1) then x1 + 1",
      MakeIf(builder.MakeNot(e), MakeBinary(OPPLUS, MakeVariable(0), n1), n0));
}

TEST_F(ExprTest, PiecewiseLinearExprPrecedence) {
  double breakpoints[] = {5, 10};
  double slopes[] = {-1, 0, 1};
  NumericExpr e = builder.MakePiecewiseLinear(
        2, breakpoints, slopes, MakeVariable(42));
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43 ^ 2",
      MakeBinary(OPPOW, e, MakeConst(2)));
}

TEST_F(ExprTest, CallExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1);
  auto f = builder.AddFunction("foo", TestFunc, -1);
  Expr args[] = {
      builder.MakeCall(f, MakeArrayRef<Expr>(0, 0)),
      MakeBinary(OPPLUS, x1, MakeConst(5)),
      MakeConst(7),
      MakeUnary(FLOOR, x2)
  };
  CHECK_WRITE("foo(foo(), x1 + 5, 7, floor(x2))", builder.MakeCall(f, args));
}

TEST_F(ExprTest, VarArgExprPrecedence) {
  NumericExpr x1 = MakeVariable(0), x2 = MakeVariable(1);
  NumericExpr e = MakeBinary(OPPLUS, x1, x2);
  NumericExpr args[] = {e, e};
  CHECK_WRITE("min(x1 + x2, x1 + x2)", builder.MakeVarArg(MINLIST, args));
  CHECK_WRITE("max(x1 + x2, x1 + x2)", builder.MakeVarArg(MAXLIST, args));
  NumericExpr args2[] = {x1}, args3[] = {x2};
  NumericExpr args4[] = {
    builder.MakeVarArg(MINLIST, args2), builder.MakeVarArg(MINLIST, args3)
  };
  CHECK_WRITE("min(min(x1), min(x2))", builder.MakeVarArg(MINLIST, args4));
  NumericExpr args5[] = {
    builder.MakeVarArg(MAXLIST, args2), builder.MakeVarArg(MAXLIST, args3)
  };
  CHECK_WRITE("max(max(x1), max(x2))", builder.MakeVarArg(MAXLIST, args5));
}

TEST_F(ExprTest, SumExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  NumericExpr args1[] = {x2, x3};
  NumericExpr args2[] = {x1, builder.MakeSum(args1)};
  CHECK_WRITE("/* sum */ (x1 + /* sum */ (x2 + x3))", builder.MakeSum(args2));
  NumericExpr args3[] = {x1, MakeBinary(OPMULT, x2, x3)};
  CHECK_WRITE("/* sum */ (x1 + x2 * x3)", builder.MakeSum(args3));
  NumericExpr args4[] = {MakeBinary(OPPLUS, x1, x2), x3};
  CHECK_WRITE("/* sum */ ((x1 + x2) + x3)", builder.MakeSum(args4));
}

TEST_F(ExprTest, CountExprPrecedence) {
  LogicalExpr e = MakeBinaryLogical(OPOR, l0, l1);
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
  NumericExpr e2 = MakeBinary(OPPLUS, x1, x2);
  NumericExpr args3[] = {e2, e2, e2};
  CHECK_WRITE("numberof x1 + x2 in (x1 + x2, x1 + x2)",
              builder.MakeNumberOf(args3));
}

TEST_F(ExprTest, NotExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), x = MakeVariable(0);
  CHECK_WRITE("if !!(x1 = 0) then 1",
      MakeIf(builder.MakeNot(
               builder.MakeNot(MakeRelational(EQ, x, n0))), n1, n0));
}

TEST_F(ExprTest, LogicalOrExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 || 1 || 1 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeBinaryLogical(OPOR, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 || (1 || 1) then 1",
      MakeIf(MakeBinaryLogical(OPOR, l0, MakeBinaryLogical(OPOR, l1, l1)),
          n1, n0));
  CHECK_WRITE("if 0 || 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(OPOR, l0, MakeBinaryLogical(OPAND, l1, l1)),
          n1, n0));
}

TEST_F(ExprTest, LogicalAndExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 && 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(OPAND, MakeBinaryLogical(OPAND, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 && (1 && 1) then 1",
      MakeIf(MakeBinaryLogical(OPAND, l0, MakeBinaryLogical(OPAND, l1, l1)),
          n1, n0));
  CHECK_WRITE("if 0 <= 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(OPAND, MakeRelational(LE, n0, n1), l1), n1, n0));
}

TEST_F(ExprTest, IffExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 <==> 1 <==> 1 then 1",
      MakeIf(MakeBinaryLogical(OP_IFF, MakeBinaryLogical(OP_IFF, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 <==> (1 <==> 1) then 1",
      MakeIf(MakeBinaryLogical(OP_IFF, l0, MakeBinaryLogical(OP_IFF, l1, l1)),
          n1, n0));
  CHECK_WRITE("if (0 <==> 1) && 1 then 1",
      MakeIf(MakeBinaryLogical(OPAND, MakeBinaryLogical(OP_IFF, l0, l1), l1),
          n1, n0));
}

TEST_F(ExprTest, RelationalExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  auto e1 = MakeBinary(OPPLUS, MakeVariable(0), n1);
  auto e2 = MakeBinary(OPPLUS, MakeVariable(1), n1);
  CHECK_WRITE("if x1 + 1 < x2 + 1 then 1",
      MakeIf(MakeRelational(LT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 <= x2 + 1 then 1",
      MakeIf(MakeRelational(LE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 = x2 + 1 then 1",
      MakeIf(MakeRelational(EQ, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 >= x2 + 1 then 1",
      MakeIf(MakeRelational(GE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 > x2 + 1 then 1",
      MakeIf(MakeRelational(GT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 != x2 + 1 then 1",
      MakeIf(MakeRelational(NE, e1, e2), n1, n0));
}

TEST_F(ExprTest, LogicalCountExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), lhs = MakeConst(42);
  LogicalExpr args[] = {MakeRelational(EQ, MakeVariable(0), MakeConst(0)), l1};
  LogicalExpr count1 =
      builder.MakeLogicalCount(OPATLEAST, lhs, builder.MakeCount(args));
  LogicalExpr args2[] = {count1, l1};
  CountExpr count2 = builder.MakeCount(args2);
  CHECK_WRITE("if atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(OPATLEAST, lhs, count2), n1, n0));
  CHECK_WRITE("if atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(OPATMOST, lhs, count2), n1, n0));
  CHECK_WRITE("if exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(OPEXACTLY, lhs, count2), n1, n0));
  CHECK_WRITE("if !atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(OPNOTATLEAST, lhs, count2), n1, n0));
  CHECK_WRITE("if !atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(OPNOTATMOST, lhs, count2), n1, n0));
  CHECK_WRITE("if !exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(OPNOTEXACTLY, lhs, count2), n1, n0));

  args2[0] = l0;
  CountExpr count = builder.MakeCount(args2);
  CHECK_WRITE("if atleast 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeLogicalCount(OPATLEAST, lhs, count),
          l0), n1, n0));
  CHECK_WRITE("if atmost 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeLogicalCount(OPATMOST, lhs, count),
          l0), n1, n0));
  CHECK_WRITE("if exactly 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeLogicalCount(OPEXACTLY, lhs, count),
          l0), n1, n0));
  CHECK_WRITE("if !atleast 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeLogicalCount(OPNOTATLEAST, lhs, count),
          l0), n1, n0));
  CHECK_WRITE("if !atmost 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeLogicalCount(OPNOTATMOST, lhs, count),
          l0), n1, n0));
  CHECK_WRITE("if !exactly 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(OPOR, MakeLogicalCount(OPNOTEXACTLY, lhs, count),
          l0), n1, n0));
}

TEST_F(ExprTest, IteratedLogicalExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  LogicalExpr args[] = {MakeBinaryLogical(OPAND, l0, l0), l0};
  CHECK_WRITE("if /* forall */ ((0 && 0) && 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ANDLIST, args), n1, n0));
  args[0] = MakeBinaryLogical(OPOR, l0, l0);
  CHECK_WRITE("if /* exists */ ((0 || 0) || 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ORLIST, args), n1, n0));
  args[0] = l0;
  LogicalExpr args2[] = {builder.MakeIteratedLogical(ANDLIST, args), l0};
  CHECK_WRITE("if /* forall */ (/* forall */ (0 && 0) && 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ANDLIST, args2), n1, n0));
  args2[0] = builder.MakeIteratedLogical(ORLIST, args);
  CHECK_WRITE("if /* exists */ (/* exists */ (0 || 0) || 0) then 1",
      MakeIf(builder.MakeIteratedLogical(ORLIST, args2), n1, n0));
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
      MakeIf(MakeImplication(l0, MakeBinaryLogical(OPOR, l1, l0), l1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 <==> 0) else 1 then 1", MakeIf(
           MakeImplication(l0, MakeBinaryLogical(OP_IFF, l1, l0), l1), n1, n0));
}

TEST_F(ExprTest, AllDiffExprPrecedence) {
  NumericConstant n0 = MakeConst(0), n1 = MakeConst(1);
  NumericExpr args[] = {MakeBinary(OPPLUS, n0, n1), MakeBinary(OPPLUS, n0, n1)};
  CHECK_WRITE("if alldiff(0 + 1, 0 + 1) then 1",
      MakeIf(builder.MakeAllDiff(args), n1, n0));
}

#ifdef HAVE_UNORDERED_MAP

using ampl::internal::HashCombine;

TEST_F(ExprTest, HashNumericConstant) {
  size_t hash = HashCombine(0, OPNUM);
  hash = HashCombine(hash, 42.0);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(builder.MakeNumericConstant(42)));
}

TEST_F(ExprTest, HashVariable) {
  size_t hash = HashCombine(0, OPVARVAL);
  hash = HashCombine(hash, 42);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(MakeVariable(42)));
}

template <Expr::Kind K, typename Base>
void CheckHash(ampl::BasicUnaryExpr<K, Base> e) {
  size_t hash = HashCombine(0, e.opcode());
  hash = HashCombine<Base>(hash, e.arg());
  EXPECT_EQ(hash, std::hash<Base>()(e));
}

TEST_F(ExprTest, HashUnaryExpr) {
  CheckHash(builder.MakeUnary(OPUMINUS, builder.MakeVariable(0)));
}

template <typename Expr, typename Base, typename Arg>
void CheckHashBinary(Expr e) {
  size_t hash = HashCombine(0, e.opcode());
  hash = HashCombine<Arg>(hash, e.lhs());
  hash = HashCombine<Arg>(hash, e.rhs());
  EXPECT_EQ(hash, std::hash<Base>()(e));
}

template <Expr::Kind K, typename Base, typename Arg>
void CheckHash(ampl::BasicBinaryExpr<K, Base, Arg> e) {
  CheckHashBinary<ampl::BasicBinaryExpr<K, Base, Arg>, Base, Arg>(e);
}

TEST_F(ExprTest, HashBinaryExpr) {
  CheckHash(builder.MakeBinary(OPPLUS, builder.MakeVariable(9), n2));
}

template <typename Base>
void CheckHash(ampl::BasicIfExpr<Base> e) {
  size_t hash = HashCombine(0, e.opcode());
  hash = HashCombine<LogicalExpr>(hash, e.condition());
  hash = HashCombine<Base>(hash, e.true_expr());
  hash = HashCombine<Base>(hash, e.false_expr());
  EXPECT_EQ(hash, std::hash<Base>()(e));
}

TEST_F(ExprTest, HashIfExpr) {
  CheckHash(builder.MakeIf(l0, builder.MakeVariable(2), n1));
}

TEST_F(ExprTest, HashPiecewiseLinearExpr) {
  size_t hash = HashCombine(0, OPPLTERM);
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

size_t HashString(const char *s) {
  size_t hash = HashCombine(0, OPHOL);
  for (; *s; ++s)
    hash = HashCombine(hash, *s);
  return hash;
}

TEST_F(ExprTest, HashCallExpr) {
  enum {NUM_ARGS = 3};
  Variable var = MakeVariable(9);
  Expr args[NUM_ARGS] = {n1, builder.MakeStringLiteral("test"), var};
  Function f = builder.AddFunction(
        "foo", TestFunc, NUM_ARGS, Function::SYMBOLIC);
  size_t hash = HashCombine(0, OPFUNCALL);
  hash = HashCombine(hash, f.name());
  hash = HashCombine<NumericExpr>(hash, n1);
  hash = HashCombine(hash, HashString("test"));
  hash = HashCombine<NumericExpr>(hash, var);
  EXPECT_EQ(hash, std::hash<NumericExpr>()(builder.MakeCall(f, args)));
}

template <typename Expr, typename Arg, typename Base>
size_t CheckHash(Expr e) {
  size_t hash = HashCombine(0, e.opcode());
  for (typename Expr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    hash = HashCombine<Arg>(hash, *i);
  EXPECT_EQ(hash, std::hash<Base>()(e));
  return hash;
}

template <Expr::Kind K>
void CheckHash(ampl::BasicIteratedExpr<K> e) {
  CheckHash<ampl::BasicIteratedExpr<K>,
      typename ampl::ExprInfo<K>::Arg, typename ampl::ExprInfo<K>::Base>(e);
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
        builder.MakeVarArg(MINLIST, args));
  CheckHash(builder.MakeSum(args));
  CheckHash(builder.MakeNumberOf(args));
}

TEST_F(ExprTest, HashCountExpr) {
  LogicalExpr args[] = {l0, l1, l1};
  CheckHash(builder.MakeCount(args));
}

TEST_F(ExprTest, HashLogicalConstant) {
  size_t hash = HashCombine(0, OPNUM);
  hash = HashCombine(hash, true);
  EXPECT_EQ(hash, std::hash<LogicalExpr>()(l1));
}

TEST_F(ExprTest, HashNotExpr) {
  CheckHash(builder.MakeNot(l1));
}

TEST_F(ExprTest, HashBinaryLogicalExpr) {
  CheckHash(builder.MakeBinaryLogical(OPOR, l1, l0));
}

TEST_F(ExprTest, HashRelationalExpr) {
  CheckHash(builder.MakeRelational(LT, builder.MakeVariable(6), n2));
}

TEST_F(ExprTest, HashLogicalCountExpr) {
  LogicalExpr args[] = {l0, l1};
  CheckHashBinary<LogicalCountExpr, LogicalExpr, NumericExpr>(
        builder.MakeLogicalCount(OPATMOST, n1, builder.MakeCount(args)));
}

TEST_F(ExprTest, HashImplicationExpr) {
  CheckHash(builder.MakeImplication(l1, l0, builder.MakeLogicalConstant(true)));
}

TEST_F(ExprTest, HashIteratedLogical) {
  LogicalExpr args[] = {l0, l1, builder.MakeLogicalConstant(false)};
  CheckHash(builder.MakeIteratedLogical(ORLIST, args));
}

TEST_F(ExprTest, HashAllDiff) {
  NumericExpr args[] = {n1, n2, builder.MakeVariable(4)};
  CheckHash(builder.MakeAllDiff(args));
}

TEST_F(ExprTest, HashStringLiteral) {
  // String literal can only occur as a function argument, so test
  // it as a part of a call expression.
  Expr args[] = {builder.MakeStringLiteral("test")};
  Function f = builder.AddFunction("foo", TestFunc, 1, Function::SYMBOLIC);
  size_t hash = HashCombine(0, OPFUNCALL);
  hash = HashCombine(hash, f.name());
  hash = HashCombine(hash, HashString("test"));
  EXPECT_EQ(hash, std::hash<NumericExpr>()(builder.MakeCall(f, args)));
}

TEST_F(ExprTest, HashNumberOfArgs) {
  size_t hash = HashCombine<NumericExpr>(0, MakeVariable(11));
  hash = HashCombine<NumericExpr>(hash, MakeConst(22));
  NumericExpr args[] = {MakeConst(42), MakeVariable(11), MakeConst(22)};
  EXPECT_EQ(hash, ampl::internal::HashNumberOfArgs()(
              builder.MakeNumberOf(args)));
}

TEST_F(ExprTest, EqualNumberOfArgs) {
  using ampl::internal::EqualNumberOfArgs;
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
  using ampl::internal::MatchNumberOfArgs;
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
}

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
#include "tests/expr-builder.h"

#include "solvers/util/problem.h"

using ampl::CallArg;
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
using ampl::ExprVisitor;

using ampl::LinearTerm;
using ampl::LinearExpr;

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

class ExprTest : public ::testing::Test, public ampl::ExprBuilder {};

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
  EXPECT_EQ(Expr::LOGICAL_END, Expr::EXPR_END);
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
    Expr e(MakeExpr(&raw));
    EXPECT_EQ(OPMINUS, e.opcode());
  }
  {
    expr raw = {reinterpret_cast<efunc*>(OPOR)};
    Expr e(MakeExpr(&raw));
    EXPECT_EQ(OPOR, e.opcode());
  }
}

TEST_F(ExprTest, SafeBool) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = {reinterpret_cast<efunc*>(42)};
  Expr e2(MakeExpr(&raw2));
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
  {OPHOL, "string", Expr::UNKNOWN},
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
      Expr e(MakeExpr(&raw));
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
  EXPECT_EQ(65, CheckExpr<Expr>(Expr::EXPR_START, Expr::EXPR_END, -1));
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
  MakeExpr<TestExpr>(&raw1);
  expr raw2 = {reinterpret_cast<efunc*>(OPUMINUS)};  // unary
  MakeExpr<TestExpr>(&raw2);
  TestAssertInCreate<TestExpr>(OPPLTERM);  // neither
}

TEST_F(ExprTest, EqualityOperator) {
  expr raw1 = {}, raw2 = {};
  EXPECT_TRUE(MakeExpr(&raw1) != Expr());
  EXPECT_TRUE(MakeExpr(&raw1) != MakeExpr(&raw2));
  EXPECT_TRUE(Expr() == Expr());
  EXPECT_TRUE(MakeExpr(&raw1) == MakeExpr(&raw1));
  EXPECT_TRUE(MakeExpr(&raw2) == MakeExpr(&raw2));
}

TEST_F(ExprTest, EqualNum) {
  EXPECT_TRUE(Equal(AddNum(0.42), AddNum(0.42)));
  EXPECT_FALSE(Equal(AddNum(0.42), AddNum(42)));
}

TEST_F(ExprTest, EqualVar) {
  EXPECT_TRUE(Equal(AddVar(0), AddVar(0)));
  EXPECT_FALSE(Equal(AddVar(0), AddVar(1)));
  EXPECT_FALSE(Equal(AddVar(0), AddNum(0)));
}

TEST_F(ExprTest, EqualUnary) {
  EXPECT_TRUE(Equal(AddUnary(OPUMINUS, AddVar(0)),
                    AddUnary(OPUMINUS, AddVar(0))));
  EXPECT_FALSE(Equal(AddUnary(OPUMINUS, AddVar(0)),
                     AddVar(0)));
  EXPECT_FALSE(Equal(AddUnary(OPUMINUS, AddVar(0)),
                     AddUnary(FLOOR, AddVar(0))));
  EXPECT_FALSE(Equal(AddUnary(OPUMINUS, AddVar(0)),
                     AddUnary(OPUMINUS, AddVar(1))));
}

TEST_F(ExprTest, EqualBinary) {
  EXPECT_TRUE(Equal(AddBinary(OPPLUS, AddVar(0), AddNum(42)),
                    AddBinary(OPPLUS, AddVar(0), AddNum(42))));
  EXPECT_FALSE(Equal(AddBinary(OPPLUS, AddVar(0), AddNum(42)),
                     AddBinary(OPMINUS, AddVar(0), AddNum(42))));
  EXPECT_FALSE(Equal(AddBinary(OPPLUS, AddVar(0), AddNum(42)),
                     AddBinary(OPPLUS, AddNum(42), AddVar(0))));
  EXPECT_FALSE(Equal(AddBinary(OPPLUS, AddVar(0), AddNum(42)),
                     AddBinary(OPPLUS, AddVar(0), AddNum(0))));
  EXPECT_FALSE(Equal(AddNum(42),
                     AddBinary(OPPLUS, AddVar(0), AddNum(42))));
}

TEST_F(ExprTest, EqualVarArg) {
  EXPECT_TRUE(Equal(
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42)),
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42)),
      AddVarArg(MINLIST, AddVar(0), AddVar(1))));
  EXPECT_FALSE(Equal(
      AddVarArg(MINLIST, AddVar(0), AddVar(1)),
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42)),
      AddVarArg(MAXLIST, AddVar(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42)),
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(0))));
  EXPECT_FALSE(Equal(
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_TRUE(Equal(AddPL(5, args, 0), AddPL(5, args, 0)));
  EXPECT_FALSE(Equal(AddPL(5, args, 0), AddPL(3, args, 0)));
  EXPECT_FALSE(Equal(AddPL(5, args, 0), AddPL(5, args, 1)));
  double args2[] = {-1, 5, 0, 11, 1};
  EXPECT_FALSE(Equal(AddPL(5, args, 0), AddPL(5, args2, 0)));
  EXPECT_FALSE(Equal(AddPL(5, args, 0), AddNum(42)));
}

TEST_F(ExprTest, EqualIf) {
  EXPECT_TRUE(Equal(
      AddIf(AddBool(0), AddVar(1), AddNum(42)),
      AddIf(AddBool(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddIf(AddBool(0), AddVar(1), AddNum(42)),
      AddSum(AddVar(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddIf(AddBool(0), AddVar(1), AddNum(42)),
      AddIf(AddBool(0), AddVar(1), AddNum(0))));
  EXPECT_FALSE(Equal(
      AddIf(AddBool(0), AddVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualSum) {
  EXPECT_TRUE(Equal(
      AddSum(AddVar(0), AddVar(1), AddNum(42)),
      AddSum(AddVar(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddSum(AddVar(0), AddVar(1), AddNum(42)),
      AddSum(AddVar(0), AddVar(1))));
  EXPECT_FALSE(Equal(
      AddSum(AddVar(0), AddVar(1)),
      AddSum(AddVar(0), AddVar(1), AddNum(42))));
  EXPECT_FALSE(Equal(
      AddSum(AddVar(0), AddVar(1), AddNum(42)),
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
  EXPECT_FALSE(Equal(
      AddSum(AddVar(0), AddVar(1), AddNum(42)),
      AddSum(AddVar(0), AddVar(1), AddNum(0))));
  EXPECT_FALSE(Equal(
      AddSum(AddVar(0), AddVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualCount) {
  EXPECT_TRUE(Equal(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
  EXPECT_FALSE(Equal(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddCount(AddBool(false), AddBool(true))));
  EXPECT_FALSE(Equal(
      AddCount(AddBool(false), AddBool(true)),
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
  EXPECT_FALSE(Equal(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddSum(AddNum(0), AddNum(1), AddNum(1))));
  EXPECT_FALSE(Equal(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddCount(AddBool(false), AddBool(true), AddBool(false))));
  EXPECT_FALSE(Equal(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddBool(true)));
}

TEST_F(ExprTest, NumericExpr) {
  EXPECT_EQ(45,
      CheckExpr<NumericExpr>(Expr::NUMERIC_START, Expr::NUMERIC_END, OPNOT));
}

TEST_F(ExprTest, LogicalExpr) {
  EXPECT_EQ(21,
      CheckExpr<LogicalExpr>(Expr::LOGICAL_START, Expr::LOGICAL_END));
}

TEST_F(ExprTest, UnaryExpr) {
  EXPECT_EQ(21, CheckExpr<UnaryExpr>(Expr::UNARY));
  NumericExpr arg(AddNum(42));
  UnaryExpr e(AddUnary(OPUMINUS, arg));
  EXPECT_EQ(arg, e.arg());
}

TEST_F(ExprTest, BinaryExpr) {
  EXPECT_EQ(14, CheckExpr<BinaryExpr>(Expr::BINARY));
  NumericExpr lhs(AddNum(42)), rhs(AddNum(43));
  BinaryExpr e(AddBinary(OPDIV, lhs, rhs));
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
}

TEST_F(ExprTest, VarArgExpr) {
  EXPECT_EQ(2, CheckExpr<VarArgExpr>(Expr::VARARG));
  NumericExpr args[] = {AddNum(42), AddNum(43), AddNum(44)};
  VarArgExpr e(AddVarArg(MINLIST, args[0], args[1], args[2]));
  int index = 0;
  VarArgExpr::iterator i = e.begin();
  for (VarArgExpr::iterator end = e.end(); i != end; ++i, ++index) {
    EXPECT_TRUE(*i);
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index].opcode(), i->opcode());
  }
  EXPECT_FALSE(*i);
  EXPECT_EQ(3, index);
  i = e.begin();
  VarArgExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

TEST_F(ExprTest, SumExpr) {
  EXPECT_EQ(1, CheckExpr<SumExpr>(Expr::SUM));
  EXPECT_EQ(0, AddSum().num_args());
  NumericExpr args[] = {AddNum(42), AddNum(43), AddNum(44)};
  SumExpr e(AddSum(args[0], args[1], args[2]));
  EXPECT_EQ(3, e.num_args());
  int index = 0;
  SumExpr::iterator i = e.begin();
  for (SumExpr::iterator end = e.end(); i != end; ++i, ++index) {
    EXPECT_TRUE(*i);
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index].opcode(), i->opcode());
  }
  EXPECT_EQ(3, index);
  i = e.begin();
  SumExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

TEST_F(ExprTest, CountExpr) {
  EXPECT_EQ(1, CheckExpr<CountExpr>(Expr::COUNT));
  LogicalExpr args[] = {AddBool(false), AddBool(true), AddBool(false)};
  EXPECT_EQ(2, AddCount(args[0], args[1]).num_args());
  CountExpr e(AddCount(args[0], args[1], args[2]));
  EXPECT_EQ(3, e.num_args());
  int index = 0;
  CountExpr::iterator i = e.begin();
  for (CountExpr::iterator end = e.end(); i != end; ++i, ++index) {
    EXPECT_TRUE(*i);
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index].opcode(), i->opcode());
  }
  EXPECT_EQ(3, index);
  i = e.begin();
  CountExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

TEST_F(ExprTest, IfExpr) {
  EXPECT_EQ(1, CheckExpr<IfExpr>(Expr::IF));
  LogicalExpr condition(AddBool(true));
  NumericExpr true_expr(AddNum(42)), false_expr(AddNum(43));
  IfExpr e(AddIf(condition, true_expr, false_expr));
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
}

TEST_F(ExprTest, PiecewiseLinearExpr) {
  EXPECT_EQ(1,
      CheckExpr<PiecewiseLinearExpr>(Expr::PLTERM, Expr::PLTERM, OPPLUS));
  double args[] = {-1, 5, 0, 10, 1};
  PiecewiseLinearExpr e(AddPL(5, args, 42));
  EXPECT_EQ(2, e.num_breakpoints());
  EXPECT_EQ(5, e.breakpoint(0));
  EXPECT_EQ(10, e.breakpoint(1));
  EXPECT_DEBUG_DEATH(e.breakpoint(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(e.breakpoint(2);, "Assertion");  // NOLINT(*)
  EXPECT_EQ(3, e.num_slopes());
  EXPECT_EQ(-1, e.slope(0));
  EXPECT_EQ(0, e.slope(1));
  EXPECT_EQ(1, e.slope(2));
  EXPECT_DEBUG_DEATH(e.slope(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(e.slope(3);, "Assertion");  // NOLINT(*)
  EXPECT_EQ(42, e.var_index());
}

TEST_F(ExprTest, NumericConstant) {
  EXPECT_EQ(1, CheckExpr<NumericConstant>(Expr::CONSTANT));
  NumericConstant e(AddNum(42));
  EXPECT_EQ(42, e.value());
}

TEST_F(ExprTest, Variable) {
  EXPECT_EQ(1, CheckExpr<Variable>(Expr::VARIABLE));
  Variable e(AddVar(42));
  EXPECT_EQ(42, e.index());
}

TEST_F(ExprTest, NumberOfExpr) {
  EXPECT_EQ(1, CheckExpr<NumberOfExpr>(Expr::NUMBEROF));
  NumericExpr value = AddNum(42);
  NumericExpr args[] = {AddNum(43), AddNum(44)};
  EXPECT_EQ(1, AddNumberOf(value, args[0]).num_args());
  NumberOfExpr e(AddNumberOf(value, args[0], args[1]));
  EXPECT_EQ(2, e.num_args());
  EXPECT_EQ(value, e.value());
  int index = 0;
  NumberOfExpr::iterator i = e.begin();
  for (NumberOfExpr::iterator end = e.end(); i != end; ++i, ++index) {
    EXPECT_TRUE(*i);
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index].opcode(), i->opcode());
    EXPECT_EQ(args[index], e[index]);
  }
  EXPECT_EQ(2, index);
  i = e.begin();
  NumberOfExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

TEST_F(ExprTest, CallExpr) {
  EXPECT_EQ(1, CheckExpr<CallExpr>(Expr::CALL));
  CallArg args[] = {
      CallArg(3, AddNum(42)), CallArg(5), CallArg(7, AddNum(44))
  };
  CallExpr e(AddCall("foo", args, args + 3));
  EXPECT_STREQ("foo", e.function().name());
  EXPECT_EQ(3, e.function().num_args());
  EXPECT_EQ(3, e.num_args());
  CallExpr e2(AddCall("bar", args, args + 2));
  EXPECT_EQ(2, e2.function().num_args());
  EXPECT_EQ(2, e2.num_args());
  EXPECT_EQ(e.function(), e.function());
  EXPECT_NE(e.function(), e2.function());
  EXPECT_FALSE(ampl::Function());
  EXPECT_EQ(2, e.num_arg_exprs());
  int index = 0, arg_expr_count = 0;
  CallExpr::arg_expr_iterator i = e.arg_expr_begin();
  for (; index != sizeof(args) / sizeof(*args); ++index) {
    NumericExpr expr = args[index].expr();
    EXPECT_EQ(args[index].constant(), e.arg_constant(index));
    if (!expr) continue;
    EXPECT_EQ(expr, *i);
    EXPECT_EQ(expr.opcode(), i->opcode());
    EXPECT_EQ(index, e.arg_index(i));
    ++arg_expr_count;
    ++i;
  }
  EXPECT_EQ(e.arg_expr_end(), i);
  EXPECT_EQ(3, index);
  EXPECT_EQ(2, arg_expr_count);
  i = e.arg_expr_begin();
  CallExpr::arg_expr_iterator i2 = i++;
  EXPECT_EQ(args[0].expr(), *i2);
  EXPECT_EQ(args[2].expr(), *i);
}

TEST_F(ExprTest, CallExprArgs) {
  NumericExpr n42 = AddNum(42), n44 = AddNum(44);
  CallExpr e(AddCall("foo", CallArg(3, n42), 5, CallArg(7, n44)));
  CallExpr::Args args(e);
  EXPECT_EQ(3, args[0].constant());
  EXPECT_EQ(n42, args[0].expr());
  EXPECT_EQ(5, args[1].constant());
  EXPECT_TRUE(!args[1].expr());
  EXPECT_EQ(7, args[2].constant());
  EXPECT_EQ(n44, args[2].expr());
}

TEST_F(ExprTest, LogicalConstant) {
  EXPECT_EQ(1, CheckExpr<LogicalConstant>(Expr::CONSTANT));
  LogicalConstant e(AddBool(true));
  EXPECT_TRUE(e.value());
  e = AddBool(false);
  EXPECT_FALSE(e.value());
}

TEST_F(ExprTest, RelationalExpr) {
  EXPECT_EQ(6, CheckExpr<RelationalExpr>(Expr::RELATIONAL));
  NumericExpr lhs(AddNum(42)), rhs(AddNum(43));
  RelationalExpr e(AddRelational(EQ, lhs, rhs));
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
}

TEST_F(ExprTest, NotExpr) {
  EXPECT_EQ(1, CheckExpr<NotExpr>(Expr::NOT));
  LogicalExpr arg(AddBool(true));
  NotExpr e(AddNot(arg));
  EXPECT_EQ(arg, e.arg());
}

TEST_F(ExprTest, LogicalCountExpr) {
  EXPECT_EQ(6, CheckExpr<LogicalCountExpr>(Expr::LOGICAL_COUNT));
  NumericExpr value(AddNum(42));
  NumericExpr n(AddNum(42));
  CountExpr count(AddCount(AddRelational(EQ, n, n), AddRelational(EQ, n, n)));
  LogicalCountExpr e(AddLogicalCount(OPATLEAST, value, count));
  EXPECT_EQ(value, e.value());
  EXPECT_EQ(count, e.count());
}

TEST_F(ExprTest, BinaryLogicalExpr) {
  EXPECT_EQ(3, CheckExpr<BinaryLogicalExpr>(Expr::BINARY_LOGICAL));
  LogicalExpr lhs(AddBool(false)), rhs(AddBool(true));
  BinaryLogicalExpr e(AddBinaryLogical(OPOR, lhs, rhs));
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
}

TEST_F(ExprTest, ImplicationExpr) {
  EXPECT_EQ(1, CheckExpr<ImplicationExpr>(Expr::IMPLICATION));
  LogicalExpr condition(AddBool(true));
  LogicalExpr true_expr(AddBool(true)), false_expr(AddBool(false));
  ImplicationExpr e(AddImplication(condition, true_expr, false_expr));
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
}

TEST_F(ExprTest, IteratedLogicalExpr) {
  EXPECT_EQ(2, CheckExpr<IteratedLogicalExpr>(Expr::ITERATED_LOGICAL));
  LogicalExpr args[] = {AddBool(false), AddBool(true), AddBool(false)};
  IteratedLogicalExpr e(AddIteratedLogical(ORLIST, args[0], args[1], args[2]));
  EXPECT_EQ(3, e.num_args());
  int index = 0;
  IteratedLogicalExpr::iterator i = e.begin();
  for (IteratedLogicalExpr::iterator end = e.end(); i != end; ++i, ++index) {
    EXPECT_TRUE(*i);
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index].opcode(), i->opcode());
    EXPECT_EQ(args[index], e[index]);
  }
  EXPECT_EQ(3, index);
  i = e.begin();
  IteratedLogicalExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

TEST_F(ExprTest, AllDiffExpr) {
  EXPECT_EQ(1, CheckExpr<AllDiffExpr>(Expr::ALLDIFF));
  NumericExpr args[] = {AddNum(42), AddNum(43), AddNum(44)};
  AllDiffExpr e(AddAllDiff(args[0], args[1], args[2]));
  EXPECT_EQ(3, e.num_args());
  int index = 0;
  AllDiffExpr::iterator i = e.begin();
  for (AllDiffExpr::iterator end = e.end(); i != end; ++i, ++index) {
    EXPECT_TRUE(*i);
    EXPECT_EQ(args[index], *i);
    EXPECT_EQ(args[index].opcode(), i->opcode());
    EXPECT_EQ(args[index], e[index]);
  }
  EXPECT_EQ(3, index);
  i = e.begin();
  AllDiffExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

struct TestResult {
  NumericExpr expr;
};

struct TestLResult {
  LogicalExpr expr;
};

// Use different classes for Result and LResult to check that it is possible.
struct FullTestVisitor : ExprVisitor<FullTestVisitor, TestResult, TestLResult> {
  TestResult Handle(NumericExpr e) {
    TestResult result = {e};
    return result;
  }

  TestLResult Handle(LogicalExpr e) {
    TestLResult result = {e};
    return result;
  }

  TestResult VisitPlus(BinaryExpr e) { return Handle(e); }
  TestResult VisitMinus(BinaryExpr e) { return Handle(e); }
  TestResult VisitMult(BinaryExpr e) { return Handle(e); }
  TestResult VisitDiv(BinaryExpr e) { return Handle(e); }
  TestResult VisitRem(BinaryExpr e) { return Handle(e); }
  TestResult VisitPow(BinaryExpr e) { return Handle(e); }
  TestResult VisitNumericLess(BinaryExpr e) { return Handle(e); }
  TestResult VisitMin(VarArgExpr e) { return Handle(e); }
  TestResult VisitMax(VarArgExpr e) { return Handle(e); }
  TestResult VisitFloor(UnaryExpr e) { return Handle(e); }
  TestResult VisitCeil(UnaryExpr e) { return Handle(e); }
  TestResult VisitAbs(UnaryExpr e) { return Handle(e); }
  TestResult VisitUnaryMinus(UnaryExpr e) { return Handle(e); }
  TestResult VisitIf(IfExpr e) { return Handle(e); }
  TestResult VisitTanh(UnaryExpr e) { return Handle(e); }
  TestResult VisitTan(UnaryExpr e) { return Handle(e); }
  TestResult VisitSqrt(UnaryExpr e) { return Handle(e); }
  TestResult VisitSinh(UnaryExpr e) { return Handle(e); }
  TestResult VisitSin(UnaryExpr e) { return Handle(e); }
  TestResult VisitLog10(UnaryExpr e) { return Handle(e); }
  TestResult VisitLog(UnaryExpr e) { return Handle(e); }
  TestResult VisitExp(UnaryExpr e) { return Handle(e); }
  TestResult VisitCosh(UnaryExpr e) { return Handle(e); }
  TestResult VisitCos(UnaryExpr e) { return Handle(e); }
  TestResult VisitAtanh(UnaryExpr e) { return Handle(e); }
  TestResult VisitAtan2(BinaryExpr e) { return Handle(e); }
  TestResult VisitAtan(UnaryExpr e) { return Handle(e); }
  TestResult VisitAsinh(UnaryExpr e) { return Handle(e); }
  TestResult VisitAsin(UnaryExpr e) { return Handle(e); }
  TestResult VisitAcosh(UnaryExpr e) { return Handle(e); }
  TestResult VisitAcos(UnaryExpr e) { return Handle(e); }
  TestResult VisitSum(SumExpr e) { return Handle(e); }
  TestResult VisitIntDiv(BinaryExpr e) { return Handle(e); }
  TestResult VisitPrecision(BinaryExpr e) { return Handle(e); }
  TestResult VisitRound(BinaryExpr e) { return Handle(e); }
  TestResult VisitTrunc(BinaryExpr e) { return Handle(e); }
  TestResult VisitCount(CountExpr e) { return Handle(e); }
  TestResult VisitNumberOf(NumberOfExpr e) { return Handle(e); }
  TestResult VisitCall(CallExpr e) { return Handle(e); }
  TestResult VisitPiecewiseLinear(PiecewiseLinearExpr e) { return Handle(e); }
  TestResult VisitPowConstExp(BinaryExpr e) { return Handle(e); }
  TestResult VisitPow2(UnaryExpr e) { return Handle(e); }
  TestResult VisitPowConstBase(BinaryExpr e) { return Handle(e); }
  TestResult VisitNumericConstant(NumericConstant e) { return Handle(e); }
  TestResult VisitVariable(Variable e) { return Handle(e); }
  TestLResult VisitOr(BinaryLogicalExpr e) { return Handle(e); }
  TestLResult VisitAnd(BinaryLogicalExpr e) { return Handle(e); }
  TestLResult VisitLess(RelationalExpr e) { return Handle(e); }
  TestLResult VisitLessEqual(RelationalExpr e) { return Handle(e); }
  TestLResult VisitEqual(RelationalExpr e) { return Handle(e); }
  TestLResult VisitGreaterEqual(RelationalExpr e) { return Handle(e); }
  TestLResult VisitGreater(RelationalExpr e) { return Handle(e); }
  TestLResult VisitNotEqual(RelationalExpr e) { return Handle(e); }
  TestLResult VisitNot(NotExpr e) { return Handle(e); }
  TestLResult VisitAtLeast(LogicalCountExpr e) { return Handle(e); }
  TestLResult VisitAtMost(LogicalCountExpr e) { return Handle(e); }
  TestLResult VisitExactly(LogicalCountExpr e) { return Handle(e); }
  TestLResult VisitNotAtLeast(LogicalCountExpr e) { return Handle(e); }
  TestLResult VisitNotAtMost(LogicalCountExpr e) { return Handle(e); }
  TestLResult VisitNotExactly(LogicalCountExpr e) { return Handle(e); }
  TestLResult VisitForAll(IteratedLogicalExpr e) { return Handle(e); }
  TestLResult VisitExists(IteratedLogicalExpr e) { return Handle(e); }
  TestLResult VisitImplication(ImplicationExpr e) { return Handle(e); }
  TestLResult VisitIff(BinaryLogicalExpr e) { return Handle(e); }
  TestLResult VisitAllDiff(AllDiffExpr e) { return Handle(e); }
  TestLResult VisitLogicalConstant(LogicalConstant e) { return Handle(e); }
};

TEST_F(ExprTest, ExprVisitorHandlesAll) {
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  for (int i = 0; i < size; ++i) {
    FullTestVisitor v;
    const OpInfo &info = OP_INFO[i];
    if (info.kind == Expr::UNKNOWN) continue;
    expr raw = {reinterpret_cast<efunc*>(info.code)};
    Expr e(MakeExpr(&raw));
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
    if (info.kind == Expr::UNKNOWN) continue;
    expr raw = {reinterpret_cast<efunc*>(info.code)};
    Expr e(MakeExpr(&raw));
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
  EXPECT_THROW(NullVisitor().Visit(AddNum(0)), UnsupportedExprError);
  EXPECT_THROW(NullVisitor().Visit(AddBool(false)), UnsupportedExprError);
}

TEST_F(ExprTest, ExprVisitorInvalidThrows) {
  expr raw = {reinterpret_cast<efunc*>(OPNUM)};
  NumericExpr ne(MakeExpr<NumericExpr>(&raw));
  LogicalExpr le(MakeExpr<LogicalExpr>(&raw));
  raw.op = reinterpret_cast<efunc*>(-1);
  EXPECT_THROW(NullVisitor().Visit(ne), InvalidNumericExprError);
  EXPECT_THROW(NullVisitor().Visit(le), InvalidLogicalExprError);
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
  NumberOfExpr e1 = AddNumberOf(AddNum(11), AddVar(0));
  NumberOfExpr e2 = AddNumberOf(AddNum(22), AddVar(1));
  map.Add(11, e1);
  map.Add(22, e2);
  map.Add(33, AddNumberOf(AddNum(33), AddVar(0)));
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
  EXPECT_TRUE(IsZero(AddNum(0)));
  EXPECT_FALSE(IsZero(AddNum(1)));
  EXPECT_FALSE(IsZero(AddNum(-1)));
  EXPECT_FALSE(IsZero(AddNum(10)));
  EXPECT_FALSE(IsZero(AddVar(0)));
}

// Checks if WriteExpr produces the expected output for expr.
static ::testing::AssertionResult CheckWrite(
    const char *, const char *expr_str,
    const std::string &expected_output, NumericExpr expr) {
  fmt::Writer w;
  WriteExpr(w, ampl::LinearObjExpr(), expr);
  auto actual_output = w.str();
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
  CHECK_WRITE("0", AddNum(0));
  CHECK_WRITE("42", AddNum(42));
  CHECK_WRITE("12.34", AddNum(12.34));
}

TEST_F(ExprTest, WriteVariable) {
  CHECK_WRITE("x1", AddVar(0));
  CHECK_WRITE("x3", AddVar(2));
}

TEST_F(ExprTest, WriteBinaryExpr) {
  Variable x1 = AddVar(0);
  NumericConstant n42 = AddNum(42);
  CHECK_WRITE("x1 + 42", AddBinary(OPPLUS, x1, n42));
  CHECK_WRITE("x1 - 42", AddBinary(OPMINUS, x1, n42));
  CHECK_WRITE("x1 * 42", AddBinary(OPMULT, x1, n42));
  CHECK_WRITE("x1 / 42", AddBinary(OPDIV, x1, n42));
  CHECK_WRITE("x1 mod 42", AddBinary(OPREM, x1, n42));
  CHECK_WRITE("x1 ^ 42", AddBinary(OPPOW, x1, n42));
  CHECK_WRITE("x1 ^ 42", AddBinary(OPPOW, x1, n42));
  CHECK_WRITE("x1 ^ 42", AddBinary(OP1POW, x1, n42));
  CHECK_WRITE("x1 ^ 42", AddBinary(OPCPOW, x1, n42));
  CHECK_WRITE("x1 less 42", AddBinary(OPLESS, x1, n42));
  CHECK_WRITE("x1 div 42", AddBinary(OPintDIV, x1, n42));
}

TEST_F(ExprTest, WriteBinaryFunc) {
  auto x1 = AddVar(0);
  auto n42 = AddNum(42);
  CHECK_WRITE("atan2(x1, 42)", AddBinary(OP_atan2, x1, n42));
  CHECK_WRITE("precision(x1, 42)", AddBinary(OPprecision, x1, n42));
  CHECK_WRITE("round(x1, 42)", AddBinary(OPround, x1, n42));
  CHECK_WRITE("trunc(x1, 42)", AddBinary(OPtrunc, x1, n42));
}

TEST_F(ExprTest, WriteUnaryExpr) {
  auto x1 = AddVar(0);
  CHECK_WRITE("-x1", AddUnary(OPUMINUS, x1));
  CHECK_WRITE("x1 ^ 2", AddUnary(OP2POW, x1));
  int count = 0;
  for (int i = 0, size = sizeof(OP_INFO) / sizeof(*OP_INFO); i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int code = info.code;
    if (info.kind != Expr::UNARY || code == OPUMINUS || code == OP2POW)
      continue;
    CHECK_WRITE(fmt::format("{}(x1)", info.str), AddUnary(code, x1));
    ++count;
  }
  EXPECT_EQ(19, count);
}

TEST_F(ExprTest, WriteVarArgExpr) {
  CHECK_WRITE("min(x1, x2, 42)",
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42)));
  CHECK_WRITE("max(x1, x2, 42)",
      AddVarArg(MAXLIST, AddVar(0), AddVar(1), AddNum(42)));
}

TEST_F(ExprTest, WriteIfExpr) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  CHECK_WRITE("if x1 = 0 then 1",
      AddIf(AddRelational(EQ, AddVar(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 0 else 1",
      AddIf(AddRelational(EQ, AddVar(0), n0), n0, n1));
}

TEST_F(ExprTest, WriteSumExpr) {
  CHECK_WRITE("/* sum */ (x1 + x2 + 42)",
    AddSum(AddVar(0), AddVar(1), AddNum(42)));
  CHECK_WRITE("/* sum */ ((x1 + x2) + 42)",
    AddSum(AddBinary(OPPLUS, AddVar(0), AddVar(1)), AddNum(42)));
}

TEST_F(ExprTest, WriteCountExpr) {
  CHECK_WRITE("count(x1 = 0, 1, 0)",
      AddCount(AddRelational(EQ, AddVar(0), AddNum(0)),
          AddBool(true), AddBool(false)));
}

TEST_F(ExprTest, WriteNumberOfExpr) {
  CHECK_WRITE("numberof 42 in (43, 44)",
      AddNumberOf(AddNum(42), AddNum(43), AddNum(44)));
}

TEST_F(ExprTest, WritePiecewiseLinearExpr) {
  double args[] = {-1, 5, 0, 10, 1};
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43", AddPL(5, args, 42));
}

TEST_F(ExprTest, WriteCallExpr) {
  CHECK_WRITE("foo(3, x1 + 5, 7, x2)",
      AddCall("foo", 3, CallArg(5, AddVar(0)), 7, CallArg(0, AddVar(1))));
  CHECK_WRITE("foo(x1 + x2)",
      AddCall("foo", AddBinary(OPPLUS, AddVar(0), AddVar(1))));
}

TEST_F(ExprTest, WriteNotExpr) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  CHECK_WRITE("if !(x1 = 0) then 1",
      AddIf(AddNot(AddRelational(EQ, AddVar(0), n0)), n1, n0));
}

TEST_F(ExprTest, WriteBinaryLogicalExpr) {
  auto e1 = AddRelational(GT, AddVar(0), AddNum(0));
  auto e2 = AddRelational(LT, AddVar(0), AddNum(10));
  CHECK_WRITE("if x1 > 0 || x1 < 10 then 1",
      AddIf(AddBinaryLogical(OPOR, e1, e2), AddNum(1), AddNum(0)));
  CHECK_WRITE("if x1 > 0 && x1 < 10 then 1",
      AddIf(AddBinaryLogical(OPAND, e1, e2), AddNum(1), AddNum(0)));
  CHECK_WRITE("if x1 > 0 <==> x1 < 10 then 1",
      AddIf(AddBinaryLogical(OP_IFF, e1, e2), AddNum(1), AddNum(0)));
}

TEST_F(ExprTest, WriteRelationalExpr) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  CHECK_WRITE("if x1 < 0 then 1",
      AddIf(AddRelational(LT, AddVar(0), n0), n1, n0));
  CHECK_WRITE("if x1 <= 0 then 1",
      AddIf(AddRelational(LE, AddVar(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 1",
      AddIf(AddRelational(EQ, AddVar(0), n0), n1, n0));
  CHECK_WRITE("if x1 >= 0 then 1",
      AddIf(AddRelational(GE, AddVar(0), n0), n1, n0));
  CHECK_WRITE("if x1 > 0 then 1",
      AddIf(AddRelational(GT, AddVar(0), n0), n1, n0));
  CHECK_WRITE("if x1 != 0 then 1",
      AddIf(AddRelational(NE, AddVar(0), n0), n1, n0));
}

TEST_F(ExprTest, WriteLogicalCountExpr) {
  auto n0 = AddNum(0), n1 = AddNum(1), value = AddNum(42);
  auto count = AddCount(
      AddRelational(EQ, AddVar(0), AddNum(0)), AddBool(true), AddBool(false));
  CHECK_WRITE("if atleast 42 (x1 = 0, 1, 0) then 1",
      AddIf(AddLogicalCount(OPATLEAST, value, count), n1, n0));
  CHECK_WRITE("if atmost 42 (x1 = 0, 1, 0) then 1",
      AddIf(AddLogicalCount(OPATMOST, value, count), n1, n0));
  CHECK_WRITE("if exactly 42 (x1 = 0, 1, 0) then 1",
      AddIf(AddLogicalCount(OPEXACTLY, value, count), n1, n0));
  CHECK_WRITE("if !atleast 42 (x1 = 0, 1, 0) then 1",
      AddIf(AddLogicalCount(OPNOTATLEAST, value, count), n1, n0));
  CHECK_WRITE("if !atmost 42 (x1 = 0, 1, 0) then 1",
      AddIf(AddLogicalCount(OPNOTATMOST, value, count), n1, n0));
  CHECK_WRITE("if !exactly 42 (x1 = 0, 1, 0) then 1",
      AddIf(AddLogicalCount(OPNOTEXACTLY, value, count), n1, n0));
}

TEST_F(ExprTest, WriteIteratedLogicalExpr) {
  auto e1 = AddRelational(EQ, AddVar(0), AddNum(0));
  auto e2 = AddBool(true), e3 = AddBool(false);
  CHECK_WRITE("if /* forall */ (x1 = 0 && 1 && 0) then 1",
      AddIf(AddIteratedLogical(ANDLIST, e1, e2, e3), AddNum(1), AddNum(0)));
  CHECK_WRITE("if /* exists */ (x1 = 0 || 1 || 0) then 1",
      AddIf(AddIteratedLogical(ORLIST, e1, e2, e3), AddNum(1), AddNum(0)));
}

TEST_F(ExprTest, WriteImplicationExpr) {
  auto e1 = AddRelational(EQ, AddVar(0), AddNum(0));
  auto e2 = AddBool(true), e3 = AddBool(false);
  CHECK_WRITE("if x1 = 0 ==> 1 then 1",
      AddIf(AddImplication(e1, e2, e3), AddNum(1), AddNum(0)));
  CHECK_WRITE("if x1 = 0 ==> 0 else 1 then 1",
      AddIf(AddImplication(e1, e3, e2), AddNum(1), AddNum(0)));
}

TEST_F(ExprTest, WriteAllDiffExpr) {
  CHECK_WRITE("if alldiff(42, 43, 44) then 1",
      AddIf(AddAllDiff(AddNum(42), AddNum(43), AddNum(44)),
          AddNum(1), AddNum(0)));
}

TEST_F(ExprTest, UnaryExprPrecedence) {
  auto x1 = AddVar(0);
  CHECK_WRITE("--x1", AddUnary(OPUMINUS, AddUnary(OPUMINUS, x1)));
  CHECK_WRITE("-(x1 ^ x1)", AddUnary(OPUMINUS, AddBinary(OPPOW, x1, x1)));
}

TEST_F(ExprTest, UnaryFuncPrecedence) {
  auto x1 = AddVar(0);
  int count = 0;
  for (int i = 0, size = sizeof(OP_INFO) / sizeof(*OP_INFO); i < size; ++i) {
    const OpInfo &info = OP_INFO[i];
    int code = info.code;
    if (info.kind != Expr::UNARY || code == OPUMINUS || code == OP2POW)
      continue;
    CHECK_WRITE(fmt::format("{0}({0}(x1))", info.str),
        AddUnary(code, AddUnary(code, x1)));
    CHECK_WRITE(fmt::format("{0}(x1 + x1)", info.str),
        AddUnary(code, AddBinary(OPPLUS, x1, x1)));
    ++count;
  }
  EXPECT_EQ(19, count);
}

TEST_F(ExprTest, Pow2Precedence) {
  auto x1 = AddVar(0);
  CHECK_WRITE("(x1 ^ 2) ^ 2", AddUnary(OP2POW, AddUnary(OP2POW, x1)));
  CHECK_WRITE("(x1 * x1) ^ 2", AddUnary(OP2POW, AddBinary(OPMULT, x1, x1)));
}

TEST_F(ExprTest, AdditiveExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1), x3 = AddVar(2);
  CHECK_WRITE("x1 + x2 + x3",
      AddBinary(OPPLUS, AddBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("x1 + x2 - x3",
      AddBinary(OPMINUS, AddBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("x1 + x2 less x3",
      AddBinary(OPLESS, AddBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("x1 + (x2 + x3)",
      AddBinary(OPPLUS, x1, AddBinary(OPPLUS, x2, x3)));
  CHECK_WRITE("(x1 + x2) * x3",
      AddBinary(OPMULT, AddBinary(OPPLUS, x1, x2), x3));
  CHECK_WRITE("if 1 then x1 else x2 + x3",
      AddIf(AddBool(true), x1, AddBinary(OPPLUS, x2, x3)));
}

TEST_F(ExprTest, MultiplicativeExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1), x3 = AddVar(2);
  CHECK_WRITE("x1 * x2 * x3",
      AddBinary(OPMULT, AddBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * x2 / x3",
      AddBinary(OPDIV, AddBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * x2 div x3",
      AddBinary(OPintDIV, AddBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * x2 mod x3",
      AddBinary(OPREM, AddBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("x1 * (x2 * x3)",
      AddBinary(OPMULT, x1, AddBinary(OPMULT, x2, x3)));
  CHECK_WRITE("(x1 * x2) ^ x3",
      AddBinary(OPPOW, AddBinary(OPMULT, x1, x2), x3));
  CHECK_WRITE("(x1 + x2) * x3",
      AddBinary(OPMULT, AddBinary(OPPLUS, x1, x2), x3));
}

TEST_F(ExprTest, ExponentiationExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1), x3 = AddVar(2);
  CHECK_WRITE("x1 ^ x2 ^ x3",
      AddBinary(OPPOW, x1, AddBinary(OPPOW, x2, x3)));
  CHECK_WRITE("x1 ^ x2 ^ 3",
      AddBinary(OPPOW, x1, AddBinary(OP1POW, x2, AddNum(3))));
  CHECK_WRITE("x1 ^ 3 ^ x2",
      AddBinary(OPPOW, x1, AddBinary(OPCPOW, AddNum(3), x2)));
  CHECK_WRITE("(x1 ^ 2) ^ 3",
      AddBinary(OP1POW, AddBinary(OP1POW, x1, AddNum(2)), AddNum(3)));
  CHECK_WRITE("-x1 ^ -x2",
      AddBinary(OPPOW, AddUnary(OPUMINUS, x1), AddUnary(OPUMINUS, x2)));
  CHECK_WRITE("x1 ^ (x2 * x3)",
      AddBinary(OPPOW, x1, AddBinary(OPMULT, x2, x3)));
}

TEST_F(ExprTest, BinaryFuncPrecedence) {
  auto x1 = AddVar(0);
  auto e = AddBinary(OPPLUS, x1, x1);
  CHECK_WRITE("atan2(atan2(x1, x1), x1 + x1)",
      AddBinary(OP_atan2, AddBinary(OP_atan2, x1, x1), e));
  CHECK_WRITE("precision(precision(x1, x1), x1 + x1)",
      AddBinary(OPprecision, AddBinary(OPprecision, x1, x1), e));
  CHECK_WRITE("round(round(x1, x1), x1 + x1)",
      AddBinary(OPround, AddBinary(OPround, x1, x1), e));
  CHECK_WRITE("trunc(trunc(x1, x1), x1 + x1)",
      AddBinary(OPtrunc, AddBinary(OPtrunc, x1, x1), e));
}

TEST_F(ExprTest, VarArgExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1);
  auto e = AddBinary(OPPLUS, x1, x2);
  CHECK_WRITE("min(x1 + x2, x1 + x2)", AddVarArg(MINLIST, e, e));
  CHECK_WRITE("max(x1 + x2, x1 + x2)", AddVarArg(MAXLIST, e, e));
  CHECK_WRITE("min(min(x1), min(x2))",
      AddVarArg(MINLIST, AddVarArg(MINLIST, x1), AddVarArg(MINLIST, x2)));
  CHECK_WRITE("max(max(x1), max(x2))",
      AddVarArg(MAXLIST, AddVarArg(MAXLIST, x1), AddVarArg(MAXLIST, x2)));
}

TEST_F(ExprTest, IfExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1), n2 = AddNum(2);
  auto e = AddBinaryLogical(OPOR, AddBool(false), AddBool(true));
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1",
      AddIf(e, AddIf(e, n1, n0), n0));
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1 else 2",
      AddIf(e, AddIf(e, n1, n2), n0));
  CHECK_WRITE("if 0 || 1 then (if 0 || 1 then 1) else 2",
      AddIf(e, AddIf(e, n1, n0), n2));
  CHECK_WRITE("if 0 || 1 then 0 else if 0 || 1 then 1 else 2",
      AddIf(e, n0, AddIf(e, n1, n2)));
  CHECK_WRITE("if !(0 || 1) then x1 + 1",
      AddIf(AddNot(e), AddBinary(OPPLUS, AddVar(0), n1), n0));
}

TEST_F(ExprTest, SumExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1), x3 = AddVar(2);
  CHECK_WRITE("/* sum */ (x1 + /* sum */ (x2 + x3))",
      AddSum(x1, AddSum(x2, x3)));
  CHECK_WRITE("/* sum */ (x1 + x2 * x3)",
      AddSum(x1, AddBinary(OPMULT, x2, x3)));
  CHECK_WRITE("/* sum */ ((x1 + x2) + x3)",
      AddSum(AddBinary(OPPLUS, x1, x2), x3));
}

TEST_F(ExprTest, CountExprPrecedence) {
  auto e = AddBinaryLogical(OPOR, AddBool(false), AddBool(false));
  CHECK_WRITE("count(0 || 0, 0 || 0)", AddCount(e, e));
}

TEST_F(ExprTest, NumberOfExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1);
  auto e = AddNumberOf(AddNum(42), x1, x2);
  CHECK_WRITE("numberof numberof 42 in (x1, x2) in ("
      "numberof 42 in (x1, x2), numberof 42 in (x1, x2))",
      AddNumberOf(e, e, e));
  auto e2 = AddBinary(OPPLUS, x1, x2);
  CHECK_WRITE("numberof x1 + x2 in (x1 + x2, x1 + x2)",
      AddNumberOf(e2, e2, e2));
}

TEST_F(ExprTest, PiecewiseLinearExprPrecedence) {
  double args[] = {-1, 5, 0, 10, 1};
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43 ^ 2",
      AddBinary(OPPOW, AddPL(5, args, 42), AddNum(2)));
}

TEST_F(ExprTest, CallExprPrecedence) {
  auto x1 = AddVar(0), x2 = AddVar(1);
  CHECK_WRITE("foo(foo(), x1 + 5, 7, floor(x2))", AddCall("foo",
      AddCall("foo", 0, 0), CallArg(5, x1), 7, AddUnary(FLOOR, x2)));
}

TEST_F(ExprTest, NotExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  CHECK_WRITE("if !!(x1 = 0) then 1",
      AddIf(AddNot(AddNot(AddRelational(EQ, AddVar(0), n0))), n1, n0));
}

TEST_F(ExprTest, LogicalOrExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  auto b0 = AddBool(false), b1 = AddBool(true);
  CHECK_WRITE("if 0 || 1 || 1 then 1",
      AddIf(AddBinaryLogical(OPOR, AddBinaryLogical(OPOR, b0, b1), b1),
          n1, n0));
  CHECK_WRITE("if 0 || (1 || 1) then 1",
      AddIf(AddBinaryLogical(OPOR, b0, AddBinaryLogical(OPOR, b1, b1)),
          n1, n0));
  CHECK_WRITE("if 0 || 1 && 1 then 1",
      AddIf(AddBinaryLogical(OPOR, b0, AddBinaryLogical(OPAND, b1, b1)),
          n1, n0));
}

TEST_F(ExprTest, LogicalAndExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  auto b0 = AddBool(false), b1 = AddBool(true);
  CHECK_WRITE("if 0 && 1 && 1 then 1",
      AddIf(AddBinaryLogical(OPAND, AddBinaryLogical(OPAND, b0, b1), b1),
          n1, n0));
  CHECK_WRITE("if 0 && (1 && 1) then 1",
      AddIf(AddBinaryLogical(OPAND, b0, AddBinaryLogical(OPAND, b1, b1)),
          n1, n0));
  CHECK_WRITE("if 0 <= 1 && 1 then 1",
      AddIf(AddBinaryLogical(OPAND, AddRelational(LE, n0, n1), b1), n1, n0));
}

TEST_F(ExprTest, IffExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  auto b0 = AddBool(false), b1 = AddBool(true);
  CHECK_WRITE("if 0 <==> 1 <==> 1 then 1",
      AddIf(AddBinaryLogical(OP_IFF, AddBinaryLogical(OP_IFF, b0, b1), b1),
          n1, n0));
  CHECK_WRITE("if 0 <==> (1 <==> 1) then 1",
      AddIf(AddBinaryLogical(OP_IFF, b0, AddBinaryLogical(OP_IFF, b1, b1)),
          n1, n0));
  CHECK_WRITE("if (0 <==> 1) && 1 then 1",
      AddIf(AddBinaryLogical(OPAND, AddBinaryLogical(OP_IFF, b0, b1), b1),
          n1, n0));
}

TEST_F(ExprTest, RelationalExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  auto e1 = AddBinary(OPPLUS, AddVar(0), n1);
  auto e2 = AddBinary(OPPLUS, AddVar(1), n1);
  CHECK_WRITE("if x1 + 1 < x2 + 1 then 1",
      AddIf(AddRelational(LT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 <= x2 + 1 then 1",
      AddIf(AddRelational(LE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 = x2 + 1 then 1",
      AddIf(AddRelational(EQ, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 >= x2 + 1 then 1",
      AddIf(AddRelational(GE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 > x2 + 1 then 1",
      AddIf(AddRelational(GT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 != x2 + 1 then 1",
      AddIf(AddRelational(NE, e1, e2), n1, n0));
}

TEST_F(ExprTest, LogicalCountExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1), value = AddNum(42);
  auto count1 = AddLogicalCount(OPATLEAST, value, AddCount(
      AddRelational(EQ, AddVar(0), AddNum(0)), AddBool(true)));
  auto count2 = AddCount(count1, AddBool(true));
  CHECK_WRITE("if atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      AddIf(AddLogicalCount(OPATLEAST, value, count2), n1, n0));
  CHECK_WRITE("if atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      AddIf(AddLogicalCount(OPATMOST, value, count2), n1, n0));
  CHECK_WRITE("if exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      AddIf(AddLogicalCount(OPEXACTLY, value, count2), n1, n0));
  CHECK_WRITE("if !atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      AddIf(AddLogicalCount(OPNOTATLEAST, value, count2), n1, n0));
  CHECK_WRITE("if !atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      AddIf(AddLogicalCount(OPNOTATMOST, value, count2), n1, n0));
  CHECK_WRITE("if !exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      AddIf(AddLogicalCount(OPNOTEXACTLY, value, count2), n1, n0));

  auto count = AddCount(AddBool(false), AddBool(true));
  CHECK_WRITE("if atleast 42 (0, 1) || 0 then 1",
      AddIf(AddBinaryLogical(OPOR, AddLogicalCount(OPATLEAST, value, count),
          AddBool(false)), n1, n0));
  CHECK_WRITE("if atmost 42 (0, 1) || 0 then 1",
      AddIf(AddBinaryLogical(OPOR, AddLogicalCount(OPATMOST, value, count),
          AddBool(false)), n1, n0));
  CHECK_WRITE("if exactly 42 (0, 1) || 0 then 1",
      AddIf(AddBinaryLogical(OPOR, AddLogicalCount(OPEXACTLY, value, count),
          AddBool(false)), n1, n0));
  CHECK_WRITE("if !atleast 42 (0, 1) || 0 then 1",
      AddIf(AddBinaryLogical(OPOR, AddLogicalCount(OPNOTATLEAST, value, count),
          AddBool(false)), n1, n0));
  CHECK_WRITE("if !atmost 42 (0, 1) || 0 then 1",
      AddIf(AddBinaryLogical(OPOR, AddLogicalCount(OPNOTATMOST, value, count),
          AddBool(false)), n1, n0));
  CHECK_WRITE("if !exactly 42 (0, 1) || 0 then 1",
      AddIf(AddBinaryLogical(OPOR, AddLogicalCount(OPNOTEXACTLY, value, count),
          AddBool(false)), n1, n0));
}

TEST_F(ExprTest, IteratedLogicalExprPrecedence) {
  auto b0 = AddBool(false);
  auto n0 = AddNum(0), n1 = AddNum(1);
  CHECK_WRITE("if /* forall */ ((0 && 0) && 0) then 1",
      AddIf(AddIteratedLogical(ANDLIST, AddBinaryLogical(OPAND, b0, b0), b0),
          n1, n0));
  CHECK_WRITE("if /* exists */ ((0 || 0) || 0) then 1",
      AddIf(AddIteratedLogical(ORLIST, AddBinaryLogical(OPOR, b0, b0), b0),
          n1, n0));
  CHECK_WRITE("if /* forall */ (/* forall */ (0 && 0) && 0) then 1",
      AddIf(AddIteratedLogical(ANDLIST,
          AddIteratedLogical(ANDLIST, b0, b0), b0), n1, n0));
  CHECK_WRITE("if /* exists */ (/* exists */ (0 || 0) || 0) then 1",
      AddIf(AddIteratedLogical(ORLIST,
          AddIteratedLogical(ORLIST, b0, b0), b0), n1, n0));
}

TEST_F(ExprTest, ImplicationExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  auto b0 = AddBool(false), b1 = AddBool(true);
  CHECK_WRITE("if 0 ==> 1 ==> 0 then 1",
      AddIf(AddImplication(AddImplication(b0, b1, b0), b0, b0), n1, n0));
  CHECK_WRITE("if 0 ==> 1 ==> 0 else 1 then 1",
      AddIf(AddImplication(AddImplication(b0, b1, b0), b0, b1), n1, n0));
  CHECK_WRITE("if 0 ==> 1 else 0 ==> 1 then 1",
      AddIf(AddImplication(b0, b1, AddImplication(b0, b1, b0)), n1, n0));
  CHECK_WRITE("if 0 ==> (1 ==> 0) else 1 then 1",
      AddIf(AddImplication(b0, AddImplication(b1, b0, b0), b1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 ==> 0 else 1) then 1",
      AddIf(AddImplication(b0, AddImplication(b1, b0, b1), b0), n1, n0));
  CHECK_WRITE("if 0 ==> 1 || 0 else 1 then 1",
      AddIf(AddImplication(b0, AddBinaryLogical(OPOR, b1, b0), b1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 <==> 0) else 1 then 1",
      AddIf(AddImplication(b0, AddBinaryLogical(OP_IFF, b1, b0), b1), n1, n0));
}

TEST_F(ExprTest, AllDiffExprPrecedence) {
  auto n0 = AddNum(0), n1 = AddNum(1);
  CHECK_WRITE("if alldiff(0 + 1, 0 + 1) then 1",
      AddIf(AddAllDiff(AddBinary(OPPLUS, n0, n1),
          AddBinary(OPPLUS, n0, n1)), n1, n0));
}

#ifdef HAVE_UNORDERED_MAP

using ampl::HashCombine;

TEST_F(ExprTest, HashNum) {
  size_t hash = 0;
  HashCombine(hash, OPNUM);
  HashCombine(hash, 42.0);
  EXPECT_EQ(hash, std::hash<Expr>()(AddNum(42)));
}

TEST_F(ExprTest, HashVar) {
  size_t hash = 0;
  HashCombine(hash, OPVARVAL);
  HashCombine(hash, 42);
  EXPECT_EQ(hash, std::hash<Expr>()(AddVar(42)));
}

TEST_F(ExprTest, HashUnary) {
  size_t hash = 0;
  HashCombine(hash, OPUMINUS);
  HashCombine<Expr>(hash, AddVar(0));
  EXPECT_EQ(hash, std::hash<Expr>()(AddUnary(OPUMINUS, AddVar(0))));
}

TEST_F(ExprTest, HashBinary) {
  size_t hash = 0;
  HashCombine(hash, OPPLUS);
  HashCombine<Expr>(hash, AddVar(11));
  HashCombine<Expr>(hash, AddNum(22));
  EXPECT_EQ(hash, std::hash<Expr>()(AddBinary(OPPLUS, AddVar(11), AddNum(22))));
}

TEST_F(ExprTest, HashVarArg) {
  size_t hash = 0;
  HashCombine(hash, MINLIST);
  HashCombine<Expr>(hash, AddVar(0));
  HashCombine<Expr>(hash, AddVar(1));
  HashCombine<Expr>(hash, AddNum(42));
  EXPECT_EQ(hash, std::hash<Expr>()(
      AddVarArg(MINLIST, AddVar(0), AddVar(1), AddNum(42))));
}

TEST_F(ExprTest, HashPiecewiseLinear) {
  size_t hash = 0;
  HashCombine(hash, OPPLTERM);
  double args[] = {-1, 5, 0, 10, 1};
  for (size_t i = 0; i < sizeof(args) / sizeof(*args); ++i)
    HashCombine(hash, args[i]);
  HashCombine<Expr>(hash, AddVar(11));
  EXPECT_EQ(hash, std::hash<Expr>()(AddPL(5, args, 11)));
}

TEST_F(ExprTest, HashIf) {
  size_t hash = 0;
  HashCombine(hash, OPIFnl);
  HashCombine<Expr>(hash, AddBool(0));
  HashCombine<Expr>(hash, AddVar(1));
  HashCombine<Expr>(hash, AddNum(42));
  EXPECT_EQ(hash, std::hash<Expr>()(AddIf(AddBool(0), AddVar(1), AddNum(42))));
}

TEST_F(ExprTest, HashSum) {
  size_t hash = 0;
  HashCombine(hash, OPSUMLIST);
  HashCombine<Expr>(hash, AddVar(0));
  HashCombine<Expr>(hash, AddVar(1));
  HashCombine<Expr>(hash, AddNum(42));
  EXPECT_EQ(hash, std::hash<Expr>()(AddSum(AddVar(0), AddVar(1), AddNum(42))));
}

TEST_F(ExprTest, HashCount) {
  size_t hash = 0;
  HashCombine(hash, OPCOUNT);
  HashCombine<Expr>(hash, AddBool(false));
  HashCombine<Expr>(hash, AddBool(true));
  HashCombine<Expr>(hash, AddBool(true));
  EXPECT_EQ(hash, std::hash<Expr>()(
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
}

TEST_F(ExprTest, HashNumberOfArgs) {
  size_t hash = 0;
  HashCombine<Expr>(hash, AddVar(11));
  HashCombine<Expr>(hash, AddNum(22));
  EXPECT_EQ(hash, ampl::HashNumberOfArgs()(
      AddNumberOf(AddNum(42), AddVar(11), AddNum(22))));
}

TEST_F(ExprTest, EqualNumberOfArgs) {
  EXPECT_TRUE(ampl::EqualNumberOfArgs()(
      AddNumberOf(AddNum(0), AddVar(11), AddNum(22)),
      AddNumberOf(AddNum(1), AddVar(11), AddNum(22))));
  EXPECT_FALSE(ampl::EqualNumberOfArgs()(
      AddNumberOf(AddNum(0), AddVar(11), AddNum(22)),
      AddNumberOf(AddNum(1), AddVar(11))));
  EXPECT_FALSE(ampl::EqualNumberOfArgs()(
      AddNumberOf(AddNum(0), AddVar(11), AddNum(22)),
      AddNumberOf(AddNum(1), AddVar(11), AddNum(33))));
}

struct TestNumberOf {
  NumberOfExpr expr;

  TestNumberOf(NumberOfExpr e) : expr(e) {}
};

TEST_F(ExprTest, MatchNumberOfArgs) {
  EXPECT_TRUE(ampl::MatchNumberOfArgs<TestNumberOf>(
      AddNumberOf(AddNum(1), AddVar(11), AddNum(22)))(
          TestNumberOf(AddNumberOf(AddNum(0), AddVar(11), AddNum(22)))));
  EXPECT_FALSE(ampl::MatchNumberOfArgs<TestNumberOf>(
      AddNumberOf(AddNum(1), AddVar(11)))(
          TestNumberOf(AddNumberOf(AddNum(0), AddVar(11), AddNum(22)))));
  EXPECT_FALSE(ampl::MatchNumberOfArgs<TestNumberOf>(
      AddNumberOf(AddNum(1), AddVar(11), AddNum(33)))(
          TestNumberOf(AddNumberOf(AddNum(0), AddVar(11), AddNum(22)))));
}
#endif
}

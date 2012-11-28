/*
 Tests of the the C++ interface to AMPL expression trees.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
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
#include "tests/expr_builder.h"

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
using ampl::PiecewiseLinearTerm;
using ampl::NumericConstant;
using ampl::Variable;
using ampl::NumberOfExpr;
using ampl::LogicalConstant;
using ampl::RelationalExpr;
using ampl::NotExpr;
using ampl::BinaryLogicalExpr;
using ampl::ImplicationExpr;
using ampl::IteratedLogicalExpr;
using ampl::AllDiffExpr;
using ampl::ExprVisitor;

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
}

namespace ampl {
namespace internal {
template <>
bool Is<TestExpr>(Expr e) {
  return e.kind() == Expr::BINARY || e.kind() == Expr::UNARY;
}
}
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
  const char *name;
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
  {OPIFnl, "if-then-else", Expr::IF},
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
  {OPATLEAST, "atleast", Expr::RELATIONAL},
  {OPATMOST,  "atmost",  Expr::RELATIONAL},
  {OPPLTERM, "pl term", Expr::PLTERM},
  {OPIFSYM,  "string if-then-else", Expr::UNKNOWN},
  {OPEXACTLY,    "exactly",     Expr::RELATIONAL},
  {OPNOTATLEAST, "not atleast", Expr::RELATIONAL},
  {OPNOTATMOST,  "not atmost",  Expr::RELATIONAL},
  {OPNOTEXACTLY, "not exactly", Expr::RELATIONAL},
  {ANDLIST, "forall", Expr::ITERATED_LOGICAL},
  {ORLIST,  "exists", Expr::ITERATED_LOGICAL},
  {OPIMPELSE, "implies else", Expr::IMPLICATION},
  {OP_IFF, "iff", Expr::BINARY_LOGICAL},
  {OPALLDIFF, "alldiff", Expr::ALLDIFF},
  {OP1POW, "1pow", Expr::BINARY},
  {OP2POW, "^2",   Expr::UNARY},
  {OPCPOW, "cpow", Expr::BINARY},
  {OPFUNCALL, "function call", Expr::UNKNOWN},
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
    const char *opname = info.name;
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
    EXPECT_STREQ(opname, e.opname());
    ++expr_count;
  }
  EXPECT_GT(expr_count, 0);
  return expr_count;
}

TEST_F(ExprTest, Expr) {
  EXPECT_EQ(64, CheckExpr<Expr>(Expr::EXPR_START, Expr::EXPR_END, -1));
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
  EXPECT_TRUE(Equal(AddPLTerm(5, args, 0), AddPLTerm(5, args, 0)));
  EXPECT_FALSE(Equal(AddPLTerm(5, args, 0), AddPLTerm(3, args, 0)));
  EXPECT_FALSE(Equal(AddPLTerm(5, args, 0), AddPLTerm(5, args, 1)));
  double args2[] = {-1, 5, 0, 11, 1};
  EXPECT_FALSE(Equal(AddPLTerm(5, args, 0), AddPLTerm(5, args2, 0)));
  EXPECT_FALSE(Equal(AddPLTerm(5, args, 0), AddNum(42)));
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
  EXPECT_EQ(44,
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
  NumericExpr args[] = {AddNum(42), AddNum(43), AddNum(44)};
  SumExpr e(AddSum(args[0], args[1], args[2]));
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

TEST_F(ExprTest, PiecewiseLinearTerm) {
  EXPECT_EQ(1,
      CheckExpr<PiecewiseLinearTerm>(Expr::PLTERM, Expr::PLTERM, OPPLUS));
  double args[] = {-1, 5, 0, 10, 1};
  PiecewiseLinearTerm e(AddPLTerm(5, args, 42));
  EXPECT_EQ(2, e.num_breakpoints());
  EXPECT_EQ(5, e.GetBreakpoint(0));
  EXPECT_EQ(10, e.GetBreakpoint(1));
  EXPECT_DEBUG_DEATH(e.GetBreakpoint(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(e.GetBreakpoint(2);, "Assertion");  // NOLINT(*)
  EXPECT_EQ(3, e.num_slopes());
  EXPECT_EQ(-1, e.GetSlope(0));
  EXPECT_EQ(0, e.GetSlope(1));
  EXPECT_EQ(1, e.GetSlope(2));
  EXPECT_DEBUG_DEATH(e.GetSlope(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(e.GetSlope(3);, "Assertion");  // NOLINT(*)
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
  }
  EXPECT_EQ(2, index);
  i = e.begin();
  NumberOfExpr::iterator i2 = i++;
  EXPECT_EQ(args[0], *i2);
  EXPECT_EQ(args[1], *i);
}

TEST_F(ExprTest, LogicalConstant) {
  EXPECT_EQ(1, CheckExpr<LogicalConstant>(Expr::CONSTANT));
  LogicalConstant e(AddBool(true));
  EXPECT_TRUE(e.value());
  e = AddBool(false);
  EXPECT_FALSE(e.value());
}

TEST_F(ExprTest, RelationalExpr) {
  EXPECT_EQ(12, CheckExpr<RelationalExpr>(Expr::RELATIONAL));
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
  TestResult VisitPLTerm(PiecewiseLinearTerm e) { return Handle(e); }
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
  TestLResult VisitAtLeast(RelationalExpr e) { return Handle(e); }
  TestLResult VisitAtMost(RelationalExpr e) { return Handle(e); }
  TestLResult VisitExactly(RelationalExpr e) { return Handle(e); }
  TestLResult VisitNotAtLeast(RelationalExpr e) { return Handle(e); }
  TestLResult VisitNotAtMost(RelationalExpr e) { return Handle(e); }
  TestLResult VisitNotExactly(RelationalExpr e) { return Handle(e); }
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

TEST_F(ExprTest, HashPLTerm) {
  size_t hash = 0;
  HashCombine(hash, OPPLTERM);
  double args[] = {-1, 5, 0, 10, 1};
  for (size_t i = 0; i < sizeof(args) / sizeof(*args); ++i)
    HashCombine(hash, args[i]);
  HashCombine<Expr>(hash, AddVar(11));
  EXPECT_EQ(hash, std::hash<Expr>()(AddPLTerm(5, args, 11)));
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

TEST_F(ExprTest, DISABLED_NumberOfMapPerformance) {
  // Results in the release mode when using linear search:
  // num_exprs  time,s
  //   10000     0.45
  //   20000     2.04
  //   40000    10.58
  ampl::NumberOfMap<Var, CreateVar> map((CreateVar()));
  NumericConstant n = AddNum(0);
  int num_exprs = 10000;
  std::vector<NumberOfExpr> exprs(num_exprs);
  for (int i = 0; i < num_exprs; ++i)
    exprs[i] = AddNumberOf(n, AddVar(i));
  clock_t start = clock();
  for (int i = 0; i < num_exprs; ++i)
    map.Add(0, exprs[i]);
  clock_t end = clock();
  std::cout << "Executed NumberOfMap.Add " << num_exprs << " times in "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s.\n";
}
}

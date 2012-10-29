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

using ampl::Expr;
using ampl::NumericExpr;
using ampl::LogicalExpr;
using ampl::UnaryExpr;
using ampl::BinaryExpr;
using ampl::VarArgExpr;
using ampl::SumExpr;
using ampl::CountExpr;
using ampl::IfExpr;

using ampl::UnsupportedExprError;

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
template <typename ExprT = Expr>
ExprT MakeExpr(expr *e) { return TestExpr::MakeExpr<ExprT>(e); }

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
    int code = reinterpret_cast<size_t>(ptrs[index]->op);
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
      EXPECT_NE(kind, kinds[j]); // Check if all different.
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
      EXPECT_NE(kind, kinds[j]); // Check if all different.
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
  EXPECT_DEBUG_DEATH(MakeExpr<ExprT>(&raw);, "Assertion");
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
      bool cast_result = ampl::Cast<ExprT>(e);
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
  expr raw1 = {reinterpret_cast<efunc*>(OPPLUS)}; // binary
  MakeExpr<TestExpr>(&raw1);
  expr raw2 = {reinterpret_cast<efunc*>(OPUMINUS)}; // unary
  MakeExpr<TestExpr>(&raw2);
  TestAssertInCreate<TestExpr>(OPPLTERM); // neither
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
  EXPECT_TRUE(AreEqual(AddNum(0.42), AddNum(0.42)));
  EXPECT_FALSE(AreEqual(AddNum(0.42), AddNum(42)));
}

TEST_F(ExprTest, EqualVar) {
  EXPECT_TRUE(AreEqual(NewVar(0), NewVar(0)));
  EXPECT_FALSE(AreEqual(NewVar(0), NewVar(1)));
  EXPECT_FALSE(AreEqual(NewVar(0), AddNum(0)));
}

TEST_F(ExprTest, EqualUnary) {
  EXPECT_TRUE(AreEqual(AddUnary(OPUMINUS, NewVar(0)),
                       AddUnary(OPUMINUS, NewVar(0))));
  EXPECT_FALSE(AreEqual(AddUnary(OPUMINUS, NewVar(0)),
                        NewVar(0)));
  EXPECT_FALSE(AreEqual(AddUnary(OPUMINUS, NewVar(0)),
                        AddUnary(FLOOR, NewVar(0))));
  EXPECT_FALSE(AreEqual(AddUnary(OPUMINUS, NewVar(0)),
                        AddUnary(OPUMINUS, NewVar(1))));
}

TEST_F(ExprTest, EqualBinary) {
  EXPECT_TRUE(AreEqual(AddBinary(OPPLUS, NewVar(0), AddNum(42)),
                       AddBinary(OPPLUS, NewVar(0), AddNum(42))));
  EXPECT_FALSE(AreEqual(AddBinary(OPPLUS, NewVar(0), AddNum(42)),
                        AddBinary(OPMINUS, NewVar(0), AddNum(42))));
  EXPECT_FALSE(AreEqual(AddBinary(OPPLUS, NewVar(0), AddNum(42)),
                        AddBinary(OPPLUS, AddNum(42), NewVar(0))));
  EXPECT_FALSE(AreEqual(AddBinary(OPPLUS, NewVar(0), AddNum(42)),
                        AddBinary(OPPLUS, NewVar(0), AddNum(0))));
  EXPECT_FALSE(AreEqual(AddNum(42),
                        AddBinary(OPPLUS, NewVar(0), AddNum(42))));
}

TEST_F(ExprTest, EqualVarArg) {
  EXPECT_TRUE(AreEqual(
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42)),
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42)),
      AddVarArg(MINLIST, NewVar(0), NewVar(1))));
  EXPECT_FALSE(AreEqual(
      AddVarArg(MINLIST, NewVar(0), NewVar(1)),
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42)),
      AddVarArg(MAXLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42)),
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(0))));
  EXPECT_FALSE(AreEqual(
      AddVarArg(MINLIST, NewVar(0), NewVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_TRUE(AreEqual(
      NewPLTerm(5, args, 0),
      NewPLTerm(5, args, 0)));
  EXPECT_FALSE(AreEqual(
      NewPLTerm(5, args, 0),
      NewPLTerm(3, args, 0)));
  EXPECT_FALSE(AreEqual(
      NewPLTerm(5, args, 0),
      NewPLTerm(5, args, 1)));
  double args2[] = {-1, 5, 0, 11, 1};
  EXPECT_FALSE(AreEqual(
      NewPLTerm(5, args, 0),
      NewPLTerm(5, args2, 0)));
  EXPECT_FALSE(AreEqual(
      NewPLTerm(5, args, 0),
      AddNum(42)));
}

TEST_F(ExprTest, EqualIf) {
  EXPECT_TRUE(AreEqual(
      AddIf(AddLogicalConstant(0), NewVar(1), AddNum(42)),
      AddIf(AddLogicalConstant(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddIf(AddLogicalConstant(0), NewVar(1), AddNum(42)),
      AddSum(NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddIf(AddLogicalConstant(0), NewVar(1), AddNum(42)),
      AddIf(AddLogicalConstant(0), NewVar(1), AddNum(0))));
  EXPECT_FALSE(AreEqual(
      AddIf(AddLogicalConstant(0), NewVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualSum) {
  EXPECT_TRUE(AreEqual(
      AddSum(NewVar(0), NewVar(1), AddNum(42)),
      AddSum(NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddSum(NewVar(0), NewVar(1), AddNum(42)),
      AddSum(NewVar(0), NewVar(1))));
  EXPECT_FALSE(AreEqual(
      AddSum(NewVar(0), NewVar(1)),
      AddSum(NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      AddSum(NewVar(0), NewVar(1), AddNum(42)),
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
  EXPECT_FALSE(AreEqual(
      AddSum(NewVar(0), NewVar(1), AddNum(42)),
      AddSum(NewVar(0), NewVar(1), AddNum(0))));
  EXPECT_FALSE(AreEqual(
      AddSum(NewVar(0), NewVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualCount) {
  EXPECT_TRUE(AreEqual(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
  EXPECT_FALSE(AreEqual(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddCount(AddBool(false), AddBool(true))));
  EXPECT_FALSE(AreEqual(
      AddCount(AddBool(false), AddBool(true)),
      AddCount(AddBool(false), AddBool(true), AddBool(true))));
  EXPECT_FALSE(AreEqual(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddSum(AddNum(0), AddNum(1), AddNum(1))));
  EXPECT_FALSE(AreEqual(
      AddCount(AddBool(false), AddBool(true), AddBool(true)),
      AddCount(AddBool(false), AddBool(true), AddBool(false))));
  EXPECT_FALSE(AreEqual(
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
  LogicalExpr args[] = {
      AddLogicalConstant(0), AddLogicalConstant(1), AddLogicalConstant(0)};
  CountExpr e(AddCount(args[0], args[1], args[2]));
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

// TODO: more tests
}

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

using ampl::UnsupportedExprError;

using ampl::internal::ExprProxy;
using ampl::internal::ExprArrayIterator;

namespace {

class TestExpr : public Expr {
public:
  TestExpr(expr *e) : Expr(e) {}
};

Expr MakeExpr(expr *e) { return TestExpr(e); }

class ExprTest : public ::testing::Test, public ampl::ExprBuilder {};

TEST_F(ExprTest, ExprProxy) {
  expr e = {reinterpret_cast<efunc*>(OPDIV)};
  ExprProxy<NumericExpr> p(&e);
  EXPECT_EQ(OPDIV, p->opcode());
}

TEST_F(ExprTest, ExprArrayIterator) {
  {
    ExprArrayIterator<NumericExpr> i;
    EXPECT_EQ(ExprArrayIterator<NumericExpr>(), i);
  }
  expr exprs[] = {
      {reinterpret_cast<efunc*>(OPDIV)},
      {reinterpret_cast<efunc*>(OPPLUS)},
      {reinterpret_cast<efunc*>(OP_atan)},
  };
  expr *const ptrs[] = {exprs, exprs + 1, exprs + 2};
  ExprArrayIterator<NumericExpr> i(ptrs);
  EXPECT_EQ(ExprArrayIterator<NumericExpr>(ptrs), i);
  EXPECT_NE(ExprArrayIterator<NumericExpr>(), i);
  EXPECT_EQ(OPDIV, (*i).opcode());
  EXPECT_EQ(OPDIV, i->opcode());

  ExprArrayIterator<NumericExpr> i2(++i);
  EXPECT_EQ(i2, i);
  EXPECT_NE(ExprArrayIterator<NumericExpr>(ptrs), i);
  EXPECT_EQ(ExprArrayIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(OPPLUS, i->opcode());

  ExprArrayIterator<NumericExpr> i3(i++);
  EXPECT_NE(i3, i);
  EXPECT_NE(ExprArrayIterator<NumericExpr>(ptrs + 1), i);
  EXPECT_EQ(ExprArrayIterator<NumericExpr>(ptrs + 2), i);
  EXPECT_EQ(OPPLUS, i3->opcode());
  EXPECT_EQ(OP_atan, i->opcode());

  int index = 0;
  for (ExprArrayIterator<NumericExpr>
      i(ptrs), e(ptrs + 3); i != e; ++i, ++index) {
    int code = reinterpret_cast<size_t>(ptrs[index]->op);
    EXPECT_EQ(code, i->opcode());
  }
  EXPECT_EQ(3, index);
  std::vector<NumericExpr> vec;
  std::copy(ExprArrayIterator<NumericExpr>(ptrs),
      ExprArrayIterator<NumericExpr>(ptrs + 3), std::back_inserter(vec));
  EXPECT_EQ(OPPLUS, vec[1].opcode());
}

TEST_F(ExprTest, ExprCtor) {
  {
    Expr e;
    EXPECT_FALSE(e);
  }
  {
    expr raw = {reinterpret_cast<efunc*>(42)};
    Expr e(MakeExpr(&raw));
    EXPECT_EQ(42, e.opcode());
  }
  {
    expr raw = {reinterpret_cast<efunc*>(N_OPS - 1)};
    Expr e(MakeExpr(&raw));
    EXPECT_EQ(N_OPS - 1, e.opcode());
  }
}

TEST_F(ExprTest, ExprOpCodeOutOfRangeInCtor) {
  {
    expr raw = {reinterpret_cast<efunc*>(-1)};
    EXPECT_DEBUG_DEATH(Expr e(MakeExpr(&raw));, "Assertion");
  }
  {
    expr raw = {reinterpret_cast<efunc*>(N_OPS)};
    EXPECT_DEBUG_DEATH(Expr e(MakeExpr(&raw));, "Assertion");
  }
  {
    expr raw = {reinterpret_cast<efunc*>(777)};
    EXPECT_DEBUG_DEATH(Expr e(MakeExpr(&raw));, "Assertion");
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
  {OPNUMBEROF,  "numberof",        Expr::COUNT},
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
  {OPVARVAL, "variable", Expr::VARIABLE}
};

TEST_F(ExprTest, Operators) {
  int known_ops = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  EXPECT_EQ(N_OPS, size);
  for (int i = 0; i < N_OPS; ++i) {
    const OpInfo &info = OP_INFO[i];
    int opcode = info.code;
    const char *opname = info.name;
    expr raw = {reinterpret_cast<efunc*>(opcode)};
    if (info.kind == Expr::UNKNOWN) {
      EXPECT_DEBUG_DEATH(MakeExpr(&raw);, "Assertion");
      continue;
    }
    Expr e(MakeExpr(&raw));
    EXPECT_EQ(opcode, e.opcode());
    EXPECT_STREQ(opname, e.opname());
    ++known_ops;
  }
  EXPECT_EQ(64, known_ops);
}

TEST_F(ExprTest, ExprOpCodeOutOfRangeInAccessors) {
  const char *message =
      "Assertion .*IsOpCodeInRange\\(\\)";
  {
    expr raw = {};
    Expr e(MakeExpr(&raw));
    raw.op = reinterpret_cast<efunc*>(-1);
    EXPECT_DEBUG_DEATH(e.opname();, message);
  }
  {
    expr raw = {};
    Expr e(MakeExpr(&raw));
    raw.op = reinterpret_cast<efunc*>(N_OPS);
    EXPECT_DEBUG_DEATH(e.opname();, message);
  }
  {
    expr raw = {};
    Expr e(MakeExpr(&raw));
    raw.op = reinterpret_cast<efunc*>(777);
    EXPECT_DEBUG_DEATH(e.opname();, message);
  }
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
      NewIf(OPIFnl, NewLogicalConstant(0), NewVar(1), AddNum(42)),
      NewIf(OPIFnl, NewLogicalConstant(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewIf(OPIFnl,  NewLogicalConstant(0), NewVar(1), AddNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewIf(OPIFnl, NewLogicalConstant(0), NewVar(1), AddNum(42)),
      NewIf(OPIFnl, NewLogicalConstant(0), NewVar(1), AddNum(0))));
  EXPECT_FALSE(AreEqual(
      NewIf(OPIFnl, NewLogicalConstant(0), NewVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualSum) {
  EXPECT_TRUE(AreEqual(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(0))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42)),
      AddNum(42)));
}

TEST_F(ExprTest, EqualCount) {
  EXPECT_TRUE(AreEqual(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPCOUNT, NewVar(0), NewVar(1)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), AddNum(42))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(0))));
  EXPECT_FALSE(AreEqual(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), AddNum(42)),
      AddNum(42)));
}

void MakeNumericExpr(int opcode) {
  expr e = {reinterpret_cast<efunc*>(opcode)};
  NumericExpr ne(&e);
}

TEST_F(ExprTest, InvalidNumericExpr) {
  int i = 0, numeric_count = 0;
  for (; i < N_OPS; ++i) {
    const OpInfo &info = OP_INFO[i];
    if (info.kind >= Expr::NUMERIC_START && info.kind <= Expr::NUMERIC_END) {
      MakeNumericExpr(info.code);
      ++numeric_count;
    } else {
      EXPECT_DEBUG_DEATH(MakeNumericExpr(info.code);, "Assertion");
    }
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(N_OPS, i);
  EXPECT_EQ(44, numeric_count);
}

void MakeLogicalExpr(int opcode) {
  expr e = {reinterpret_cast<efunc*>(opcode)};
  LogicalExpr le(&e);
}

TEST_F(ExprTest, InvalidLogicalExpr) {
  int i = 0, logical_count = 0;
  for (; i < N_OPS; ++i) {
    const OpInfo &info = OP_INFO[i];
    if (info.kind >= Expr::LOGICAL_START && info.kind <= Expr::LOGICAL_END) {
      MakeLogicalExpr(info.code);
      ++logical_count;
    } else {
      EXPECT_DEBUG_DEATH(MakeLogicalExpr(info.code);, "Assertion");
    }
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(N_OPS, i);
  EXPECT_EQ(21, logical_count);
}

TEST_F(ExprTest, UnaryExpr) {
  UnaryExpr null_expr;
  EXPECT_FALSE(null_expr);
  NumericExpr arg(AddNum(42));
  UnaryExpr e(AddUnary(OPUMINUS, arg));
  EXPECT_EQ(arg, e.arg());
  EXPECT_DEBUG_DEATH(AddUnary(OPPLUS, arg);,
    "Assertion .*HasKind\\(UNARY\\)");
}

TEST_F(ExprTest, BinaryExpr) {
  BinaryExpr null_expr;
  EXPECT_FALSE(null_expr);
  NumericExpr lhs(AddNum(42)), rhs(AddNum(43));
  BinaryExpr e(AddBinary(OPDIV, lhs, rhs));
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_DEBUG_DEATH(AddBinary(OPUMINUS, lhs, rhs);,
    "Assertion .*HasKind\\(BINARY\\)");
}

TEST_F(ExprTest, VarArgExpr) {
  VarArgExpr null_expr;
  EXPECT_FALSE(null_expr);
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
  EXPECT_DEBUG_DEATH(AddVarArg(OPUMINUS, args[0], args[1], args[2]);,
    "Assertion .*HasKind\\(VARARG\\)");
}

// TODO: more tests
}

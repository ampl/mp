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

using ampl::OPTYPE_UNARY;
using ampl::OPTYPE_BINARY;
using ampl::OPTYPE_VARARG;
using ampl::OPTYPE_PLTERM;
using ampl::OPTYPE_IF;
using ampl::OPTYPE_SUM;
using ampl::OPTYPE_FUNCALL;
using ampl::OPTYPE_STRING;
using ampl::OPTYPE_NUMBER;
using ampl::OPTYPE_VARIABLE;
using ampl::OPTYPE_COUNT;

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
  const char *message =
      "Assertion .*!expr_ \\|\\| IsOpCodeInRange\\(\\)";
  {
    expr raw = {reinterpret_cast<efunc*>(-1)};
    EXPECT_DEBUG_DEATH(Expr e(MakeExpr(&raw));, message);
  }
  {
    expr raw = {reinterpret_cast<efunc*>(N_OPS)};
    EXPECT_DEBUG_DEATH(Expr e(MakeExpr(&raw));, message);
  }
  {
    expr raw = {reinterpret_cast<efunc*>(777)};
    EXPECT_DEBUG_DEATH(Expr e(MakeExpr(&raw));, message);
  }
}

TEST_F(ExprTest, SafeBool) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = {reinterpret_cast<efunc*>(42)};
  Expr e2(MakeExpr(&raw2));
  EXPECT_TRUE(e2);
}

enum {LOGICAL = 1, UNSUPPORTED = 2};

struct OpInfo {
  int code;
  const char *name;
  int optype;
  int type;
};

const OpInfo OP_INFO[] = {
  {OPPLUS,  "+",    OPTYPE_BINARY},
  {OPMINUS, "-",    OPTYPE_BINARY},
  {OPMULT,  "*",    OPTYPE_BINARY},
  {OPDIV,   "/",    OPTYPE_BINARY},
  {OPREM,   "mod",  OPTYPE_BINARY},
  {OPPOW,   "^",    OPTYPE_BINARY},
  {OPLESS,  "less", OPTYPE_BINARY},
  { 7, "unknown"},
  { 8, "unknown"},
  { 9, "unknown"},
  {10, "unknown"},
  {MINLIST,  "min",     OPTYPE_VARARG},
  {MAXLIST,  "max",     OPTYPE_VARARG},
  {FLOOR,    "floor",   OPTYPE_UNARY},
  {CEIL,     "ceil",    OPTYPE_UNARY},
  {ABS,      "abs",     OPTYPE_UNARY},
  {OPUMINUS, "unary -", OPTYPE_UNARY},
  {17, "unknown"},
  {18, "unknown"},
  {19, "unknown"},
  {OPOR,  "||", OPTYPE_BINARY, LOGICAL},
  {OPAND, "&&", OPTYPE_BINARY, LOGICAL},
  {LT,    "<",  OPTYPE_BINARY, LOGICAL},
  {LE,    "<=", OPTYPE_BINARY, LOGICAL},
  {EQ,    "=",  OPTYPE_BINARY, LOGICAL},
  {25, "unknown"},
  {26, "unknown"},
  {27, "unknown"},
  {GE, ">=", OPTYPE_BINARY, LOGICAL},
  {GT, ">",  OPTYPE_BINARY, LOGICAL},
  {NE, "!=", OPTYPE_BINARY, LOGICAL},
  {31, "unknown"},
  {32, "unknown"},
  {33, "unknown"},
  {OPNOT, "!", OPTYPE_UNARY, LOGICAL},
  {OPIFnl, "if-then-else", OPTYPE_IF},
  {36, "unknown"},
  {OP_tanh,  "tanh",  OPTYPE_UNARY},
  {OP_tan,   "tan",   OPTYPE_UNARY},
  {OP_sqrt,  "sqrt",  OPTYPE_UNARY},
  {OP_sinh,  "sinh",  OPTYPE_UNARY},
  {OP_sin,   "sin",   OPTYPE_UNARY},
  {OP_log10, "log10", OPTYPE_UNARY},
  {OP_log,   "log",   OPTYPE_UNARY},
  {OP_exp,   "exp",   OPTYPE_UNARY},
  {OP_cosh,  "cosh",  OPTYPE_UNARY},
  {OP_cos,   "cos",   OPTYPE_UNARY},
  {OP_atanh, "atanh", OPTYPE_UNARY},
  {OP_atan2, "atan2", OPTYPE_BINARY},
  {OP_atan,  "atan",  OPTYPE_UNARY},
  {OP_asinh, "asinh", OPTYPE_UNARY},
  {OP_asin,  "asin",  OPTYPE_UNARY},
  {OP_acosh, "acosh", OPTYPE_UNARY},
  {OP_acos,  "acos",  OPTYPE_UNARY},
  {OPSUMLIST, "sum",  OPTYPE_SUM},
  {OPintDIV,  "div",  OPTYPE_BINARY},
  {OPprecision, "precision", OPTYPE_BINARY},
  {OPround,     "round",     OPTYPE_BINARY},
  {OPtrunc,     "trunc",     OPTYPE_BINARY},
  {OPCOUNT,     "count",           OPTYPE_COUNT},
  {OPNUMBEROF,  "numberof",        OPTYPE_COUNT},
  {OPNUMBEROFs, "string numberof", OPTYPE_COUNT, UNSUPPORTED},
  {OPATLEAST, "atleast", OPTYPE_BINARY, LOGICAL},
  {OPATMOST,  "atmost",  OPTYPE_BINARY, LOGICAL},
  {OPPLTERM, "pl term", OPTYPE_PLTERM},
  {OPIFSYM,  "string if-then-else", OPTYPE_IF, UNSUPPORTED},
  {OPEXACTLY,    "exactly",     OPTYPE_BINARY, LOGICAL},
  {OPNOTATLEAST, "not atleast", OPTYPE_BINARY, LOGICAL},
  {OPNOTATMOST,  "not atmost",  OPTYPE_BINARY, LOGICAL},
  {OPNOTEXACTLY, "not exactly", OPTYPE_BINARY, LOGICAL},
  {ANDLIST, "forall", OPTYPE_SUM, LOGICAL},
  {ORLIST,  "exists", OPTYPE_SUM, LOGICAL},
  {OPIMPELSE, "implies else", OPTYPE_IF, LOGICAL},
  {OP_IFF, "iff", OPTYPE_BINARY, LOGICAL},
  {OPALLDIFF, "alldiff", OPTYPE_COUNT, LOGICAL},
  {OP1POW, "1pow", OPTYPE_UNARY},
  {OP2POW, "^2",   OPTYPE_UNARY},
  {OPCPOW, "cpow", OPTYPE_UNARY},
  {OPFUNCALL, "function call", OPTYPE_FUNCALL, UNSUPPORTED},
  {OPNUM, "number", OPTYPE_NUMBER, LOGICAL},
  {OPHOL, "string", OPTYPE_STRING, UNSUPPORTED},
  {OPVARVAL, "variable", OPTYPE_VARIABLE}
};

TEST_F(ExprTest, Operators) {
  int known_ops = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  EXPECT_EQ(N_OPS, size);
  for (int i = 0; i < N_OPS; ++i) {
    int opcode = OP_INFO[i].code;
    const char *opname = OP_INFO[i].name;
    expr raw = {reinterpret_cast<efunc*>(opcode)};
    Expr e(MakeExpr(&raw));
    EXPECT_EQ(opcode, e.opcode());
    EXPECT_STREQ(opname, e.opname());
    EXPECT_EQ(OP_INFO[i].optype, e.optype()) << opname;
    if (strcmp(opname, "unknown"))
      ++known_ops;
  }
  EXPECT_EQ(68, known_ops);
}

TEST_F(ExprTest, ExprOpCodeOutOfRangeInAccessors) {
  const char *message =
      "Assertion .*IsOpCodeInRange\\(\\)";
  {
    expr raw = {};
    Expr e(MakeExpr(&raw));
    raw.op = reinterpret_cast<efunc*>(-1);
    EXPECT_DEBUG_DEATH(e.opname();, message);
    EXPECT_DEBUG_DEATH(e.optype();, message);
  }
  {
    expr raw = {};
    Expr e(MakeExpr(&raw));
    raw.op = reinterpret_cast<efunc*>(N_OPS);
    EXPECT_DEBUG_DEATH(e.opname();, message);
    EXPECT_DEBUG_DEATH(e.optype();, message);
  }
  {
    expr raw = {};
    Expr e(MakeExpr(&raw));
    raw.op = reinterpret_cast<efunc*>(777);
    EXPECT_DEBUG_DEATH(e.opname();, message);
    EXPECT_DEBUG_DEATH(e.optype();, message);
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
  const char *message = "Assertion .*IsValid\\(\\)";
  for (; i < N_OPS; ++i) {
    const OpInfo &info = OP_INFO[i];
    if (info.optype != 0 && (info.type == 0 || info.code == OPNUM)) {
      MakeNumericExpr(info.code);
      ++numeric_count;
    } else {
      EXPECT_DEBUG_DEATH(MakeNumericExpr(info.code);, message);
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
  const char *message = "Assertion .*IsValid\\(\\)";
  for (; i < N_OPS; ++i) {
    const OpInfo &info = OP_INFO[i];
    if ((info.type & LOGICAL) != 0) {
      MakeLogicalExpr(info.code);
      ++logical_count;
    } else {
      EXPECT_DEBUG_DEATH(MakeLogicalExpr(info.code);, message);
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
    "Assertion .*HasTypeOrNull\\(OPTYPE_UNARY\\)");
}

TEST_F(ExprTest, BinaryExpr) {
  BinaryExpr null_expr;
  EXPECT_FALSE(null_expr);
  NumericExpr lhs(AddNum(42)), rhs(AddNum(43));
  BinaryExpr e(AddBinary(OPDIV, lhs, rhs));
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_DEBUG_DEATH(AddBinary(OPUMINUS, lhs, rhs);,
    "Assertion .*HasTypeOrNull\\(OPTYPE_BINARY\\)");
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
    "Assertion .*HasTypeOrNull\\(OPTYPE_VARARG\\)");
}

// TODO: more tests
}

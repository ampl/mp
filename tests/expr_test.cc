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

#include "gtest/gtest.h"
#include "solvers/util/expr.h"

using ampl::Expr;
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

namespace {

TEST(ExprTest, ExprCtor) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = {reinterpret_cast<efunc*>(42)};
  Expr e2(&raw2);
  EXPECT_EQ(42, e2.opcode());
  expr raw3 = {reinterpret_cast<efunc*>(77)};
  Expr e3(&raw3);
  EXPECT_EQ(77, e3.opcode());
}

TEST(ExprTest, SafeBool) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = {reinterpret_cast<efunc*>(42)};
  Expr e2(&raw2);
  EXPECT_TRUE(e2);
}

struct OpInfo {
  int code;
  const char *name;
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
  {OPOR,  "||", OPTYPE_BINARY},
  {OPAND, "&&", OPTYPE_BINARY},
  {LT,    "<",  OPTYPE_BINARY},
  {LE,    "<=", OPTYPE_BINARY},
  {EQ,    "=",  OPTYPE_BINARY},
  {25, "unknown"},
  {26, "unknown"},
  {27, "unknown"},
  {GE, ">=", OPTYPE_BINARY},
  {GT, ">",  OPTYPE_BINARY},
  {NE, "!=", OPTYPE_BINARY},
  {31, "unknown"},
  {32, "unknown"},
  {33, "unknown"},
  {OPNOT, "!", OPTYPE_UNARY},
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
  {OPNUMBEROFs, "string numberof", OPTYPE_COUNT},
  {OPATLEAST, "atleast", OPTYPE_BINARY},
  {OPATMOST,  "atmost",  OPTYPE_BINARY},
  {OPPLTERM, "pl term", OPTYPE_PLTERM},
  {OPIFSYM,  "string if-then-else", OPTYPE_IF},
  {OPEXACTLY,    "exactly",     OPTYPE_BINARY},
  {OPNOTATLEAST, "not atleast", OPTYPE_BINARY},
  {OPNOTATMOST,  "not atmost",  OPTYPE_BINARY},
  {OPNOTEXACTLY, "not exactly", OPTYPE_BINARY},
  {ANDLIST, "forall", OPTYPE_SUM},
  {ORLIST,  "exists", OPTYPE_SUM},
  {OPIMPELSE, "implies else", OPTYPE_IF},
  {OP_IFF, "iff", OPTYPE_BINARY},
  {OPALLDIFF, "alldiff", OPTYPE_COUNT},
  {OP1POW, "1pow", OPTYPE_UNARY},
  {OP2POW, "^2",   OPTYPE_UNARY},
  {OPCPOW, "cpow", OPTYPE_UNARY},
  {OPFUNCALL, "function call", OPTYPE_FUNCALL},
  {OPNUM, "number", OPTYPE_NUMBER},
  {OPHOL, "string", OPTYPE_STRING},
  {OPVARVAL, "variable", OPTYPE_VARIABLE}
};

TEST(ExprTest, Operators) {
  int known_ops = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  EXPECT_EQ(N_OPS, size);
  for (int i = 0; i < N_OPS; ++i) {
    int opcode = OP_INFO[i].code;
    const char *opname = OP_INFO[i].name;
    expr raw = {reinterpret_cast<efunc*>(opcode)};
    Expr e(&raw);
    EXPECT_EQ(opcode, e.opcode());
    EXPECT_STREQ(opname, e.opname());
    EXPECT_EQ(OP_INFO[i].type, e.type()) << opname;
    if (strcmp(opname, "unknown"))
      ++known_ops;
  }
  EXPECT_EQ(68, known_ops);
}

// TODO
}

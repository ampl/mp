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

namespace {

TEST(ExprTest, ExprCtor) {
  Expr e1;
  EXPECT_FALSE(e1);
  expr raw2 = {reinterpret_cast<efunc*>(42)};
  Expr e2(&raw2);
  EXPECT_EQ(42, e2.opcode());
  expr raw3 = {reinterpret_cast<efunc*>(777)};
  Expr e3(&raw3);
  EXPECT_EQ(777, e3.opcode());
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
};

const OpInfo OP_INFO[] = {
  {OPPLUS, "+"},
  {OPMINUS, "-"},
  {OPMULT, "*"},
  {OPDIV, "/"},
  {OPREM, "mod"},
  {OPPOW, "^"},
  {OPLESS, "less"},
  {MINLIST, "min"},
  {MAXLIST, "max"},
  {FLOOR, "floor"},
  {CEIL, "ceil"},
  {ABS, "abs"},
  {OPUMINUS, "unary -"},
  {OPOR, "||"},
  {OPAND, "&&"},
  {LT, "<"},
  {LE, "<="},
  {EQ, "="},
  {GE, ">="},
  {GT, ">"},
  {NE, "!="},
  {OPNOT, "!"},
  {OPIFnl, "if-then-else"},
  {OP_tanh, "tanh"},
  {OP_tan, "tan"},
  {OP_sqrt, "sqrt"},
  {OP_sinh, "sinh"},
  {OP_sin, "sin"},
  {OP_log10, "log10"},
  {OP_log, "log"},
  {OP_exp, "exp"},
  {OP_cosh, "cosh"},
  {OP_cos, "cos"},
  {OP_atanh, "atanh"},
  {OP_atan2, "atan2"},
  {OP_atan, "atan"},
  {OP_asinh, "asinh"},
  {OP_asin, "asin"},
  {OP_acosh, "acosh"},
  {OP_acos, "acos"},
  {OPSUMLIST, "sum"},
  {OPintDIV, "div"},
  {OPprecision, "precision"},
  {OPround, "round"},
  {OPtrunc, "trunc"},
  {OPCOUNT, "count"},
  {OPNUMBEROF, "numberof"},
  {OPNUMBEROFs, "string numberof"},
  {OPATLEAST, "atleast"},
  {OPATMOST, "atmost"},
  {OPPLTERM, "pl term"},
  {OPIFSYM, "string if-then-else"},
  {OPEXACTLY, "exactly"},
  {OPNOTATLEAST, "not atleast"},
  {OPNOTATMOST, "not atmost"},
  {OPNOTEXACTLY, "not exactly"},
  {ANDLIST, "forall"},
  {ORLIST, "exists"},
  {OPIMPELSE, "implies else"},
  {OP_IFF, "iff"},
  {OPALLDIFF, "alldiff"},
  {OP1POW, "1pow"},
  {OP2POW, "^2"},
  {OPCPOW, "cpow"},
  {OPFUNCALL, "function call"},
  {OPNUM, "number"},
  {OPHOL, "string"},
  {OPVARVAL, "variable"},
  {7, "unknown"},
  {8, "unknown"},
  {9, "unknown"},
  {10, "unknown"},
  {17, "unknown"},
  {18, "unknown"},
  {19, "unknown"},
  {25, "unknown"},
  {26, "unknown"},
  {27, "unknown"},
  {31, "unknown"},
  {32, "unknown"},
  {33, "unknown"},
  {36, "unknown"},
  {N_OPS, "unknown"},
  {-1, "unknown"},
  {500, "unknown"}
};

TEST(ExprTest, OpNames) {
  int known_ops = 0;
  int size = sizeof(OP_INFO) / sizeof(*OP_INFO);
  EXPECT_EQ(N_OPS + 3, size);
  for (int i = 0; i < size; ++i) {
    int opcode = OP_INFO[i].code;
    const char *opname = OP_INFO[i].name;
    expr raw = {reinterpret_cast<efunc*>(opcode)};
    Expr e(&raw);
    EXPECT_EQ(opcode, e.opcode());
    EXPECT_STREQ(opname, e.opname());
    if (strcmp(opname, "unknown"))
      ++known_ops;
  }
  EXPECT_EQ(68, known_ops);
}

// TODO
}

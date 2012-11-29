/*
 A C++ interface to AMPL expression trees.

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

#include "solvers/util/expr.h"

#include <sstream>
#include <cstdio>

using ampl::Expr;

namespace {
// An operation type.
// Numeric values for the operation types should be in sync with the ones in
// op_type.hd.
enum OpType {
  OPTYPE_UNARY    =  1,  // Unary operation
  OPTYPE_BINARY   =  2,  // Binary operation
  OPTYPE_VARARG   =  3,  // Variable-argument function such as min or max
  OPTYPE_PLTERM   =  4,  // Piecewise-linear term
  OPTYPE_IF       =  5,  // The if-then-else expression
  OPTYPE_SUM      =  6,  // The sum expression
  OPTYPE_FUNCALL  =  7,  // Function call
  OPTYPE_STRING   =  8,  // String
  OPTYPE_NUMBER   =  9,  // Number
  OPTYPE_VARIABLE = 10,  // Variable
  OPTYPE_COUNT    = 11   // The count expression
};
}

#ifdef HAVE_UNORDERED_MAP
std::size_t std::hash<ampl::Expr>::operator()(Expr expr) const {
  std::size_t hash = 0;
  ampl::HashCombine(hash, expr.opcode());

  struct expr *e = expr.expr_;
  switch (optype[expr.opcode()]) {
    case OPTYPE_UNARY:
      HashCombine(hash, Expr(e->L.e));
      break;

    case OPTYPE_BINARY:
      HashCombine(hash, Expr(e->L.e));
      HashCombine(hash, Expr(e->R.e));
      break;

    case OPTYPE_VARARG:
      for (de *d = reinterpret_cast<const expr_va*>(e)->L.d; d->e; d++)
        HashCombine(hash, Expr(d->e));
      break;

    case OPTYPE_PLTERM: {
      plterm *p = e->L.p;
      real *pce = p->bs;
      for (int i = 0, n = p->n * 2 - 1; i < n; i++)
        ampl::HashCombine(hash, pce[i]);
      HashCombine(hash, Expr(e->R.e));
      break;
    }

    case OPTYPE_IF: {
      const expr_if *eif = reinterpret_cast<const expr_if*>(e);
      HashCombine(hash, Expr(eif->e));
      HashCombine(hash, Expr(eif->T));
      HashCombine(hash, Expr(eif->F));
      break;
    }

    case OPTYPE_SUM:
    case OPTYPE_COUNT: {
      struct expr **ep = e->L.ep;
      for (; ep < e->R.ep; ep++)
        HashCombine(hash, Expr(*ep));
      break;
    }

    case OPTYPE_NUMBER:
      ampl::HashCombine(hash, reinterpret_cast<const expr_n*>(e)->v);
      break;

    case OPTYPE_VARIABLE:
      ampl::HashCombine(hash, e->a);
      break;

    default:
      throw ampl::UnsupportedExprError(expr.opname());
  }
  return hash;
}
#endif

namespace ampl {

const Expr::Kind Expr::KINDS[N_OPS] = {
    Expr::BINARY,  // OPPLUS
    Expr::BINARY,  // OPMINUS
    Expr::BINARY,  // OPMULT
    Expr::BINARY,  // OPDIV
    Expr::BINARY,  // OPREM
    Expr::BINARY,  // OPPOW
    Expr::BINARY,  // OPLESS
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::VARARG,  // MINLIST
    Expr::VARARG,  // MAXLIST
    Expr::UNARY,  // FLOOR
    Expr::UNARY,  // CEIL
    Expr::UNARY,  // ABS
    Expr::UNARY,  // OPUMINUS
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::BINARY_LOGICAL,  // OPOR
    Expr::BINARY_LOGICAL,  // OPAND
    Expr::RELATIONAL,  // LT
    Expr::RELATIONAL,  // LE
    Expr::RELATIONAL,  // EQ
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::RELATIONAL,  // GE
    Expr::RELATIONAL,  // GT
    Expr::RELATIONAL,  // NE
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::UNKNOWN,
    Expr::NOT,  // OPNOT
    Expr::IF,  // OPIFnl
    Expr::UNKNOWN,
    Expr::UNARY,  // OP_tanh
    Expr::UNARY,  // OP_tan
    Expr::UNARY,  // OP_sqrt
    Expr::UNARY,  // OP_sinh
    Expr::UNARY,  // OP_sin
    Expr::UNARY,  // OP_log10
    Expr::UNARY,  // OP_log
    Expr::UNARY,  // OP_exp
    Expr::UNARY,  // OP_cosh
    Expr::UNARY,  // OP_cos
    Expr::UNARY,  // OP_atanh
    Expr::BINARY,  // OP_atan2
    Expr::UNARY,  // OP_atan
    Expr::UNARY,  // OP_asinh
    Expr::UNARY,  // OP_asin
    Expr::UNARY,  // OP_acosh
    Expr::UNARY,  // OP_acos
    Expr::SUM,  // OPSUMLIST
    Expr::BINARY,  // OPintDIV
    Expr::BINARY,  // OPprecision
    Expr::BINARY,  // OPround
    Expr::BINARY,  // OPtrunc
    Expr::COUNT,  // OPCOUNT
    Expr::NUMBEROF,  // OPNUMBEROF
    Expr::UNKNOWN,  // OPNUMBEROFs - not supported yet
    Expr::LOGICAL_COUNT,  // OPATLEAST
    Expr::LOGICAL_COUNT,  // OPATMOST
    Expr::PLTERM,  // OPPLTERM
    Expr::UNKNOWN,  // OPIFSYM - not supported yet
    Expr::LOGICAL_COUNT,  // OPEXACTLY
    Expr::LOGICAL_COUNT,  // OPNOTATLEAST
    Expr::LOGICAL_COUNT,  // OPNOTATMOST
    Expr::LOGICAL_COUNT,  // OPNOTEXACTLY
    Expr::ITERATED_LOGICAL,  // ANDLIST
    Expr::ITERATED_LOGICAL,  // ORLIST
    Expr::IMPLICATION,  // OPIMPELSE
    Expr::BINARY_LOGICAL,  // OP_IFF
    Expr::ALLDIFF,  // OPALLDIFF
    Expr::BINARY,  // OP1POW
    Expr::UNARY,  // f_OP2POW
    Expr::BINARY,  // f_OPCPOW
    Expr::UNKNOWN,  // OPFUNCALL - not supported yet
    Expr::CONSTANT,  // OPNUM
    Expr::UNKNOWN,  // OPHOL - not supported yet
    Expr::VARIABLE  // OPVARVAL
};

// Operator names indexed by opcodes which are defined in opcode.hd.
const char *const Expr::OP_NAMES[N_OPS] = {
  "+",
  "-",
  "*",
  "/",
  "mod",
  "^",
  "less",
  "unknown",
  "unknown",
  "unknown",
  "unknown",
  "min",
  "max",
  "floor",
  "ceil",
  "abs",
  "unary -",
  "unknown",
  "unknown",
  "unknown",
  "||",
  "&&",
  "<",
  "<=",
  "=",
  "unknown",
  "unknown",
  "unknown",
  ">=",
  ">",
  "!=",
  "unknown",
  "unknown",
  "unknown",
  "!",
  "if-then-else",
  "unknown",
  "tanh",
  "tan",
  "sqrt",
  "sinh",
  "sin",
  "log10",
  "log",
  "exp",
  "cosh",
  "cos",
  "atanh",
  "atan2",
  "atan",
  "asinh",
  "asin",
  "acosh",
  "acos",
  "sum",
  "div",
  "precision",
  "round",
  "trunc",
  "count",
  "numberof",
  "string numberof",
  "atleast",
  "atmost",
  "pl term",
  "string if-then-else",
  "exactly",
  "not atleast",
  "not atmost",
  "not exactly",
  "forall",
  "exists",
  "implies else",
  "iff",
  "alldiff",
  "1pow",
  "^2",
  "cpow",
  "function call",
  "number",
  "string",
  "variable"
};

const de VarArgExpr::END = {0};

bool Equal(Expr expr1, Expr expr2) {
  if (expr1.opcode() != expr2.opcode())
    return false;

  expr *e1 = expr1.expr_;
  expr *e2 = expr2.expr_;
  switch (optype[expr1.opcode()]) {
    case OPTYPE_UNARY:
      return Equal(Expr(e1->L.e), Expr(e2->L.e));

    case OPTYPE_BINARY:
      return Equal(Expr(e1->L.e), Expr(e2->L.e)) &&
             Equal(Expr(e1->R.e), Expr(e2->R.e));

    case OPTYPE_VARARG: {
      de *d1 = reinterpret_cast<const expr_va*>(e1)->L.d;
      de *d2 = reinterpret_cast<const expr_va*>(e2)->L.d;
      for (; d1->e && d2->e; d1++, d2++)
        if (!Equal(Expr(d1->e), Expr(d2->e)))
          return false;
      return !d1->e && !d2->e;
    }

    case OPTYPE_PLTERM: {
      plterm *p1 = e1->L.p, *p2 = e2->L.p;
      if (p1->n != p2->n)
        return false;
      real *pce1 = p1->bs, *pce2 = p2->bs;
      for (int i = 0, n = p1->n * 2 - 1; i < n; i++) {
        if (pce1[i] != pce2[i])
          return false;
      }
      return Equal(Expr(e1->R.e), Expr(e2->R.e));
    }

    case OPTYPE_IF: {
      const expr_if *eif1 = reinterpret_cast<const expr_if*>(e1);
      const expr_if *eif2 = reinterpret_cast<const expr_if*>(e2);
      return Equal(Expr(eif1->e), Expr(eif2->e)) &&
             Equal(Expr(eif1->T), Expr(eif2->T)) &&
             Equal(Expr(eif1->F), Expr(eif2->F));
    }

    case OPTYPE_SUM:
    case OPTYPE_COUNT: {
      expr **ep1 = e1->L.ep;
      expr **ep2 = e2->L.ep;
      for (; ep1 < e1->R.ep && ep2 < e2->R.ep; ep1++, ep2++)
        if (!Equal(Expr(*ep1), Expr(*ep2)))
          return false;
      return ep1 == e1->R.ep && ep2 == e2->R.ep;
    }

    case OPTYPE_NUMBER:
      return reinterpret_cast<const expr_n*>(e1)->v ==
             reinterpret_cast<const expr_n*>(e2)->v;

    case OPTYPE_VARIABLE:
      return e1->a == e2->a;

    default:
      throw UnsupportedExprError(expr1.opname());
  }
}

std::string internal::FormatOpCode(Expr e) {
  char buffer[64];
  snprintf(buffer, sizeof(buffer), "%d", e.opcode());
  return buffer;
}

#ifdef HAVE_UNORDERED_MAP
std::size_t HashNumberOfArgs::operator()(const NumberOfExpr &e) const {
  std::size_t hash = 0;
  for (NumberOfExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    HashCombine(hash, static_cast<Expr>(*i));
  return hash;
}
#endif

bool EqualNumberOfArgs::operator()(
    const NumberOfExpr &lhs, const NumberOfExpr &rhs) const {
  if (lhs.num_args() != rhs.num_args())
    return false;
  for (NumberOfExpr::iterator i = lhs.begin(), end = lhs.end(),
       j = rhs.begin(); i != end; ++i, ++j) {
    if (!Equal(*i, *j))
      return false;
  }
  return true;
}
}

/*
 AMPL solver utilities.

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

#include "solvers/util/util.h"
#include "solvers/util/expr.h"

#include <cstddef>
#include <algorithm>
#include <sstream>

#include "asl.h"
#include "nlp.h"
#include "opcode.hd"

using std::size_t;

namespace {

struct OpInfo {
  int code;
  const char *name;
};

#define FUNOP(op) {OP_##op, #op}

// Information about operators sorted by code.
const OpInfo OP_INFO[] = {
  {OPPLUS,   "+"},
  {OPMINUS,  "-"},
  {OPMULT,   "*"},
  {OPDIV,    "/"},
  {OPREM,    "mod"},
  {OPPOW,    "^"},
  {OPLESS,   "less"},
  {MINLIST,  "min"},
  {MAXLIST,  "max"},
  {FLOOR,    "floor"},
  {CEIL,     "ceil"},
  {ABS,      "abs"},
  {OPUMINUS, "unary -"},
  {OPOR,     "||"},
  {OPAND,    "&&"},
  {LT,       "<"},
  {LE,       "<="},
  {EQ,       "="},
  {GE,       ">="},
  {GT,       ">"},
  {NE,       "!="},
  {OPNOT,    "!"},
  {OPIFnl,   "if-then-else"},
  FUNOP(tanh),
  FUNOP(tan),
  FUNOP(sqrt),
  FUNOP(sinh),
  FUNOP(sin),
  FUNOP(log10),
  FUNOP(log),
  FUNOP(exp),
  FUNOP(cosh),
  FUNOP(cos),
  FUNOP(atanh),
  FUNOP(atan2),
  FUNOP(atan),
  FUNOP(asinh),
  FUNOP(asin),
  FUNOP(acosh),
  FUNOP(acos),
  {OPSUMLIST,    "sum"},
  {OPintDIV,     "div"},
  {OPprecision,  "precision"},
  {OPround,      "round"},
  {OPtrunc,      "trunc"},
  {OPCOUNT,      "count"},
  {OPNUMBEROF,   "numberof"},
  {OPNUMBEROFs,  "string numberof"},
  {OPATLEAST,    "atleast"},
  {OPATMOST,     "atmost"},
  {OPPLTERM,     "pl term"},
  {OPIFSYM,      "string if-then-else"},
  {OPEXACTLY,    "exactly"},
  {OPNOTATLEAST, "not atleast"},
  {OPNOTATMOST,  "not atmost"},
  {OPNOTEXACTLY, "not exactly"},
  {ANDLIST,      "forall"},
  {ORLIST,       "exists"},
  {OPIMPELSE,    "implies else"},
  {OP_IFF,       "iff"},
  {OPALLDIFF,    "alldiff"},
  {OP1POW,       "1pow"},
  {OP2POW,       "^2"},
  {OPCPOW,       "cpow"},
  {OPFUNCALL,    "function call"},
  {OPNUM,        "number"},
  {OPHOL,        "string"},
  {OPVARVAL,     "variable"}
};

struct OpCodeLess {
  bool operator()(const OpInfo &lhs, const OpInfo &rhs) const {
    return lhs.code < rhs.code;
  }
};
}

namespace ampl {

const char *GetOpName(int opcode) {
  OpInfo info = { opcode, 0 };
  const OpInfo *end = OP_INFO + sizeof(OP_INFO) / sizeof(*OP_INFO);
  const OpInfo *p = std::lower_bound(OP_INFO, end, info, OpCodeLess());
  return p != end && p->code == opcode ? p->name : "unknown";
}

bool Equal(const expr *e1, const expr *e2) {
  size_t opnum = reinterpret_cast<size_t>(e1->op);
  if (opnum != reinterpret_cast<size_t>(e2->op))
    return false;

  int type = optype[opnum];
  switch (type) {
  case OPTYPE_UNARY:
    return Equal(e1->L.e, e2->L.e);

  case OPTYPE_BINARY:
    return Equal(e1->L.e, e2->L.e) && Equal(e1->R.e, e2->R.e);

  case OPTYPE_VARARG: {
    de *d1 = reinterpret_cast<const expr_va*>(e1)->L.d;
    de *d2 = reinterpret_cast<const expr_va*>(e2)->L.d;
    for (; d1->e && d2->e; d1++, d2++)
      if (!Equal(d1->e, d2->e))
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
    return Equal(e1->R.e, e2->R.e);
  }

  case OPTYPE_IF: {
    const expr_if *eif1 = reinterpret_cast<const expr_if*>(e1);
    const expr_if *eif2 = reinterpret_cast<const expr_if*>(e2);
    return Equal(eif1->e, eif2->e) && Equal(eif1->T, eif2->T)
        && Equal(eif1->F, eif2->F);
  }

  case OPTYPE_SUM:
  case OPTYPE_COUNT: {
    expr **ep1 = e1->L.ep;
    expr **ep2 = e2->L.ep;
    for (; ep1 < e1->R.ep && ep2 < e2->R.ep; ep1++, ep2++)
      if (!Equal(*ep1, *ep2))
        return false;
    return ep1 == e1->R.ep && ep2 == e2->R.ep;
  }

  case OPTYPE_FUNCALL:
  case OPTYPE_STRING:
    throw UnsupportedExprError(GetOpName(opnum));

  case OPTYPE_NUMBER:
    return reinterpret_cast<const expr_n*>(e1)->v
        == reinterpret_cast<const expr_n*>(e2)->v;

  case OPTYPE_VARIABLE:
    return e1->a == e2->a;

  default: {
    std::ostringstream oss;
    oss << "unknown operator type " << type;
    throw Error(oss.str());
  }
  }
}
}

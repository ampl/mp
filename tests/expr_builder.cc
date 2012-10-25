/*
 AMPL expression builder.

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

#include "tests/expr_builder.h"

namespace ampl {

ExprBuilder::~ExprBuilder() {
  for (std::vector<expr*>::const_iterator
       i = exprs_.begin(), end = exprs_.end(); i != end; ++i) {
    expr *e = *i;
    if (NumericExpr(e).opcode() >= N_OPS)
      continue;
    switch (NumericExpr(e).type()) {
      case OPTYPE_VARARG:
        delete reinterpret_cast<expr_va*>(e);
        break;
      case OPTYPE_PLTERM:
        std::free(e->L.p);
        delete e;
        break;
      case OPTYPE_IF:
        delete reinterpret_cast<expr_if*>(e);
        break;
      case OPTYPE_NUMBER:
        delete reinterpret_cast<expr_n*>(e);
        break;
      default:
        delete e;
        break;
    }
  }
  exprs_.clear();
}

NumericExpr ExprBuilder::NewVarArg(int opcode,
    NumericExpr e1, NumericExpr e2, NumericExpr e3) {
  expr_va e = {reinterpret_cast<efunc*>(opcode), 0, {0}, {0}, 0, 0, 0};
  expr_va *copy = new expr_va(e);
  expr *result(reinterpret_cast<expr*>(copy));
  de *args = new de[4];
  args[0] = MakeDE(e1);
  args[1] = MakeDE(e2);
  args[2] = MakeDE(e3);
  args[3] = MakeDE(NumericExpr());
  copy->L.d = args;
  return NewExpr(result);
}

NumericExpr ExprBuilder::NewPLTerm(
    int size, const double *args, int var_index) {
  expr e = {reinterpret_cast<efunc*>(OPPLTERM), 0, 0, {0}, {0}, 0};
  NumericExpr pl(NewExpr(new expr(e)));
  pl.expr_->L.p = static_cast<plterm*>(
      std::calloc(1, sizeof(plterm) + sizeof(real) * (size - 1)));
  pl.expr_->L.p->n = (size + 1) / 2;
  real *bs = pl.expr_->L.p->bs;
  for (int i = 0; i < size; i++)
    bs[i] = args[i];
  pl.expr_->R.e = NewVar(var_index).expr_;
  return pl;
}

NumericExpr ExprBuilder::NewIf(int opcode, NumericExpr condition,
    NumericExpr true_expr, NumericExpr false_expr) {
  expr_if e = {reinterpret_cast<efunc*>(opcode), 0, condition.expr_,
               true_expr.expr_, false_expr.expr_,
               0, 0, 0, 0, {0}, {0}, 0, 0};
  return NewExpr(reinterpret_cast<expr*>(new expr_if(e)));
}

NumericExpr ExprBuilder::NewSum(int opcode,
    NumericExpr arg1, NumericExpr arg2, NumericExpr arg3) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {0}, {0}, 0};
  NumericExpr sum(NewExpr(new expr(e)));
  expr** args = sum.expr_->L.ep = new expr*[3];
  sum.expr_->R.ep = args + (arg3.expr_ ? 3 : 2);
  args[0] = arg1.expr_;
  args[1] = arg2.expr_;
  args[2] = arg3.expr_;
  return sum;
}
}

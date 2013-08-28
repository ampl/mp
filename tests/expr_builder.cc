/*
 AMPL expression builder.

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

#include "tests/expr_builder.h"

namespace ampl {

expr *ExprBuilder::Call::Init(const char *name,
    const CallArg *arg_begin, const CallArg *arg_end) {
  name_ = name;
  info_.name = name_.c_str();
  int num_args = info_.nargs = arg_end - arg_begin;
  expr_.op = reinterpret_cast<efunc*>(OPFUNCALL);
  expr_.fi = &info_;
  expr_.al = &args_;
  args_.n = num_args;
  constants_.resize(num_args);
  args_.ra = &constants_[0];
  expr_args_.reserve(info_.nargs);
  for (int i = 0; i < num_args; ++i) {
    const CallArg &arg = arg_begin[i];
    constants_[i] = arg.constant();
    if (!arg.expr()) continue;
    argpair ap = {GetImpl(arg.expr())};
    ap.u.v = args_.ra + i;
    expr_args_.push_back(ap);
  }
  expr_.ap = &expr_args_[0];
  expr_.ape = expr_.ap + expr_args_.size();
  return reinterpret_cast<expr*>(&expr_);
}

ExprBuilder::~ExprBuilder() {
  for (std::vector<expr*>::const_iterator
       i = exprs_.begin(), end = exprs_.end(); i != end; ++i) {
    expr *e = *i;
    if (Expr(e).opcode() >= N_OPS)
      continue;
    switch (Expr(e).kind()) {
      case Expr::VARARG:
        delete reinterpret_cast<expr_va*>(e);
        break;
      case Expr::PLTERM:
        std::free(e->L.p);
        delete e;
        break;
      case Expr::IF:
        delete reinterpret_cast<expr_if*>(e);
        break;
      case Expr::CONSTANT:
        delete reinterpret_cast<expr_n*>(e);
        break;
      default:
        delete e;
        break;
    }
  }
  exprs_.clear();
}

VarArgExpr ExprBuilder::AddVarArg(int opcode,
    NumericExpr e1, NumericExpr e2, NumericExpr e3) {
  assert(e1);
  expr_va e = {reinterpret_cast<efunc*>(opcode), 0, {0}, {0}, 0, 0, 0};
  expr_va *copy = new expr_va(e);
  expr *result(reinterpret_cast<expr*>(copy));
  de *args = new de[4];
  args[0] = MakeDE(e1);
  args[1] = MakeDE(e2);
  args[2] = MakeDE(e3);
  args[3] = MakeDE(NumericExpr());
  copy->L.d = args;
  return AddExpr<VarArgExpr>(result);
}

PiecewiseLinearTerm ExprBuilder::AddPLTerm(
    int size, const double *args, int var_index) {
  expr e = {reinterpret_cast<efunc*>(OPPLTERM), 0, 0, {0}, {0}, 0};
  PiecewiseLinearTerm pl(AddExpr<PiecewiseLinearTerm>(new expr(e)));
  pl.expr_->L.p = static_cast<plterm*>(
      std::calloc(1, sizeof(plterm) + sizeof(double) * (size - 1)));
  pl.expr_->L.p->n = (size + 1) / 2;
  double *bs = pl.expr_->L.p->bs;
  for (int i = 0; i < size; i++)
    bs[i] = args[i];
  pl.expr_->R.e = AddVar(var_index).expr_;
  return pl;
}

CallExpr ExprBuilder::AddCall(const char *func_name,
    const CallArg *arg_begin, const CallArg *arg_end) {
  calls_.resize(calls_.size() + 1);
  return Expr::Create<CallExpr>(Expr(
      calls_.back().Init(func_name, arg_begin, arg_end)));
}
}

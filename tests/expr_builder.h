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

#ifndef TESTS_EXPR_BUILDER_H_
#define TESTS_EXPR_BUILDER_H_

#include <vector>
#include "solvers/util/expr.h"

namespace ampl {

class ExprBuilder {
private:
  std::vector<expr*> exprs_;

  template <typename T>
  T NewExpr(expr *e) {
    exprs_.push_back(e);
    return T(e);
  }

public:
  ~ExprBuilder();

  // Creates an expression representing a number.
  NumericExpr NewNum(double n) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), n};
    return NewExpr<NumericExpr>(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates a logical constant.
  LogicalExpr NewLogicalConstant(bool value) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), value ? 1. : 0.};
    return NewExpr<LogicalExpr>(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates an expression representing a variable reference.
  NumericExpr NewVar(int var_index) {
    expr e = {reinterpret_cast<efunc*>(OPVARVAL), var_index, 0, {0}, {0}, 0};
    return NewExpr<NumericExpr>(new expr(e));
  }

  // Creates a unary expression.
  template <typename T>
  T NewUnary(int opcode, T arg) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {arg.expr_}, {0}, 0};
    return NewExpr<T>(new expr(e));
  }

  // Creates a binary expression.
  template <typename T>
  T NewBinary(int opcode, T lhs, T rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return NewExpr<T>(new expr(e));
  }

  // Creates a relational expression.
  LogicalExpr NewRelational(int opcode, NumericExpr lhs, NumericExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return NewExpr<LogicalExpr>(new expr(e));
  }

  static de MakeDE(NumericExpr e) {
    de result = {e.expr_, 0, {0}};
    return result;
  }

  // Creates a variable-argument expression with up to 3 arguments.
  NumericExpr NewVarArg(int opcode, NumericExpr e1,
      NumericExpr e2, NumericExpr e3 = NumericExpr());

  NumericExpr NewPLTerm(int size, const double *args, int var_index);

  // Creates an expression representing if-then-else.
  template <typename T>
  T NewIf(int opcode, LogicalExpr condition, T true_expr, T false_expr);

  // Creates an expression representing a sum with up to 3 arguments.
  template <typename Result, typename T>
  Result NewSum(int opcode, T arg1, T arg2, T arg3 = T());

  template <typename T>
  NumericExpr NewSum(int opcode, T arg1, T arg2, T arg3 = T()) {
    return NewSum<NumericExpr, T>(opcode, arg1, arg2, arg3);
  }

  LogicalExpr NewIterated(int opcode, LogicalExpr arg1,
      LogicalExpr arg2, LogicalExpr arg3 = LogicalExpr()) {
    return NewSum<LogicalExpr, LogicalExpr>(opcode, arg1, arg2, arg3);
  }
};

template <typename T>
T ExprBuilder::NewIf(int opcode,
    LogicalExpr condition, T true_expr, T false_expr) {
  expr_if e = {reinterpret_cast<efunc*>(opcode), 0, condition.expr_,
               true_expr.expr_, false_expr.expr_,
               0, 0, 0, 0, {0}, {0}, 0, 0};
  return NewExpr<T>(reinterpret_cast<expr*>(new expr_if(e)));
}

template <typename Result, typename T>
Result ExprBuilder::NewSum(int opcode, T arg1, T arg2, T arg3) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {0}, {0}, 0};
  Result sum(NewExpr<Result>(new expr(e)));
  expr** args = sum.expr_->L.ep = new expr*[3];
  sum.expr_->R.ep = args + (arg3.expr_ ? 3 : 2);
  args[0] = arg1.expr_;
  args[1] = arg2.expr_;
  args[2] = arg3.expr_;
  return sum;
}
}

#endif  // TESTS_EXPR_BUILDER_H_

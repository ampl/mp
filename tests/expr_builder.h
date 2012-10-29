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
  T AddExpr(expr *e) {
    exprs_.push_back(e);
    return Expr::Create<T>(Expr(e));
  }

  static de MakeDE(NumericExpr e) {
    de result = {e.expr_, 0, {0}};
    return result;
  }

public:
  ~ExprBuilder();

  // Adds a new unary numeric expression.
  UnaryExpr AddUnary(int opcode, NumericExpr arg) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {arg.expr_}, {0}, 0};
    return AddExpr<UnaryExpr>(new expr(e));
  }

  // Adds a new binary numeric expression.
  BinaryExpr AddBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return AddExpr<BinaryExpr>(new expr(e));
  }

  // Adds a new variable-argument expression with up to 3 arguments.
  VarArgExpr AddVarArg(int opcode, NumericExpr e1,
      NumericExpr e2, NumericExpr e3 = NumericExpr());

  NumericExpr NewPLTerm(int size, const double *args, int var_index);

  // Adds a new numeric constant and returns it.
  NumericExpr AddNum(double n) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), n};
    return AddExpr<NumericExpr>(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates an expression representing a variable reference.
  NumericExpr NewVar(int var_index) {
    expr e = {reinterpret_cast<efunc*>(OPVARVAL), var_index, 0, {0}, {0}, 0};
    return AddExpr<NumericExpr>(new expr(e));
  }

  LogicalExpr AddBinaryLogical(int opcode, LogicalExpr lhs, LogicalExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return LogicalExpr(AddExpr<LogicalExpr>(new expr(e)));
  }

  // Creates a logical constant.
  LogicalExpr NewLogicalConstant(bool value) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), value ? 1. : 0.};
    return AddExpr<LogicalExpr>(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates a relational expression.
  LogicalExpr NewRelational(int opcode, NumericExpr lhs, NumericExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return AddExpr<LogicalExpr>(new expr(e));
  }

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

  // Creates a logical NOT expression adding it to this builder.
  LogicalExpr NewNot(LogicalExpr arg) {
    expr e = {reinterpret_cast<efunc*>(OPNOT), 0, 0, {arg.expr_}, {0}, 0};
    return AddExpr<LogicalExpr>(new expr(e));
  }
};

template <typename T>
T ExprBuilder::NewIf(int opcode,
    LogicalExpr condition, T true_expr, T false_expr) {
  expr_if e = {reinterpret_cast<efunc*>(opcode), 0, condition.expr_,
               true_expr.expr_, false_expr.expr_,
               0, 0, 0, 0, {0}, {0}, 0, 0};
  return AddExpr<T>(reinterpret_cast<expr*>(new expr_if(e)));
}

template <typename Result, typename T>
Result ExprBuilder::NewSum(int opcode, T arg1, T arg2, T arg3) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {0}, {0}, 0};
  Result sum(AddExpr<Result>(new expr(e)));
  expr** args = sum.expr_->L.ep = new expr*[3];
  sum.expr_->R.ep = args + (arg3.expr_ ? 3 : 2);
  args[0] = arg1.expr_;
  args[1] = arg2.expr_;
  args[2] = arg3.expr_;
  return sum;
}
}

#endif  // TESTS_EXPR_BUILDER_H_

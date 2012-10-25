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

  NumericExpr NewExpr(expr *e) {
    exprs_.push_back(e);
    return NumericExpr(e);
  }

public:
  ~ExprBuilder();

  static expr *GetExprPtr(Expr e) {
    return e.expr_;
  }

  // Creates an ASL expression representing a number.
  NumericExpr NewNum(double n) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), n};
    return NewExpr(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates an ASL expression representing a variable.
  NumericExpr NewVar(int var_index) {
    expr e = {reinterpret_cast<efunc*>(OPVARVAL), var_index, 0, {0}, {0}, 0};
    return NewExpr(new expr(e));
  }

  // Creates an unary ASL expression.
  NumericExpr NewUnary(int opcode, NumericExpr arg) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {arg.expr_}, {0}, 0};
    return NewExpr(new expr(e));
  }

  // Creates a binary ASL expression.
  NumericExpr NewBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return NewExpr(new expr(e));
  }

  static de MakeDE(NumericExpr e) {
    de result = {e.expr_, 0, {0}};
    return result;
  }

  // Creates a variable-argument ASL expression with up to 3 arguments.
  NumericExpr NewVarArg(int opcode, NumericExpr e1,
      NumericExpr e2, NumericExpr e3 = NumericExpr());

  NumericExpr NewPLTerm(int size, const double *args, int var_index);

  // Creates an ASL expression representing if-then-else.
  NumericExpr NewIf(int opcode, NumericExpr condition,
      NumericExpr true_expr, NumericExpr false_expr);

  // Creates an ASL expression representing a sum with up to 3 arguments.
  NumericExpr NewSum(int opcode, NumericExpr arg1,
      NumericExpr arg2, NumericExpr arg3 = NumericExpr());
};
}

#endif  // TESTS_EXPR_BUILDER_H_

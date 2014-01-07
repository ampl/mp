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

#ifndef TESTS_EXPR_BUILDER_H_
#define TESTS_EXPR_BUILDER_H_

#include <deque>
#include <vector>
#include "solvers/util/expr.h"

namespace ampl {

class CallArg {
 private:
  double constant_;
  NumericExpr expr_;

 public:
  CallArg(double constant, NumericExpr e = NumericExpr())
  : constant_(constant), expr_(e) {}
  CallArg(NumericExpr e) : constant_(0), expr_(e) {}

  double constant() const { return constant_; }
  NumericExpr expr() const { return expr_; }
};

class ExprBuilder {
 private:
  std::vector<expr*> exprs_;

  template <typename T>
  T AddExpr(expr *e) {
    exprs_.push_back(e);
    return Expr::Create<T>(Expr(e));
  }

  static expr *GetImpl(Expr e) { return e.expr_; }

  class Call {
   private:
    std::string name_;
    func_info info_;
    arglist args_;
    std::vector<double> constants_;
    std::vector<argpair> expr_args_;
    expr_f expr_;

   public:
    Call() : info_(), args_(), expr_() {}

    expr *Init(const char *name,
        const CallArg *arg_begin, const CallArg *arg_end);
  };

  std::deque<Call> calls_;

  static de MakeDE(NumericExpr e) {
    de result = {e.expr_, 0, {0}};
    return result;
  }

  // Adds a new if-then-else expression.
  template <typename Result, typename Arg>
  Result AddIf(int opcode, LogicalExpr condition,
      Arg true_expr, Arg false_expr);

  // Adds a new iterated expression with up to 3 arguments.
  template <typename Result, typename Arg>
  Result AddIterated(int opcode, Arg arg1, Arg arg2, Arg arg3 = Arg());

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
      NumericExpr e2 = NumericExpr(), NumericExpr e3 = NumericExpr());

  // Adds a new sum expression with up to 3 arguments.
  SumExpr AddSum(NumericExpr arg1 = NumericExpr(),
      NumericExpr arg2 = NumericExpr(), NumericExpr arg3 = NumericExpr()) {
    return AddIterated<SumExpr, NumericExpr>(OPSUMLIST, arg1, arg2, arg3);
  }

  // Adds a new count expression with up to 3 arguments.
  CountExpr AddCount(LogicalExpr arg1, LogicalExpr arg2,
      LogicalExpr arg3 = LogicalExpr()) {
    return AddIterated<CountExpr, LogicalExpr>(OPCOUNT, arg1, arg2, arg3);
  }

  // Adds a new if-then-else expression.
  IfExpr AddIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return AddIf<IfExpr, NumericExpr>(OPIFnl, condition, true_expr, false_expr);
  }

  // Adds a new piecewise-linear term.
  PiecewiseLinearExpr AddPL(int size, const double *args, int var_index);

  // Adds a new numeric constant.
  NumericConstant AddNum(double n) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), n};
    return AddExpr<NumericConstant>(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Adds a new variable reference.
  Variable AddVar(int var_index) {
    expr e = {reinterpret_cast<efunc*>(OPVARVAL), var_index, 0, {0}, {0}, 0};
    return AddExpr<Variable>(new expr(e));
  }

  // Adds a new numberof expression with up to 2 arguments.
  NumberOfExpr AddNumberOf(NumericExpr value, NumericExpr arg1,
      NumericExpr arg2 = NumericExpr()) {
    return AddIterated<NumberOfExpr, NumericExpr>(
        OPNUMBEROF, value, arg1, arg2);
  }

  // Adds a new call expression with up to 3 arguments.
  CallExpr AddCall(const char *func_name,
      const CallArg* arg_begin, const CallArg *arg_end);
  CallExpr AddCall(const char *func_name) { return AddCall(func_name, 0, 0); }
  CallExpr AddCall(const char *func_name, const CallArg &arg);
  CallExpr AddCall(const char *func_name,
        const CallArg &arg1, const CallArg &arg2);
  CallExpr AddCall(const char *func_name,
      const CallArg &arg1, const CallArg &arg2, const CallArg &arg3);
  CallExpr AddCall(const char *func_name,
      const CallArg &arg1, const CallArg &arg2,
      const CallArg &arg3, const CallArg &arg4);
  CallExpr AddCall(const char *func_name,
      const CallArg &arg1, const CallArg &arg2, const CallArg &arg3,
      const CallArg &arg4, const CallArg &arg5);

  // Adds a new alldiff expression with up to 3 arguments.
  AllDiffExpr AddAllDiff(NumericExpr arg1, NumericExpr arg2,
      NumericExpr arg3 = NumericExpr()) {
    return AddIterated<AllDiffExpr, NumericExpr>(OPALLDIFF, arg1, arg2, arg3);
  }

  // Adds a new logical constant.
  LogicalConstant AddBool(bool value) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), value ? 1. : 0.};
    return AddExpr<LogicalConstant>(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Adds a new relational expression.
  RelationalExpr AddRelational(int opcode, NumericExpr lhs, NumericExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return AddExpr<RelationalExpr>(new expr(e));
  }

  // Adds a new logical NOT expression.
  NotExpr AddNot(LogicalExpr arg) {
    expr e = {reinterpret_cast<efunc*>(OPNOT), 0, 0, {arg.expr_}, {0}, 0};
    return AddExpr<NotExpr>(new expr(e));
  }

  // Adds a new logical count expression.
  LogicalCountExpr AddLogicalCount(
      int opcode, NumericExpr value, CountExpr count) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {value.expr_}, {count.expr_}, 0};
    return AddExpr<LogicalCountExpr>(new expr(e));
  }

  // Adds a new binary logical expression.
  BinaryLogicalExpr AddBinaryLogical(
      int opcode, LogicalExpr lhs, LogicalExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.expr_}, {rhs.expr_}, 0};
    return AddExpr<BinaryLogicalExpr>(new expr(e));
  }

  // Adds a new implication expression.
  ImplicationExpr AddImplication(LogicalExpr condition,
      LogicalExpr true_expr, LogicalExpr false_expr) {
    return AddIf<ImplicationExpr, LogicalExpr>(
        OPIMPELSE, condition, true_expr, false_expr);
  }

  // Adds a new iterated logical expression.
  IteratedLogicalExpr AddIteratedLogical(int opcode,
      LogicalExpr arg1, LogicalExpr arg2, LogicalExpr arg3 = LogicalExpr()) {
    return AddIterated<IteratedLogicalExpr, LogicalExpr>(
        opcode, arg1, arg2, arg3);
  }
};

template <typename Result, typename Arg>
Result ExprBuilder::AddIterated(int opcode, Arg arg1, Arg arg2, Arg arg3) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {0}, {0}, 0};
  Result sum(AddExpr<Result>(new expr(e)));
  expr** args = sum.expr_->L.ep = new expr*[3];
  args[0] = arg1.expr_;
  args[1] = arg2.expr_;
  args[2] = arg3.expr_;
  sum.expr_->R.ep = args;
  for (int i = 2; i >= 0; --i) {
    if (args[i]) {
      sum.expr_->R.ep = args + i + 1;
      break;
    }
  }
  return sum;
}

template <typename Result, typename Arg>
Result ExprBuilder::AddIf(int opcode,
    LogicalExpr condition, Arg true_expr, Arg false_expr) {
  expr_if e = {reinterpret_cast<efunc*>(opcode), 0, condition.expr_,
               true_expr.expr_, false_expr.expr_,
               0, 0, 0, 0, {0}, {0}, 0, 0};
  return AddExpr<Result>(reinterpret_cast<expr*>(new expr_if(e)));
}
}

#endif  // TESTS_EXPR_BUILDER_H_

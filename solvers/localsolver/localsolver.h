/*
 AMPL solver interface to LocalSolver.

 Copyright (C) 2014 AMPL Optimization Inc

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

#ifndef AMPL_SOLVERS_LOCALSOLVER_H
#define AMPL_SOLVERS_LOCALSOLVER_H

#include <cmath>
#include <localsolver.h>
#include "solvers/util/solver.h"

namespace ampl {

namespace ls = localsolver;

// Converter of optimization problems from NL to LocalSolver format.
class NLToLocalSolverConverter :
  public ExprVisitor<NLToLocalSolverConverter,
                     ls::LSExpression*, ls::LSExpression*> {
 private:
  ls::LSModel &model_;
  std::vector<ls::LSExpression*> vars_;

  template<typename Term>
  ls::LSExpression *ConvertExpr(LinearExpr<Term> linear, NumericExpr nonlinear);

  // Converts an unary expression.
  template <typename Expr>
  ls::LSExpression *ConvertUnary(ls::LSOperator op, Expr e) {
    return model_.createExpression(op, Visit(e.arg()));
  }

  // Converts a binary expression.
  template <typename Expr>
  ls::LSExpression *ConvertBinary(ls::LSOperator op, Expr e) {
    return model_.createExpression(op, Visit(e.lhs()), Visit(e.rhs()));
  }

  // Converts a variable argument expression (VarArgExpr or SumExpr).
  template <typename Expr>
  ls::LSExpression *ConvertVarArg(ls::LSOperator op, Expr e) {
    ls::LSExpression *result = model_.createExpression(op);
    for (typename Expr::iterator i = e.begin(); *i; ++i)
      result->addOperand(Visit(*i));
    return result;
  }

  static ls::lsint MakeConst(int value) { return value; }

public:
  NLToLocalSolverConverter(ls::LSModel &model) : model_(model) {}

  void Convert(const Problem &p);

  ls::LSExpression *const *vars() const { return &vars_[0]; }

  // The methods below perform conversion of AMPL NL expressions into
  // equivalent LocalSolver expressions. LocalSolver doesn't support the
  // following expressions/functions:
  // * hyperbolic functions
  // * atan, asin, acos, atan2
  // * piecewise linear
  // TODO
  // * log, log10, exp

  ls::LSExpression *VisitUnaryMinus(UnaryExpr e) {
    return model_.createExpression(ls::O_Sub, MakeConst(0), Visit(e.arg()));
  }
  ls::LSExpression *VisitPow2(UnaryExpr e) {
    return model_.createExpression(ls::O_Pow, Visit(e.arg()), MakeConst(2));
  }
  ls::LSExpression *VisitFloor(UnaryExpr e) {
    return ConvertUnary(ls::O_Floor, e);
  }
  ls::LSExpression *VisitCeil(UnaryExpr e) {
    return ConvertUnary(ls::O_Ceil, e);
  }
  ls::LSExpression *VisitAbs(UnaryExpr e) {
    return ConvertUnary(ls::O_Abs, e);
  }
  ls::LSExpression *VisitTan(UnaryExpr e) {
    return ConvertUnary(ls::O_Tan, e);
  }
  ls::LSExpression *VisitSqrt(UnaryExpr e) {
    return ConvertUnary(ls::O_Sqrt, e);
  }
  ls::LSExpression *VisitSin(UnaryExpr e) {
    return ConvertUnary(ls::O_Sin, e);
  }
  ls::LSExpression *VisitLog10(UnaryExpr e) {
    return model_.createExpression(
        ls::O_Div, ConvertUnary(ls::O_Log, e), M_LN10);
  }
  ls::LSExpression *VisitLog(UnaryExpr e) {
    return ConvertUnary(ls::O_Log, e);
  }
  ls::LSExpression *VisitExp(UnaryExpr e) {
    return ConvertUnary(ls::O_Exp, e);
  }
  ls::LSExpression *VisitCos(UnaryExpr e) {
    return ConvertUnary(ls::O_Cos, e);
  }

  ls::LSExpression *VisitPlus(BinaryExpr e) {
    return ConvertBinary(ls::O_Sum, e);
  }
  ls::LSExpression *VisitMinus(BinaryExpr e) {
    return ConvertBinary(ls::O_Sub, e);
  }
  ls::LSExpression *VisitMult(BinaryExpr e) {
    return ConvertBinary(ls::O_Prod, e);
  }
  ls::LSExpression *VisitDiv(BinaryExpr e) {
    return ConvertBinary(ls::O_Div, e);
  }
  ls::LSExpression *VisitRem(BinaryExpr e) {
    return ConvertBinary(ls::O_Mod, e);
  }
  ls::LSExpression *VisitPow(BinaryExpr e) {
    return ConvertBinary(ls::O_Pow, e);
  }
  ls::LSExpression *VisitPowConstExp(BinaryExpr e) {
    return model_.createExpression(ls::O_Pow, Visit(e.lhs()),
        Cast<NumericConstant>(e.rhs()).value());
  }
  ls::LSExpression *VisitPowConstBase(BinaryExpr e) {
    return model_.createExpression(ls::O_Pow,
        Cast<NumericConstant>(e.rhs()).value(), Visit(e.lhs()));
  }
  ls::LSExpression *VisitNumericLess(BinaryExpr e) {
    return model_.createExpression(ls::O_Max,
        ConvertBinary(ls::O_Sub, e), MakeConst(0));
  }
  ls::LSExpression *VisitIntDiv(BinaryExpr e) {
    ls::LSExpression *rem = VisitRem(e);
    return model_.createExpression(O_div,
        model_.createExpression(ls::O_Sub, rem->getOperand(0), rem),
        rem->getOperand(1));
  }

  // TODO
  /*ls::LSExpression *VisitRound(BinaryExpr e) {
    // round does nothing because Gecode supports only integer expressions.
    RequireZeroRHS(e, "round");
    return Visit(e.lhs());
  }

  ls::LSExpression *VisitTrunc(BinaryExpr e) {
    // trunc does nothing because Gecode supports only integer expressions.
    RequireZeroRHS(e, "trunc");
    return Visit(e.lhs());
  }

  ls::LSExpression *VisitCount(CountExpr e);

  ls::LSExpression *VisitNumberOf(NumberOfExpr e);*/

  ls::LSExpression *VisitMin(VarArgExpr e) {
    return ConvertVarArg(ls::O_Min, e);
  }
  ls::LSExpression *VisitMax(VarArgExpr e) {
    return ConvertVarArg(ls::O_Max, e);
  }

  template <typename Expr>
  ls::LSExpression *VisitIf(Expr e) {
    return model_.createExpression(ls::O_If, Visit(e.condition()),
        Visit(e.true_expr()), Visit(e.false_expr()));
  }

  ls::LSExpression *VisitSum(SumExpr e) {
    return ConvertVarArg(ls::O_Sum, e);
  }

  ls::LSExpression *VisitNumericConstant(NumericConstant c) {
    return model_.createConstant(c.value());
  }

  ls::LSExpression *VisitVariable(Variable v) { return vars_[v.index()]; }

  ls::LSExpression *VisitNot(NotExpr e) {
    return ConvertUnary(ls::O_Not, e);
  }

  ls::LSExpression *VisitOr(BinaryLogicalExpr e) {
    return ConvertBinary(ls::O_Or, e);
  }

  ls::LSExpression *VisitAnd(BinaryLogicalExpr e) {
    return ConvertBinary(ls::O_And, e);
  }

  ls::LSExpression *VisitIff(BinaryLogicalExpr e) {
    return ConvertBinary(ls::O_Eq, e);
  }

  ls::LSExpression *VisitLess(RelationalExpr e) {
    return ConvertBinary(ls::O_Lt, e);
  }

  ls::LSExpression *VisitLessEqual(RelationalExpr e) {
    return ConvertBinary(ls::O_Leq, e);
  }

  ls::LSExpression *VisitEqual(RelationalExpr e) {
    return ConvertBinary(ls::O_Eq, e);
  }

  ls::LSExpression *VisitGreaterEqual(RelationalExpr e) {
    return ConvertBinary(ls::O_Geq, e);
  }

  ls::LSExpression *VisitGreater(RelationalExpr e) {
    return ConvertBinary(ls::O_Gt, e);
  }

  ls::LSExpression *VisitNotEqual(RelationalExpr e) {
    return ConvertBinary(ls::O_Neq, e);
  }

  ls::LSExpression *VisitForAll(IteratedLogicalExpr e) {
    return ConvertVarArg(ls::O_And, e);
  }

  ls::LSExpression *VisitExists(IteratedLogicalExpr e) {
    return ConvertVarArg(ls::O_Or, e);
  }

  ls::LSExpression *VisitImplication(ImplicationExpr e) {
    return VisitIf(e);
  }

  ls::LSExpression *VisitAllDiff(AllDiffExpr e);

  ls::LSExpression *VisitLogicalConstant(LogicalConstant c) {
    ls::lsint value = c.value();
    return model_.createConstant(value);
  }
};

class LocalSolver : public Solver {
 private:
  ls::LocalSolver solver_;
  int timelimit_;

  int GetTimeLimit(const SolverOption &) const {
    return timelimit_;
  }

  void SetTimeLimit(const SolverOption &opt, int value) {
    if (value <= 0)
      throw InvalidOptionValue(opt, value);
    timelimit_ = value;
  }

 protected:
  void DoSolve(Problem &p);

 public:
  LocalSolver();
};
}

#endif // AMPL_SOLVERS_LOCALSOLVER_H

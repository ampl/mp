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

#ifndef MP_SOLVERS_LOCALSOLVER_H_
#define MP_SOLVERS_LOCALSOLVER_H_

#include <localsolver.h>
#include "asl/aslsolver.h"
#include "asl/aslexpr-visitor.h"

namespace mp {

namespace ls = localsolver;

// Converter of optimization problems from NL to LocalSolver format.
class NLToLocalSolverConverter :
  public asl::ExprConverter<NLToLocalSolverConverter, ls::LSExpression> {
 private:
  ls::LSModel &model_;
  std::vector<ls::LSExpression> vars_;

  template<typename Term>
  ls::LSExpression ConvertExpr(asl::LinearExpr<Term> linear,
                               asl::NumericExpr nonlinear);

  // Converts an unary expression.
  template <typename Expr>
  ls::LSExpression ConvertUnary(ls::LSOperator op, Expr e) {
    return model_.createExpression(op, Visit(e.arg()));
  }

  // Converts a binary expression.
  template <typename Expr>
  ls::LSExpression ConvertBinary(ls::LSOperator op, Expr e) {
    return model_.createExpression(op, Visit(e.lhs()), Visit(e.rhs()));
  }

  // Converts a variable argument expression (VarArgExpr or SumExpr).
  template <typename Expr>
  ls::LSExpression ConvertVarArg(ls::LSOperator op, Expr e) {
    ls::LSExpression result = model_.createExpression(op);
    for (typename Expr::iterator i = e.begin(); *i; ++i)
      result.addOperand(Visit(*i));
    return result;
  }

  static ls::lsint MakeConst(int value) { return value; }

public:
  NLToLocalSolverConverter(ls::LSModel &model) : model_(model) {}

  void Convert(const ASLProblem &p);

  ls::LSExpression const *vars() const { return &vars_[0]; }

  // The methods below perform conversion of AMPL NL expressions into
  // equivalent LocalSolver expressions. LocalSolver doesn't support the
  // following expressions/functions:
  // * hyperbolic functions
  // * atan, asin, acos, atan2
  // * piecewise linear
  // TODO

  ls::LSExpression VisitMinus(asl::UnaryExpr e) {
    return model_.createExpression(ls::O_Sub, MakeConst(0), Visit(e.arg()));
  }
  ls::LSExpression VisitPow2(asl::UnaryExpr e) {
    return model_.createExpression(ls::O_Pow, Visit(e.arg()), MakeConst(2));
  }
  ls::LSExpression VisitFloor(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Floor, e);
  }
  ls::LSExpression VisitCeil(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Ceil, e);
  }
  ls::LSExpression VisitAbs(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Abs, e);
  }
  ls::LSExpression VisitTan(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Tan, e);
  }
  ls::LSExpression VisitSqrt(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Sqrt, e);
  }
  ls::LSExpression VisitSin(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Sin, e);
  }
  ls::LSExpression VisitLog10(asl::UnaryExpr e);
  ls::LSExpression VisitLog(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Log, e);
  }
  ls::LSExpression VisitExp(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Exp, e);
  }
  ls::LSExpression VisitCos(asl::UnaryExpr e) {
    return ConvertUnary(ls::O_Cos, e);
  }

  ls::LSExpression VisitAdd(asl::BinaryExpr e) {
    return ConvertBinary(ls::O_Sum, e);
  }
  ls::LSExpression VisitSub(asl::BinaryExpr e) {
    return ConvertBinary(ls::O_Sub, e);
  }
  ls::LSExpression VisitMul(asl::BinaryExpr e) {
    return ConvertBinary(ls::O_Prod, e);
  }
  ls::LSExpression VisitDiv(asl::BinaryExpr e) {
    return ConvertBinary(ls::O_Div, e);
  }
  ls::LSExpression VisitMod(asl::BinaryExpr e) {
    return ConvertBinary(ls::O_Mod, e);
  }
  ls::LSExpression VisitPow(asl::BinaryExpr e) {
    return ConvertBinary(ls::O_Pow, e);
  }
  ls::LSExpression VisitPowConstExp(asl::BinaryExpr e) {
    return model_.createExpression(ls::O_Pow, Visit(e.lhs()),
        asl::Cast<asl::NumericConstant>(e.rhs()).value());
  }
  ls::LSExpression VisitPowConstBase(asl::BinaryExpr e) {
    return model_.createExpression(ls::O_Pow,
        asl::Cast<asl::NumericConstant>(e.rhs()).value(), Visit(e.lhs()));
  }
  ls::LSExpression VisitLess(asl::BinaryExpr e) {
    return model_.createExpression(ls::O_Max,
        ConvertBinary(ls::O_Sub, e), MakeConst(0));
  }
  ls::LSExpression VisitIntDiv(asl::BinaryExpr e) {
    ls::LSExpression rem = VisitMod(e);
    return model_.createExpression(ls::O_Div,
        model_.createExpression(ls::O_Sub, rem.getOperand(0), rem),
        rem.getOperand(1));
  }

  // TODO
  /*ls::LSExpression VisitRound(BinaryExpr e) {
    // round does nothing because Gecode supports only integer expressions.
    RequireZeroRHS(e, "round");
    return Visit(e.lhs());
  }

  ls::LSExpression VisitTrunc(BinaryExpr e) {
    // trunc does nothing because Gecode supports only integer expressions.
    RequireZeroRHS(e, "trunc");
    return Visit(e.lhs());
  }

  ls::LSExpression VisitCount(CountExpr e);

  ls::LSExpression VisitNumberOf(NumberOfExpr e);*/

  ls::LSExpression VisitMin(asl::VarArgExpr e) {
    return ConvertVarArg(ls::O_Min, e);
  }
  ls::LSExpression VisitMax(asl::VarArgExpr e) {
    return ConvertVarArg(ls::O_Max, e);
  }

  template <typename Expr>
  ls::LSExpression VisitIf(Expr e) {
    return model_.createExpression(ls::O_If, Visit(e.condition()),
        Visit(e.true_expr()), Visit(e.false_expr()));
  }

  ls::LSExpression VisitSum(asl::SumExpr e) {
    return ConvertVarArg(ls::O_Sum, e);
  }

  ls::LSExpression VisitNumericConstant(asl::NumericConstant c) {
    return model_.createConstant(c.value());
  }

  ls::LSExpression VisitVariable(asl::Variable v) { return vars_[v.index()]; }

  ls::LSExpression VisitNot(asl::NotExpr e) {
    return ConvertUnary(ls::O_Not, e);
  }

  ls::LSExpression VisitOr(asl::BinaryLogicalExpr e) {
    return ConvertBinary(ls::O_Or, e);
  }

  ls::LSExpression VisitAnd(asl::BinaryLogicalExpr e) {
    return ConvertBinary(ls::O_And, e);
  }

  ls::LSExpression VisitIff(asl::BinaryLogicalExpr e) {
    return ConvertBinary(ls::O_Eq, e);
  }

  ls::LSExpression VisitLT(asl::RelationalExpr e) {
    return ConvertBinary(ls::O_Lt, e);
  }

  ls::LSExpression VisitLE(asl::RelationalExpr e) {
    return ConvertBinary(ls::O_Leq, e);
  }

  ls::LSExpression VisitEQ(asl::RelationalExpr e) {
    return ConvertBinary(ls::O_Eq, e);
  }

  ls::LSExpression VisitGE(asl::RelationalExpr e) {
    return ConvertBinary(ls::O_Geq, e);
  }

  ls::LSExpression VisitGT(asl::RelationalExpr e) {
    return ConvertBinary(ls::O_Gt, e);
  }

  ls::LSExpression VisitNE(asl::RelationalExpr e) {
    return ConvertBinary(ls::O_Neq, e);
  }

  ls::LSExpression VisitForAll(asl::IteratedLogicalExpr e) {
    return ConvertVarArg(ls::O_And, e);
  }

  ls::LSExpression VisitExists(asl::IteratedLogicalExpr e) {
    return ConvertVarArg(ls::O_Or, e);
  }

  ls::LSExpression VisitImplication(asl::ImplicationExpr e) {
    return VisitIf(e);
  }

  ls::LSExpression VisitAllDiff(asl::PairwiseExpr e);

  ls::LSExpression VisitLogicalConstant(asl::LogicalConstant c) {
    ls::lsint value = c.value();
    return model_.createConstant(value);
  }
};

class LocalSolver : public ASLSolver {
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
  void DoSolve(ASLProblem &p, SolutionHandler &sh);

 public:
  LocalSolver();
};
}

#endif  // MP_SOLVERS_LOCALSOLVER_H_

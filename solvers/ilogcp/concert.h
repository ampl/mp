/*
 AMPL to IBM/ILOG Concert interface.

 Copyright (C) 2013 AMPL Optimization Inc

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

#ifndef MP_SOLVERS_ILOGCP_CONCERT_H_
#define MP_SOLVERS_ILOGCP_CONCERT_H_

#ifdef __APPLE__
#include <limits.h>
#include <string.h>
#endif

#include <ilconcert/ilomodel.h>

#include <iostream>
#include <memory>

#include "mp/solver.h"
#include "asl/aslexpr-visitor.h"

namespace mp {

class Problem;

class NLToConcertConverter;

typedef asl::ExprConverter<
  NLToConcertConverter, IloExpr, IloConstraint> Converter;

// Converter of optimization problems from NL to Concert format.
class NLToConcertConverter : public Converter {
 private:
  IloEnv env_;
  IloModel model_;
  IloNumVarArray vars_;
  IloRangeArray cons_;
  int objno_;
  unsigned flags_;

  class CreateVar {
   private:
    IloEnv env_;

   public:
    CreateVar(IloEnv env) : env_(env) {}

    IloIntVar operator()() const {
      return IloIntVar(env_, IloIntMin, IloIntMax);
    }
  };

  typedef asl::NumberOfMap<IloIntVar, CreateVar> IlogNumberOfMap;
  IlogNumberOfMap numberofs_;

  IloNumExprArray ConvertArgs(asl::VarArgExpr e);
  IloIntVar ConvertArg(asl::CallExpr call, int index,
                       IloInt lb = IloIntMin, IloInt ub = IloIntMax);

  bool ConvertGlobalConstraint(asl::CallExpr expr, IloConstraint &con);

  // Converts a pairwise expression (alldiff or !alldiff).
  template <typename Constraint, bool negate>
  IloConstraint Convert(asl::PairwiseExpr e) {
    int n = e.num_args();
    std::vector<IloExpr> args(n);
    int index = 0;
    for (asl::PairwiseExpr::iterator
         i = e.begin(), end = e.end(); i != end; ++i) {
      args[index++] = Visit(*i);
    }
    Constraint alldiff(env_);
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j)
        alldiff.add(negate ? args[i] == args[j] : args[i] != args[j]);
    }
    return alldiff;
  }

 public:
  // Flags.
  enum {
    USENUMBEROF = 1,
    DEBUG       = 2
  };
  NLToConcertConverter(IloEnv env, unsigned flags);

  IloModel model() const { return model_; }
  IloNumVarArray vars() const { return vars_; }
  IloRangeArray cons() const { return cons_; }

  IloExpr Visit(asl::NumericExpr e) {
    if ((flags_ & DEBUG) != 0)
      fmt::print("{}\n", e.opstr());
    return Converter::Visit(e);
  }

  IloConstraint Visit(asl::LogicalExpr e) {
    if ((flags_ & DEBUG) != 0)
      fmt::print("{}\n", e.opstr());
    return Converter::Visit(e);
  }

  IloExpr VisitAdd(asl::BinaryExpr e) {
    return Visit(e.lhs()) + Visit(e.rhs());
  }

  IloExpr VisitSub(asl::BinaryExpr e) {
    return Visit(e.lhs()) - Visit(e.rhs());
  }

  IloExpr VisitMul(asl::BinaryExpr e) {
    return Visit(e.lhs()) * Visit(e.rhs());
  }

  IloExpr VisitDiv(asl::BinaryExpr e) {
    return Visit(e.lhs()) / Visit(e.rhs());
  }

  IloExpr VisitMod(asl::BinaryExpr e) {
    IloNumExpr lhs(Visit(e.lhs())), rhs(Visit(e.rhs()));
    return lhs - IloTrunc(lhs / rhs) * rhs;
  }

  IloExpr VisitPow(asl::BinaryExpr e) {
    return IloPower(Visit(e.lhs()), Visit(e.rhs()));
  }

  IloExpr VisitLess(asl::BinaryExpr e) {
    return IloMax(Visit(e.lhs()) - Visit(e.rhs()), 0.0);
  }

  IloExpr VisitMin(asl::VarArgExpr e) {
    return IloMin(ConvertArgs(e));
  }

  IloExpr VisitMax(asl::VarArgExpr e) {
    return IloMax(ConvertArgs(e));
  }

  IloExpr VisitMinus(asl::UnaryExpr e) {
    return -Visit(e.arg());
  }

  IloExpr VisitAbs(asl::UnaryExpr e) {
    return IloAbs(Visit(e.arg()));
  }

  IloExpr VisitFloor(asl::UnaryExpr e) {
    return IloFloor(Visit(e.arg()));
  }

  IloExpr VisitCeil(asl::UnaryExpr e) {
    return IloCeil(Visit(e.arg()));
  }

  IloExpr VisitIf(asl::IfExpr e);

  IloExpr VisitTanh(asl::UnaryExpr e) {
    IloNumExpr exp(IloExponent(2 * Visit(e.arg())));
    return (exp - 1) / (exp + 1);
  }

  IloExpr VisitTan(asl::UnaryExpr e) {
    return IloTan(Visit(e.arg()));
  }

  IloExpr VisitSqrt(asl::UnaryExpr e) {
    return IloPower(Visit(e.arg()), 0.5);
  }

  IloExpr VisitSinh(asl::UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return (IloExponent(arg) - IloExponent(-arg)) * 0.5;
  }

  IloExpr VisitSin(asl::UnaryExpr e) {
    return IloSin(Visit(e.arg()));
  }

  IloExpr VisitLog10(asl::UnaryExpr e) {
    return IloLog10(Visit(e.arg()));
  }

  IloExpr VisitLog(asl::UnaryExpr e) {
    return IloLog(Visit(e.arg()));
  }

  IloExpr VisitExp(asl::UnaryExpr e) {
    return IloExponent(Visit(e.arg()));
  }

  IloExpr VisitCosh(asl::UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return (IloExponent(arg) + IloExponent(-arg)) * 0.5;
  }

  IloExpr VisitCos(asl::UnaryExpr e) {
    return IloCos(Visit(e.arg()));
  }

  IloExpr VisitAtanh(asl::UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog((1 + arg) / (1 - arg)) * 0.5;
  }

  IloExpr VisitAtan2(asl::BinaryExpr e);

  IloExpr VisitAtan(asl::UnaryExpr e) {
    return IloArcTan(Visit(e.arg()));
  }

  IloExpr VisitAsinh(asl::UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog(arg + IloPower(IloSquare(arg) + 1, 0.5));
  }

  IloExpr VisitAsin(asl::UnaryExpr e) {
    return IloArcSin(Visit(e.arg()));
  }

  IloExpr VisitAcosh(asl::UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog(arg + IloPower(arg + 1, 0.5) * IloPower(arg - 1, 0.5));
  }

  IloExpr VisitAcos(asl::UnaryExpr e) {
    return IloArcCos(Visit(e.arg()));
  }

  IloExpr VisitSum(asl::SumExpr e);

  IloExpr VisitIntDiv(asl::BinaryExpr e) {
    return IloTrunc(Visit(e.lhs()) / Visit(e.rhs()));
  }

  IloExpr VisitRound(asl::BinaryExpr e);

  IloExpr VisitTrunc(asl::BinaryExpr e);

  IloExpr VisitCount(asl::CountExpr e);

  IloExpr VisitNumberOf(asl::NumberOfExpr e);

  IloExpr VisitPowConstExp(asl::BinaryExpr e) {
    return IloPower(Visit(e.lhs()),
                    asl::Cast<asl::NumericConstant>(e.rhs()).value());
  }

  IloExpr VisitPow2(asl::UnaryExpr e) {
    return IloSquare(Visit(e.arg()));
  }

  IloExpr VisitPowConstBase(asl::BinaryExpr e) {
    return IloPower(asl::Cast<asl::NumericConstant>(
                      e.lhs()).value(), Visit(e.rhs()));
  }

  IloExpr VisitPLTerm(asl::PiecewiseLinearExpr e);

  IloExpr VisitCall(asl::CallExpr e);

  IloExpr VisitNumericConstant(asl::NumericConstant n) {
    return IloExpr(env_, n.value());
  }

  IloExpr VisitVariable(asl::Variable v) {
    return vars_[v.index()];
  }

  IloConstraint VisitLogicalConstant(asl::LogicalConstant c) {
    return IloNumVar(env_, 1, 1) == c.value();
  }

  IloConstraint VisitLT(asl::RelationalExpr e) {
    return Visit(e.lhs()) < Visit(e.rhs());
  }

  IloConstraint VisitLE(asl::RelationalExpr e) {
    return Visit(e.lhs()) <= Visit(e.rhs());
  }

  IloConstraint VisitEQ(asl::RelationalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  IloConstraint VisitGE(asl::RelationalExpr e) {
    return Visit(e.lhs()) >= Visit(e.rhs());
  }

  IloConstraint VisitGT(asl::RelationalExpr e) {
    return Visit(e.lhs()) > Visit(e.rhs());
  }

  IloConstraint VisitNE(asl::RelationalExpr e) {
    return Visit(e.lhs()) != Visit(e.rhs());
  }

  IloConstraint VisitOr(asl::BinaryLogicalExpr e) {
    return IloIfThen(env_, !Visit(e.lhs()), Visit(e.rhs()));
  }

  IloConstraint VisitExists(asl::IteratedLogicalExpr e);

  IloConstraint VisitAnd(asl::BinaryLogicalExpr e) {
    return Visit(e.lhs()) && Visit(e.rhs());
  }

  IloConstraint VisitForAll(asl::IteratedLogicalExpr e);

  IloConstraint VisitNot(asl::NotExpr e) {
    return !Visit(e.arg());
  }

  IloConstraint VisitIff(asl::BinaryLogicalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  IloConstraint VisitImplication(asl::ImplicationExpr e);

  IloConstraint VisitAllDiff(asl::PairwiseExpr e) {
    return Convert<IloAnd, false>(e);
  }

  IloConstraint VisitNotAllDiff(asl::PairwiseExpr e) {
    return Convert<IloOr, true>(e);
  }

  // Combines 'numberof' operators into IloDistribute constraints
  // which are much more useful to the solution procedure.
  void FinishBuildingNumberOf();

  void Convert(const Problem &p);
};
}

#endif  // MP_SOLVERS_ILOGCP_CONCERT_H_

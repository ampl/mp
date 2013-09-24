/*
 Converter of optimization problems from NL to Concert format.

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

#include "solvers/ilogcp/converter.h"

namespace {
inline IloInt CastToInt(double value) {
  IloInt int_value = static_cast<int>(value);
  if (int_value != value) {
    throw ampl::Error(str(
        fmt::Format("value {} can't be represented as int") << value));
  }
  return int_value;
}

void RequireNonzeroConstRHS(ampl::BinaryExpr e, const std::string &func_name) {
  if (!IsZero(e.rhs())) {
    throw ampl::UnsupportedExprError::CreateFromExprString(
        func_name + " with nonzero second parameter");
  }
}
}

namespace ampl {

IloNumExprArray NLToConcertConverter::ConvertArgs(VarArgExpr e) {
  IloNumExprArray args(env_);
  for (VarArgExpr::iterator i = e.begin(); *i; ++i)
    args.add(Visit(*i));
  return args;
}

NLToConcertConverter::NLToConcertConverter(
    IloEnv env, IloNumVarArray vars, bool use_numberof, bool debug)
: env_(env), model_(env), vars_(vars), use_numberof_(use_numberof),
  debug_(debug), numberofs_(CreateVar(env)) {
}

IloExpr NLToConcertConverter::VisitIf(IfExpr e) {
  IloConstraint condition(Visit(e.condition()));
  IloNumVar var(env_, -IloInfinity, IloInfinity);
  model_.add(IloIfThen(env_, condition, var == Visit(e.true_expr())));
  model_.add(IloIfThen(env_, !condition, var == Visit(e.false_expr())));
  return var;
}

IloExpr NLToConcertConverter::VisitAtan2(BinaryExpr e) {
  IloNumExpr y(Visit(e.lhs())), x(Visit(e.rhs()));
  IloNumExpr atan(IloArcTan(y / x));
  IloNumVar result(env_, -IloInfinity, IloInfinity);
  model_.add(IloIfThen(env_, x >= 0, result == atan));
  model_.add(IloIfThen(env_, x <= 0 && y >= 0, result == atan + M_PI));
  model_.add(IloIfThen(env_, x <= 0 && y <= 0, result == atan - M_PI));
  return result;
}

IloExpr NLToConcertConverter::VisitSum(SumExpr e) {
  IloExpr sum(env_);
  for (SumExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += Visit(*i);
  return sum;
}

IloExpr NLToConcertConverter::VisitRound(BinaryExpr e) {
  RequireNonzeroConstRHS(e, "round");
  // Note that IloOplRound rounds half up.
  return IloOplRound(Visit(e.lhs()));
}

IloExpr NLToConcertConverter::VisitTrunc(BinaryExpr e) {
  RequireNonzeroConstRHS(e, "trunc");
  return IloTrunc(Visit(e.lhs()));
}

IloExpr NLToConcertConverter::VisitCount(CountExpr e) {
  IloExpr sum(env_);
  for (CountExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += Visit(*i);
  return sum;
}

IloExpr NLToConcertConverter::VisitNumberOf(NumberOfExpr e) {
  NumericExpr value = e.value();
  NumericConstant num = Cast<NumericConstant>(value);
  if (num && use_numberof_)
    return numberofs_.Add(num.value(), e);
  IloExpr sum(env_);
  IloExpr concert_value(Visit(value));
  for (NumberOfExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += (Visit(*i) == concert_value);
  return sum;
}

IloExpr NLToConcertConverter::VisitPiecewiseLinear(PiecewiseLinearExpr e) {
  IloNumArray slopes(env_), breakpoints(env_);
  int num_breakpoints = e.num_breakpoints();
  for (int i = 0; i < num_breakpoints; ++i) {
    slopes.add(e.slope(i));
    breakpoints.add(e.breakpoint(i));
  }
  slopes.add(e.slope(num_breakpoints));
  return IloPiecewiseLinear(vars_[e.var_index()], breakpoints, slopes, 0, 0);
}

IloConstraint NLToConcertConverter::VisitExists(IteratedLogicalExpr e) {
  IloOr disjunction(env_);
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i) {
    disjunction.add(Visit(*i));
  }
  return disjunction;
}

IloConstraint NLToConcertConverter::VisitForAll(IteratedLogicalExpr e) {
  IloAnd conjunction(env_);
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i) {
    conjunction.add(Visit(*i));
  }
  return conjunction;
}

IloConstraint NLToConcertConverter::VisitImplication(ImplicationExpr e) {
  IloConstraint condition(Visit(e.condition()));
  return IloIfThen(env_,  condition, Visit(e.true_expr())) &&
      IloIfThen(env_, !condition, Visit(e.false_expr()));
}

IloConstraint NLToConcertConverter::VisitAllDiff(AllDiffExpr e) {
  IloIntVarArray vars(env_);
  for (AllDiffExpr::iterator i = e.begin(), end = e.end(); i != end; ++i) {
    if (Variable v = Cast<Variable>(*i)) {
      vars.add(vars_[v.index()]);
      continue;
    }
    IloIntVar var(env_, IloIntMin, IloIntMax);
    model_.add(var == Visit(*i));
    vars.add(var);
  }
  return IloAllDiff(env_, vars);
}

void NLToConcertConverter::FinishBuildingNumberOf() {
  for (IlogNumberOfMap::iterator
      i = numberofs_.begin(), end = numberofs_.end(); i != end; ++i) {
    int index = 0;
    const IlogNumberOfMap::ValueMap &val_map = i->values;
    IloIntVarArray cards(env_, val_map.size());
    IloIntArray values(env_, val_map.size());
    for (IlogNumberOfMap::ValueMap::const_iterator j = val_map.begin(),
        val_end = val_map.end(); j != val_end; ++j, ++index) {
      values[index] = CastToInt(j->first);
      cards[index] = j->second;
    }

    index = 0;
    NumberOfExpr expr = i->expr;
    IloIntVarArray vars(env_, expr.num_args());
    for (NumberOfExpr::iterator
        j = expr.begin(), expr_end = expr.end(); j != expr_end; ++j, ++index) {
      IloIntVar var(env_, IloIntMin, IloIntMax);
      vars[index] = var;
      model_.add(var == Visit(*j));
    }

    model_.add(IloDistribute(env_, cards, values, vars));
  }
}
}

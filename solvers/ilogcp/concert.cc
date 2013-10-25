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

#include "solvers/ilogcp/concert.h"

#include <ilconcert/ilotupleset.h>

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

namespace {

inline IloInt ConvertToInt(double value) {
  IloInt int_value = static_cast<int>(value);
  if (int_value != value)
    ampl::ThrowError("value {} can't be represented as int") << value;
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

IloIntVar NLToConcertConverter::ConvertArg(
    const CallExpr::Arg &arg, IloInt lb, IloInt ub) {
  NumericExpr expr = arg.expr();
  if (!expr) {
    int constant = ConvertToInt(arg.constant());
    return IloIntVar(env_, constant, constant);
  }
  Variable var = Cast<Variable>(expr);
  if (var && arg.constant() == 0)
    return vars_[var.index()];
  IloExpr ilo_expr = Visit(expr);
  if (double constant = arg.constant())
    ilo_expr += constant;
  IloIntVar ilo_var(env_, lb, ub);
  model_.add(ilo_var == ilo_expr);
  return ilo_var;
}

NLToConcertConverter::NLToConcertConverter(IloEnv env, unsigned flags)
: env_(env), model_(env), vars_(env), cons_(env), flags_(flags),
  numberofs_(CreateVar(env)) {
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
  if (num && (flags_ & USENUMBEROF) != 0)
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

IloExpr NLToConcertConverter::VisitCall(CallExpr e) {
  const char *function_name = e.function().name();
  CallExpr::Args args(e);
  int num_args = e.num_args();
  if (std::strcmp(function_name, "element") == 0) {
    if (num_args < 2)
      ThrowError("{}: too few arguments") << function_name;
    CallExpr::Arg last_arg = args[num_args - 1];
    if (!last_arg.expr()) {
      // Index is constant - return the argument at specified index.
      IloInt index = ConvertToInt(last_arg.constant());
      if (index < 0 || index >= num_args - 1)
        ThrowError("{}: index {} is out of bounds") << function_name << index;
      CallExpr::Arg result = args[index];
      if (!result.expr())
        return IloExpr(env_, result.constant());
      IloExpr expr = Visit(result.expr());
      if (double constant = result.constant())
        expr += constant;
      return expr;
    }
    IloIntVar index_var = ConvertArg(last_arg, 0, num_args - 2);
    if (e.num_arg_exprs() > 1) {
      // Some elements are expressions - build IloIntVarArray.
      IloIntVarArray elements(env_, num_args - 1);
      for (int i = 0, n = num_args - 1; i < n; ++i)
        elements[i] = ConvertArg(args[i]);
      return IloElement(elements, index_var);
    }
    // All elements are constants - build IloNumArray.
    IloNumArray elements(env_, num_args - 1);
    for (int i = 0, n = num_args - 1; i < n; ++i) {
      assert(!args[i].expr());
      elements[i] = args[i].constant();
    }
    return IloElement(elements, index_var);
  } else if (std::strcmp(function_name, "in_relation") == 0)
    throw UnsupportedExprError::CreateFromExprString("nested 'in_relation'");
  throw UnsupportedExprError::CreateFromMessage(
      fmt::Format("unsupported function: {}") << function_name);
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
      values[index] = ConvertToInt(j->first);
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

bool NLToConcertConverter::ConvertGlobalConstraint(
    CallExpr expr, IloConstraint &con) {
  const char *function_name = expr.function().name();
  if (std::strcmp(function_name, "in_relation") != 0)
    return false;
  CallExpr::Args args(expr);
  int num_args = expr.num_args();
  if (num_args < 1)
    ThrowError("{}: too few arguments") << function_name;
  int arity = 0;
  for (; arity < num_args && args[arity].expr(); ++arity)
    ;  // Count variables.
  IloIntVarArray vars(env_, arity);
  for (int i = 0; i < arity; ++i)
    vars[i] = ConvertArg(args[i]);
  IloIntTupleSet set(env_, arity);
  if (num_args % arity != 0) {
    ThrowError("{}: the number of arguments {} is not a multiple of arity {}")
        << function_name << num_args << arity;
  }
  for (int i = arity; i < num_args; i += arity) {
    IloIntArray tuple(env_, arity);
    for (int j = 0; j < arity; ++j) {
      const CallExpr::Arg &arg = args[i + j];
      if (arg.expr()) {
        ThrowError("{}: argument {} is not constant")
          << function_name << (i + j + 1);
      }
      tuple[j] = ConvertToInt(arg.constant());
    }
    set.add(tuple);
  }
  con = IloAllowedAssignments(env_, vars, set);
  return true;
}

void NLToConcertConverter::Convert(const Problem &p) {
  int num_continuous_vars = p.num_continuous_vars();

  // Set up optimization problem using the Concert API.
  int num_vars = p.num_vars();
  vars_.setSize(num_vars);
  for (int j = 0; j < num_continuous_vars; j++)
    vars_[j] = IloNumVar(env_, p.var_lb(j), p.var_ub(j), ILOFLOAT);
  for (int j = num_continuous_vars; j < num_vars; j++)
    vars_[j] = IloNumVar(env_, p.var_lb(j), p.var_ub(j), ILOINT);

  int num_objs = p.num_objs();
  if (num_objs > 0) {
    if ((flags_ & MULTIOBJ) == 0)
      num_objs = 1;
    ObjType main_obj_type = p.obj_type(0);
    IloNumExprArray objs(env_);
    for (int i = 0; i < num_objs; ++i) {
      NumericExpr expr(p.nonlinear_obj_expr(i));
      NumericConstant constant(Cast<NumericConstant>(expr));
      IloExpr ilo_expr(env_, constant ? constant.value() : 0);
      if (p.num_nonlinear_objs() > 0 && !constant)
        ilo_expr += Visit(expr);
      LinearObjExpr linear = p.linear_obj_expr(i);
      for (LinearObjExpr::iterator
          j = linear.begin(), end = linear.end(); j != end; ++j) {
        ilo_expr += j->coef() * vars_[j->var_index()];
      }
      objs.add(p.obj_type(i) == main_obj_type ? ilo_expr : -ilo_expr);
    }
    IloObjective::Sense sense =
        main_obj_type == MIN ? IloObjective::Minimize : IloObjective::Maximize;
    IloAdd(model_, objs.getSize() == 1 ?
        IloObjective(env_, objs[0], sense) :
        IloObjective(env_, IloStaticLex(env_, objs), sense));
  }

  if (int n_cons = p.num_cons()) {
    cons_.setSize(n_cons);
    for (int i = 0; i < n_cons; ++i) {
      IloExpr expr(env_);
      LinearConExpr linear = p.linear_con_expr(i);
      for (LinearConExpr::iterator
          j = linear.begin(), end = linear.end(); j != end; ++j) {
        expr += j->coef() * vars_[j->var_index()];
      }
      if (i < p.num_nonlinear_cons())
        expr += Visit(p.nonlinear_con_expr(i));
      cons_[i] = (p.con_lb(i) <= expr <= p.con_ub(i));
    }
    model_.add(cons_);
  }

  if (int n_lcons = p.num_logical_cons()) {
    IloConstraintArray cons(env_, n_lcons);
    for (int i = 0; i < n_lcons; ++i) {
      LogicalExpr expr = p.logical_con_expr(i);
      if (expr.opcode() == NE) {
        RelationalExpr rel = Cast<RelationalExpr>(expr);
        NumericConstant const_rhs = Cast<NumericConstant>(rel.rhs());
        if (const_rhs && const_rhs.value() == 0) {
          CallExpr call = Cast<CallExpr>(rel.lhs());
          if (call && ConvertGlobalConstraint(call, cons[i]))
            continue;
        }
      }
      cons[i] = Visit(expr);
    }
    model_.add(cons);
  }

  FinishBuildingNumberOf();
}
}

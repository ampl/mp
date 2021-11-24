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

#include "concert.h"

#include <ilconcert/ilotupleset.h>

#include <algorithm>
#include <functional>

#include "mp/problem.h"
#include "mp/flat/backend.h"

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

using mp::Error;
using mp::NumericExpr;

namespace {

template <typename T>
inline T ConvertTo(double value) {
  T int_value = static_cast<T>(value);
  if (int_value != value)
    throw Error("value {} can't be represented as int", value);
  return int_value;
}

NumericExpr GetArg(mp::CallExpr e, int index) {
  NumericExpr result = mp::Cast<NumericExpr>(e.arg(index));
  if (!result) {
    throw Error("{}: argument {} is not numeric",
        e.function().name(), index + 1);
  }
  return result;
}

void RequireNonzeroConstRHS(mp::BinaryExpr e, fmt::StringRef func_name) {
  if (!IsZero(e.rhs())) {
    throw mp::MakeUnsupportedError(
          "{} with nonzero second parameter", func_name);
  }
}
}  // namespace

namespace mp {

bool EqualNumberOfArgs::operator()(IteratedExpr lhs, IteratedExpr rhs) const {
  int num_args = lhs.num_args();
  if (num_args != rhs.num_args())
    return false;
  for (int i = 1; i < num_args; ++i) {
    if (!Equal(lhs.arg(i), rhs.arg(i)))
      return false;
  }
  return true;
}

#ifdef MP_USE_HASH
size_t HashNumberOfArgs::operator()(IteratedExpr e) const {
  size_t hash = 0;
  for (int i = 1, n = e.num_args(); i < n; ++i)
    hash = internal::HashCombine<Expr>(hash, e.arg(i));
  return hash;
}
#elif defined(_MSC_VER)
# pragma message("warning: unordered_map not available, numberof may be slow")
#else
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wpedantic"
# warning "unordered_map not available, numberof may be slow"
# pragma clang diagnostic pop
#endif

IloNumExprArray MPToConcertConverter::ConvertArgs(VarArgExpr e) {
  IloNumExprArray args(env_);
  for (VarArgExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    args.add(Visit(*i));
  return args;
}

IloIntVar MPToConcertConverter::ConvertArg(
    CallExpr call, int index, IloInt lb, IloInt ub) {
  NumericExpr arg = GetArg(call, index);
  if (Variable var = Cast<Variable>(arg))
    return vars_[var.index()];
  IloIntVar ilo_var(env_, lb, ub);
  model_.add(ilo_var == Visit(arg));
  return ilo_var;
}

MPToConcertConverter::MPToConcertConverter(IloEnv env, unsigned flags)
: env_(env), model_(env), vars_(env), cons_(env), flags_(flags), objs_(env),
  nNumberofsFinalized_(0), numberofs_(CreateVar(env)) {
}

IloExpr MPToConcertConverter::VisitIf(IfExpr e) {
  IloConstraint condition(Visit(e.condition()));
  IloNumVar var(env_, -IloInfinity, IloInfinity);
  model_.add(IloIfThen(env_, condition, var == Visit(e.then_expr())));
  model_.add(IloIfThen(env_, !condition, var == Visit(e.else_expr())));
  return var;
}

IloExpr MPToConcertConverter::VisitAtan2(BinaryExpr e) {
  IloNumExpr y(Visit(e.lhs())), x(Visit(e.rhs()));
  IloNumExpr atan(IloArcTan(y / x));
  IloNumVar result(env_, -IloInfinity, IloInfinity);
  model_.add(IloIfThen(env_, x >= 0, result == atan));
  model_.add(IloIfThen(env_, x <= 0 && y >= 0, result == atan + M_PI));
  model_.add(IloIfThen(env_, x <= 0 && y <= 0, result == atan - M_PI));
  return result;
}

IloExpr MPToConcertConverter::VisitSum(SumExpr e) {
  IloExpr sum(env_);
  for (SumExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += Visit(*i);
  return sum;
}

IloExpr MPToConcertConverter::VisitRound(BinaryExpr e) {
  RequireNonzeroConstRHS(e, "round");
  // Note that IloOplRound rounds half up.
  return IloOplRound(Visit(e.lhs()));
}

IloExpr MPToConcertConverter::VisitTrunc(BinaryExpr e) {
  RequireNonzeroConstRHS(e, "trunc");
  return IloTrunc(Visit(e.lhs()));
}

IloExpr MPToConcertConverter::VisitCount(CountExpr e) {
  IloExpr sum(env_);
  for (CountExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += Visit(*i);
  return sum;
}

IloExpr MPToConcertConverter::VisitNumberOf(NumberOfExpr e) {
  NumericExpr value = e.arg(0);
  NumericConstant num = Cast<NumericConstant>(value);
  if (num && (flags_ & USENUMBEROF) != 0)
    return numberofs_.Add(num.value(), e);
  IloExpr sum(env_);
  IloExpr concert_value(Visit(value));
  for (int i = 1, n = e.num_args(); i < n; ++i)
    sum += (Visit(e.arg(i)) == concert_value);
  return sum;
}

IloExpr MPToConcertConverter::VisitPLTerm(PLTerm e) {
  IloNumArray slopes(env_), breakpoints(env_);
  int num_breakpoints = e.num_breakpoints();
  for (int i = 0; i < num_breakpoints; ++i) {
    slopes.add(e.slope(i));
    breakpoints.add(e.breakpoint(i));
  }
  slopes.add(e.slope(num_breakpoints));
  Variable var = Cast<Variable>(e.arg());
  return IloPiecewiseLinear(vars_[var.index()], breakpoints, slopes, 0, 0);
}

IloExpr MPToConcertConverter::VisitCall(CallExpr e) {
  const char *function_name = e.function().name();
  int num_args = e.num_args();
  if (std::strcmp(function_name, "element") == 0) {
    if (num_args < 2)
      throw Error("{}: too few arguments", function_name);
    if (NumericConstant num = Cast<NumericConstant>(e.arg(num_args - 1))) {
      // Index is constant - return the argument at specified index.
      int index = ConvertTo<int>(num.value());
      if (index < 0 || index >= num_args - 1)
        throw Error("{}: index {} is out of bounds", function_name, index);
      return Visit(GetArg(e, index));
    }
    IloIntVar index_var = ConvertArg(e, num_args - 1, num_args - 2);
    bool const_args = true;
    for (int i = 0, n = num_args - 1; i < n; ++i) {
      if (!Cast<NumericConstant>(e.arg(i))) {
        const_args = false;
        break;
      }
    }
    if (!const_args) {
      // Some elements are expressions - build IloIntVarArray.
      IloIntVarArray elements(env_, num_args - 1);
      for (int i = 0, n = num_args - 1; i < n; ++i)
        elements[i] = ConvertArg(e, i);
      return IloElement(elements, index_var);
    }
    // All elements are constants - build IloNumArray.
    IloNumArray elements(env_, num_args - 1);
    for (int i = 0, n = num_args - 1; i < n; ++i)
      elements[i] = Cast<NumericConstant>(e.arg(i)).value();
    return IloElement(elements, index_var);
  } else if (std::strcmp(function_name, "in_relation") == 0)
    throw MakeUnsupportedError("nested 'in_relation'");
  throw UnsupportedError("unsupported function: {}", function_name);
}

IloConstraint MPToConcertConverter::LogicalExprConverter::VisitExists(
    IteratedLogicalExpr e) {
  IloOr disjunction(converter_.env_);
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i) {
    disjunction.add(Visit(*i));
  }
  return disjunction;
}

IloConstraint MPToConcertConverter::LogicalExprConverter::VisitForAll(
    IteratedLogicalExpr e) {
  IloAnd conjunction(converter_.env_);
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i) {
    conjunction.add(Visit(*i));
  }
  return conjunction;
}

IloConstraint MPToConcertConverter::LogicalExprConverter::VisitImplication(
    ImplicationExpr e) {
  IloConstraint condition(Visit(e.condition()));
  return IloIfThen(converter_.env_,  condition, Visit(e.then_expr())) &&
      IloIfThen(converter_.env_, !condition, Visit(e.else_expr()));
}

void MPToConcertConverter::FinishBuildingNumberOf() {
  if (nNumberofsFinalized_ != 0)
    throw Error("Incrementality with numberofs not implemented");
  for (IlogNumberOfMap::iterator
      i = numberofs_.begin(), end = numberofs_.end(); i != end; ++i) {
    ++nNumberofsFinalized_;
    int index = 0;
    const IlogNumberOfMap::ValueMap &val_map = i->values;
    IloIntVarArray cards(env_, val_map.size());
    IloIntArray values(env_, val_map.size());
    for (IlogNumberOfMap::ValueMap::const_iterator j = val_map.begin(),
        val_end = val_map.end(); j != val_end; ++j, ++index) {
      values[index] = ConvertTo<IloInt>(j->first);
      cards[index] = j->second;
    }

    index = 0;
    NumberOfExpr expr = i->expr;
    int num_args = expr.num_args();
    IloIntVarArray vars(env_, num_args - 1);
    for (int j = 1; j < num_args; ++j) {
      IloIntVar var(env_, IloIntMin, IloIntMax);
      vars[j - 1] = var;
      model_.add(var == Visit(expr.arg(j)));
    }

    model_.add(IloDistribute(env_, cards, values, vars));
  }
}

bool MPToConcertConverter::ConvertGlobalConstraint(
    CallExpr expr, IloConstraint &con) {
  const char *function_name = expr.function().name();
  if (std::strcmp(function_name, "in_relation") != 0)
    return false;
  int num_args = expr.num_args();
  if (num_args < 1)
    throw Error("{}: too few arguments", function_name);
  int arity = 0;
  for (; arity < num_args && !Cast<NumericConstant>(expr.arg(arity)); ++arity)
    ;  // Count variables.
  IloIntVarArray vars(env_, arity);
  for (int i = 0; i < arity; ++i)
    vars[i] = ConvertArg(expr, i);
  IloIntTupleSet set(env_, arity);
  if (num_args % arity != 0) {
    throw Error("{}: the number of arguments {} is not a multiple of arity {}",
        function_name, num_args, arity);
  }
  for (int i = arity; i < num_args; i += arity) {
    IloIntArray tuple(env_, arity);
    for (int j = 0; j < arity; ++j) {
      NumericConstant num = Cast<NumericConstant>(expr.arg(i + j));
      if (!num) {
        throw Error("{}: argument {} is not constant",
            function_name, (i + j + 1));
      }
      tuple[j] = ConvertTo<IloInt>(num.value());
    }
    set.add(tuple);
  }
  con = IloAllowedAssignments(env_, vars, set);
  return true;
}

IloExpr MPToConcertConverter::ConvertExpr(
    const LinearExpr &linear, NumericExpr nonlinear) {
  IloExpr ilo_expr(env_, 0);
  if (nonlinear)
    ilo_expr = Visit(nonlinear);
  for (LinearExpr::const_iterator i = linear.begin(), e = linear.end(); i != e; ++i)
    ilo_expr += i->coef() * vars_[i->var_index()];
  return ilo_expr;
}

void MPToConcertConverter::Convert(const Problem &p) {
  // Set up optimization problem using the Concert API.
  p.PushModelTo(*this);
}

void MPToConcertConverter::FinishProblemModificationPhase() {
  FinishBuildingNumberOf();
  FinishBuildingObjectives();
  FinishBuildingAlgebraicConstraints();
}

void MPToConcertConverter::FinishBuildingObjectives() {
  IloObjective::Sense sense = main_obj_type_ == obj::MIN ?
      IloObjective::Minimize : IloObjective::Maximize;
  if (objs_.getSize()>0)
    IloAdd(model_, objs_.getSize() == 1 ?
      IloObjective(env_, objs_[0], sense) :
      IloObjective(env_, IloStaticLex(env_, objs_), sense));
}

void MPToConcertConverter::FinishBuildingAlgebraicConstraints() {
}

void MPToConcertConverter::AddVariable(Problem::Variable var) {
  vars_.add(IloNumVar(env_, var.lb(), var.ub(),
                      var.type() == mp::var::CONTINUOUS ? ILOFLOAT : ILOINT));
}

void MPToConcertConverter::AddCommonExpression(Problem::CommonExpr cexpr) {
  common_exprs_.push_back(
        ConvertExpr(cexpr.linear_expr(), cexpr.nonlinear_expr()) );
}

void MPToConcertConverter::AddObjective(Problem::Objective obj) {
  if (0==objs_.getSize())
    main_obj_type_ = obj.type();
  IloExpr ilo_expr = ConvertExpr(obj.linear_expr(), obj.nonlinear_expr());
  objs_.add(obj.type() == main_obj_type_ ? ilo_expr : -ilo_expr);
}

void MPToConcertConverter::AddAlgebraicConstraint(Problem::AlgebraicCon con) {
  IloExpr expr(env_);
  const LinearExpr &linear = con.linear_expr();
  for (LinearExpr::const_iterator
       j = linear.begin(), end = linear.end(); j != end; ++j) {
    expr += j->coef() * vars_[j->var_index()];
  }
  if (NumericExpr e = con.nonlinear_expr())
    expr += Visit(e);
  IloRange rng(env_, con.lb(), expr, con.ub());
  cons_.add(rng);
  model_.add(rng);
}

void MPToConcertConverter::AddLogicalConstraint(Problem::LogicalCon lcon) {
  IloConstraint cons;
  do {
    LogicalExpr expr = lcon.expr();
    if (expr.kind() == expr::NE) {
      RelationalExpr rel = Cast<RelationalExpr>(expr);
      NumericConstant const_rhs = Cast<NumericConstant>(rel.rhs());
      if (const_rhs && const_rhs.value() == 0) {
        CallExpr call = Cast<CallExpr>(rel.lhs());
        if (call && ConvertGlobalConstraint(call, cons))
          break;
      }
    } else if (expr.kind() == mp::expr::ALLDIFF) {
      PairwiseExpr alldiff = Cast<PairwiseExpr>(expr);
      // IloAllDiff is a global constraint that cannot be used
      // as a subexpression (the Concert API allows this, but
      // IloAlgorithm::extract throws CannotExtractException).
      IloIntVarArray vars(env_);
      for (PairwiseExpr::iterator
           j = alldiff.begin(), end = alldiff.end(); j != end; ++j) {
        if (j->kind() == expr::VARIABLE) {
          vars.add(vars_[Cast<Variable>(*j).index()]);
          continue;
        }
        IloIntVar var(env_, IloIntMin, IloIntMax);
        model_.add(var == Visit(*j));
        vars.add(var);
      }
      cons = IloAllDiff(env_, vars);
      break;
    }
    cons = Visit(expr);
  } while(0);
  model_.add(cons);
}

}  // namespace mp

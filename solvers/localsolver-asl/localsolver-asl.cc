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

#include "localsolver-asl/localsolver-asl.h"

#include <cmath>
#include "mp/clock.h"
#include "asl/aslproblem.h"

using namespace mp::asl;

namespace {
// Returns the value of an expression.
inline double GetValue(localsolver::LSExpression e) {
  return e.isDouble() ? e.getDoubleValue() : e.getValue();
}
}

namespace mp {

template <typename Term>
ls::LSExpression NLToLocalSolverConverter::ConvertExpr(
    LinearExpr<Term> linear, NumericExpr nonlinear) {
  ls::LSExpression result;
  typename LinearExpr<Term>::iterator i = linear.begin(), end = linear.end();
  bool has_linear_part = i != end;
  if (has_linear_part) {
    result = model_.createExpression(ls::O_Sum);
    for (; i != end; ++i) {
      ls::LSExpression term = vars_[i->var_index()];
      double coef = i->coef();
      if (coef != 1)
        term = model_.createExpression(ls::O_Prod, coef, term);
      result.addOperand(term);
    }
  }
  if (nonlinear) {
    ls::LSExpression nl = Visit(nonlinear);
    result = has_linear_part ?
        model_.createExpression(ls::O_Sum, result, nl) : nl;
  }
  return result;
}

void NLToLocalSolverConverter::Convert(const ASLProblem &p) {
  int num_vars = p.num_vars();
  vars_.resize(num_vars);

  // Convert continuous variables.
  int num_continuous_vars = p.num_continuous_vars();
  for (int j = 0; j < num_continuous_vars; ++j) {
    ls::LSExpression var =
        model_.createExpression(ls::O_Float, p.var_lb(j), p.var_ub(j));
    vars_[j] = var;
  }

  // Convert discrete variables.
  for (int j = num_continuous_vars; j < num_vars; j++) {
    // TODO: generate several bool vars for a general integer and handle bounds
    ls::LSExpression var = model_.createExpression(ls::O_Bool);
    vars_[j] = var;
  }

  // Convert objective.
  if (p.num_objs() != 0) {
    model_.addObjective(
        ConvertExpr(p.linear_obj_expr(0), p.nonlinear_obj_expr(0)),
        p.obj_type(0) == obj::MIN ? ls::OD_Minimize : ls::OD_Maximize);
  }

  // Convert constraints.
  for (int i = 0, n = p.num_cons(); i < n; ++i) {
    ls::LSExpression expr =
        ConvertExpr(p.linear_con_expr(i), p.nonlinear_con_expr(i));
    double lb = p.con_lb(i), ub = p.con_ub(i);
    if (lb <= negInfinity) {
      expr = model_.createExpression(ls::O_Leq, expr, ub);
    } else if (ub >= Infinity) {
      expr = model_.createExpression(ls::O_Geq, expr, lb);
    } else if (lb == ub) {
      expr = model_.createExpression(ls::O_Eq, expr, lb);
    } else {
      expr = model_.createExpression(ls::O_Geq, expr, lb);
      expr = model_.createExpression(ls::O_Leq, expr, ub);
    }
    model_.addConstraint(expr);
  }

  // Convert logical constraints.
/*  int num_logical_cons = p.num_logical_cons();
  for (int i = 0; i < num_logical_cons; ++i) {
    LogicalExpr e = p.logical_con_expr(i);
    AllDiffExpr alldiff = Cast<AllDiffExpr>(e);
    ICLSetter icl_setter(icl_, GetICL(p.num_cons() + i));
    if (!alldiff) {
      rel(problem_, Visit(e), icl_);
      continue;
    }
    int num_args = alldiff.num_args();
    IntVarArgs args(num_args);
    for (int i = 0; i < num_args; ++i) {
      NumericExpr arg(alldiff[i]);
      if (Variable var = ampl::Cast<Variable>(arg))
        args[i] = vars[var.index()];
      else
        args[i] = Gecode::expr(problem_, Visit(arg), icl_);
    }
    distinct(problem_, args, icl_);
  }*/
}

ls::LSExpression NLToLocalSolverConverter::VisitLog10(UnaryExpr e) {
  return model_.createExpression(
      ls::O_Div, ConvertUnary(ls::O_Log, e), std::log(10.0));
}

LocalSolver::LocalSolver()
  : ASLSolver("localsolver", 0, 20140710), timelimit_(0) {
  std::string version = fmt::format("{}.{}",
      localsolver::LSVersion::getMajorVersionNumber(),
      localsolver::LSVersion::getMinorVersionNumber());
  set_long_name("localsolver " + version);
  set_version("LocalSolver " + version);

  set_option_header(
      "LocalSolver Options for AMPL\n"
      "----------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to "
      "the AMPL option ``localsolver_options``. For example::\n"
      "\n"
      "  ampl: option localsolver_options 'version timelimit=30;\n");

  AddIntOption("timelimit",
      "Time limit in seconds (positive integer). Default = no limit.",
      &LocalSolver::GetTimeLimit, &LocalSolver::SetTimeLimit);
}

ls::LSExpression NLToLocalSolverConverter::VisitAllDiff(PairwiseExpr e) {
  ls::LSExpression result = model_.createExpression(ls::O_And);
  int num_args = e.num_args();
  std::vector<ls::LSExpression > args(num_args);
  for (int i = 0; i < num_args; ++i)
    args[i] = Visit(e[i]);
  for (int i = 0; i < num_args; ++i) {
    for (int j = i + 1; j < num_args; ++j)
      result.addOperand(model_.createExpression(ls::O_Neq, args[i], args[j]));
  }
  return result;
}

void LocalSolver::DoSolve(ASLProblem &p, SolutionHandler &sh) {
  steady_clock::time_point time = steady_clock::now();

  // Set up an optimization problem in LocalSolver.
  ls::LSModel model = solver_.getModel();
  NLToLocalSolverConverter converter(model);
  converter.Convert(p);
  model.close();

  // Set options. LS requires this to be done after the model is closed.
  ls::LSPhase phase = solver_.createPhase();
  if (timelimit_ != 0)
    phase.setTimeLimit(timelimit_);

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  solver_.solve();

  // Convert solution status.
  int solve_code = 0;
  ls::LSSolution sol = solver_.getSolution();
  const char *status = "unknown";
  switch (sol.getStatus()) {
  case ls::SS_Inconsistent:
    solve_code = sol::INFEASIBLE;
    status = "infeasible problem";
    break;
  case ls::SS_Infeasible:
    // Solution is infeasible, but problem may be feasible.
    // This can only happen if stopped by a limit.
    solve_code = sol::LIMIT;
    status = "infeasible solution";
    break;
  case ls::SS_Feasible:
    solve_code = sol::UNSOLVED;
    status = "feasible solution";
    break;
  case ls::SS_Optimal:
    solve_code = sol::SOLVED;
    status = "optimal solution";
    break;
  default:
    solve_code = sol::FAILURE;
    status = "unknown solution status";
    break;
  }

  int num_vars = p.num_vars();;
  ls::LSExpression const *vars = converter.vars();
  std::vector<double> solution(num_vars);
  int num_continuous_vars = p.num_continuous_vars();
  for (int i = 0; i < num_continuous_vars; ++i)
    solution[i] = vars[i].getDoubleValue();
  for (int i = num_continuous_vars; i < num_vars; ++i)
    solution[i] = vars[i].getValue();
  double solution_time = GetTimeAndReset(time);

  fmt::MemoryWriter w;
  w.write("{}: {}\n", long_name(), status);
  w.write("{}", solver_.getStatistics().toString());
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  if (p.num_objs() != 0) {
    obj_val = GetValue(model.getObjective(0));
    w.write("objective {}", FormatObjValue(obj_val));
  }
  sh.HandleSolution(solve_code, w.c_str(),
      solution.empty() ? 0 : solution.data(), 0, obj_val);
  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          setup_time, solution_time, output_time);
  }
}

SolverPtr CreateSolver(const char *) { return SolverPtr(new LocalSolver()); }
}

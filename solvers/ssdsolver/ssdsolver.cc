/*
 AMPL solver for problems with second-order stochastic dominance (SSD)
 constraints.

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

#include "ssdsolver/ssdsolver.h"
#include "asl/problem.h"

#ifdef _WIN32
# define putenv _putenv
#endif

namespace {

struct ValueScenario {
  double value;
  int scenario;
};

struct ValueLess {
  bool operator()(const ValueScenario &lhs, const ValueScenario &rhs) const {
    return lhs.value < rhs.value;
  }
};
}

namespace mp {

SSDSolver::SSDSolver() : ASLSolver("ssdsolver", 0, SSDSOLVER_VERSION),
  output_(false), scaled_(false), abs_tolerance_(1e-5), solver_name_("cplex") {
  set_version("SSD Solver");
  set_read_flags(Problem::READ_INITIAL_VALUES);
  AddIntOption("outlev", "0 or 1 (default 0):  Whether to print solution log.",
      &SSDSolver::GetBoolOption, &SSDSolver::SetBoolOption, &output_);
  AddIntOption("scaled", "0 or 1 (default 0):  Whether to use a scaled model.",
      &SSDSolver::GetBoolOption, &SSDSolver::SetBoolOption, &scaled_);
  AddDblOption("abs_tolerance", "Absolute tolerance. Default = 1e-5.",
      &SSDSolver::GetAbsTolerance, &SSDSolver::SetAbsTolerance);
  AddStrOption("solver", "Solver to use for subproblems. Default = cplex.",
      &SSDSolver::GetSolverName, &SSDSolver::SetSolverName);
}

void SSDSolver::DoSolve(Problem &p, SolutionHandler &sh) {
  Function ssd_uniform;
  int num_scenarios = p.num_logical_cons();
  int num_vars = p.num_vars();
  SSDExtractor extractor(num_scenarios, num_vars);
  for (int i = 0; i < num_scenarios; ++i) {
    LogicalExpr logical_expr = p.logical_con_expr(i);
    RelationalExpr rel_expr = Cast<RelationalExpr>(logical_expr);
    if (!rel_expr || rel_expr.kind() != expr::NE ||
        Cast<NumericConstant>(rel_expr.rhs()).value() != 0) {
      throw UnsupportedExprError::CreateFromExprString(logical_expr.opstr());
    }
    CallExpr call = Cast<CallExpr>(rel_expr.lhs());
    if (!call)
      throw UnsupportedExprError::CreateFromExprString(rel_expr.lhs().opstr());
    Function f = call.function();
    if (f == ssd_uniform)
      ; // Do nothing.
    else if (!ssd_uniform && std::strcmp(f.name(), "ssd_uniform") == 0)
      ssd_uniform = f;
    else
      throw UnsupportedExprError::CreateFromExprString(f.name());
    extractor.Extract(call);
  }

  if (p.num_objs() != 0)
    throw Error("SSD solver doesn't support user-defined objectives");

  ProblemChanges pc(p);
  int dominance_var = pc.AddVar(-Infinity, Infinity);
  double coef = 1;
  pc.AddObj(obj::MAX, 1, &coef, &dominance_var);

  // Compute the tails of the reference distribution.
  std::vector<double> ref_tails(extractor.rhs());
  std::sort(ref_tails.begin(), ref_tails.end());
  for (int i = 1; i < num_scenarios; ++i)
    ref_tails[i] += ref_tails[i - 1];

  // Get the initial solution.
  std::vector<double> solution;
  if (const double *initial_values = p.initial_values())
    solution.assign(initial_values, initial_values + num_vars);
  else
    solution.assign(num_vars, 0);

  // Disable solver output.
  char solver_msg[] = "solver_msg=0";
  putenv(solver_msg);

  // Solve the problem using a cutting-plane method.
  Solution sol;
  double dominance_lb = -Infinity;
  double dominance_ub =  Infinity;
  std::vector<double> cut_coefs(num_vars + 1);
  const double *coefs = extractor.coefs();
  std::vector<ValueScenario> tails(num_scenarios);
  int iteration = 1;
  printf("\nItn          Gap\n") ;
  for (; ; ++iteration) {
    // Compute the tails of the distribution.
    for (int i = 0; i < num_scenarios; ++i) {
      double value = 0;
      const double *row = coefs + i * num_vars;
      for (int j = 0; j < num_vars; ++j)
        value += row[j] * solution[j];
      tails[i].value = value;
      tails[i].scenario = i;
    }
    std::sort(tails.begin(), tails.end(), ValueLess());
    for (int i = 1; i < num_scenarios; ++i)
      tails[i].value += tails[i - 1].value;

    // Compute violation and minimal tail difference.
    double min_tail_diff = Infinity;
    double max_rel_violation = 0;
    int max_rel_violation_scen = -1;
    for (int i = 0; i < num_scenarios; ++i) {
      double scaling = scaled_ ? (i + 1.0) / num_scenarios : 1;
      double scaled_dominance = dominance_ub * scaling;
      double rel_violation =
          (scaled_dominance + ref_tails[i] + i + 1) / (tails[i].value + i + 1);
      if (rel_violation > max_rel_violation) {
        max_rel_violation = rel_violation;
        max_rel_violation_scen = i;
      }
      double tail_diff = (tails[i].value - ref_tails[i]) / scaling;
      if (tail_diff < min_tail_diff)
        min_tail_diff = tail_diff;
    }

    double scaling = scaled_ ?
        (max_rel_violation_scen + 1.0) / num_scenarios : 1;

    // Update the lower bound for the objective which by definition is a
    // minimum of tail differences (possibly scaled). Don't update in the
    // first iteration because the solution may not be feasible.
    if (min_tail_diff > dominance_lb && iteration != 1)
      dominance_lb = min_tail_diff;

    fmt::print("{:3} {:>12}\n", iteration, dominance_ub - dominance_lb);

    if ((dominance_ub - dominance_lb) * scaling <= abs_tolerance_) {
      fmt::print("Absolute tolerance reached.\n");
      break;
    }

    // Add a cut.
    for (int i = 0; i < num_vars; ++i) {
      double coef = 0;
      for (int j = 0; j <= max_rel_violation_scen; ++j)
        coef += coefs[tails[j].scenario * num_vars + i];
      cut_coefs[i] = coef;
    }
    cut_coefs[dominance_var] = -scaling;
    pc.AddCon(&cut_coefs[0], ref_tails[max_rel_violation_scen], Infinity);

    p.Solve(solver_name_, sol, &pc, Problem::IGNORE_FUNCTIONS);
    if (sol.status() != SOLVED) break;
    dominance_ub = sol.value(dominance_var);
    const double *values = sol.values();
    solution.assign(values, values + num_vars);
  }

  // Convert solution status.
  const char *message = 0;
  switch (sol.status()) {
  case SOLVED:
    message = "optimal solution";
    break;
  case INFEASIBLE:
    message = "infeasible problem";
    break;
  case UNBOUNDED:
    message = "unbounded problem";
    break;
  default:
    message = "error";
    break;
  }
  p.set_solve_code(sol.solve_code());

  fmt::Writer w;
  w.write("{}: {}", long_name(), message);
  if (sol.status() == SOLVED)
    w.write("; dominance {}", dominance_ub);
  w.write("\n{} iteration(s)", iteration);
  sh.HandleSolution(w.c_str(), solution.data(), 0, 0);
}

SolverPtr CreateSolver(const char *) { return SolverPtr(new SSDSolver()); }
}

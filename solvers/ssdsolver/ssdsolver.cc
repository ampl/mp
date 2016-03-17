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
#include "mp/problem.h"
#include "ilogcp/ilogcp.h"

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

SSDSolver::SSDSolver()
: SolverImpl<Problem>("ssdsolver", 0, SSDSOLVER_VERSION),
  output_(false), scaled_(false), abs_tolerance_(1e-5), solver_name_("cplex") {
  set_version("SSD Solver");
  AddIntOption("outlev", "0 or 1 (default 0):  Whether to print solution log.",
      &SSDSolver::GetBoolOption, &SSDSolver::SetBoolOption, &output_);
  AddIntOption("scaled", "0 or 1 (default 0):  Whether to use a scaled model.",
      &SSDSolver::GetBoolOption, &SSDSolver::SetBoolOption, &scaled_);
  AddDblOption("abs_tolerance", "Absolute tolerance. Default = 1e-5.",
      &SSDSolver::GetAbsTolerance, &SSDSolver::SetAbsTolerance);
  AddStrOption("solver", "Solver to use for subproblems. Default = cplex.",
      &SSDSolver::GetSolverName, &SSDSolver::SetSolverName);
}

void SSDSolver::Solve(Problem &p, SolutionHandler &outer_sh) {
  Function ssd_uniform;
  int num_scenarios = p.num_logical_cons();
  int num_vars = p.num_vars();
  SSDExtractor extractor(num_scenarios, num_vars);
  for (int i = 0; i < num_scenarios; ++i) {
    LogicalExpr logical_expr = p.logical_con(i).expr();
    RelationalExpr rel_expr = Cast<RelationalExpr>(logical_expr);
    if (!rel_expr || rel_expr.kind() != expr::NE ||
        Cast<NumericConstant>(rel_expr.rhs()).value() != 0) {
      throw MakeUnsupportedError(str(logical_expr.kind()));
    }
    CallExpr call = Cast<CallExpr>(rel_expr.lhs());
    if (!call)
      throw MakeUnsupportedError(str(rel_expr.lhs().kind()));
    Function f = call.function();
    if (f == ssd_uniform)
      ; // Do nothing.
    else if (!ssd_uniform && std::strcmp(f.name(), "ssd_uniform") == 0)
      ssd_uniform = f;
    else
      throw UnsupportedError("unsupported function: {}", f.name());
    extractor.Extract(call);
  }

  if (p.num_objs() != 0)
    throw Error("SSD solver doesn't support user-defined objectives");

  double inf = std::numeric_limits<double>::infinity();
  int dominance_var = p.AddVar(-inf, inf).index();
  p.AddObj(obj::MAX, 1).AddTerm(1, dominance_var);

  // Compute the tails of the reference distribution.
  std::vector<double> ref_tails(extractor.rhs());
  std::sort(ref_tails.begin(), ref_tails.end());
  for (int i = 1; i < num_scenarios; ++i)
    ref_tails[i] += ref_tails[i - 1];

  struct SSDSolutionHandler : BasicSolutionHandler {
    int num_vars;
    int status;
    std::vector<double> solution;
    explicit SSDSolutionHandler(int num_vars) : num_vars(num_vars), status(0) {}

    void HandleSolution(int status, fmt::CStringRef,
        const double *values, const double *, double) {
      this->status = status;
      solution.assign(values, values + num_vars + 1);
    }
  } sh(p.num_vars());

  // Get the initial solution.
  sh.solution.reserve(p.num_vars());
  auto vars = p.vars();
  for (auto i = vars.begin(), end = vars.end(); i != end; ++i)
    sh.solution.push_back(i->value());

  // Disable solver output.
  char solver_msg[] = "solver_msg=0";
  putenv(solver_msg);

  // TODO: convert problem

  // Solve the problem using a cutting-plane method.
  IlogCPSolver solver;
  solver.SetStrOption("optimizer", "cplex");
  double dominance_lb = -inf;
  double dominance_ub =  inf;
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
        value += row[j] * sh.solution[j];
      tails[i].value = value;
      tails[i].scenario = i;
    }
    std::sort(tails.begin(), tails.end(), ValueLess());
    for (int i = 1; i < num_scenarios; ++i)
      tails[i].value += tails[i - 1].value;

    // Compute violation and minimal tail difference.
    double min_tail_diff = inf;
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
    p.AddCon(ref_tails[max_rel_violation_scen], inf);
    Problem::LinearConBuilder cut = p.algebraic_con(p.num_algebraic_cons() - 1).
        set_linear_expr(num_vars + 1);
    for (int i = 0; i < num_vars; ++i) {
      double coef = 0;
      for (int j = 0; j <= max_rel_violation_scen; ++j)
        coef += coefs[tails[j].scenario * num_vars + i];
      cut.AddTerm(i, coef);
    }
    cut.AddTerm(dominance_var, -scaling);

    solver.Solve(p, sh);
    if (sh.status != sol::SOLVED) break;
    dominance_ub = sh.solution[dominance_var];
  }

  // Convert solution status.
  const char *message = 0;
  switch (sh.status) {
  case sol::SOLVED:
    message = "optimal solution";
    break;
  case sol::INFEASIBLE:
    message = "infeasible problem";
    break;
  case sol::UNBOUNDED:
    message = "unbounded problem";
    break;
  default:
    message = "error";
    break;
  }

  fmt::MemoryWriter w;
  w.write("{}: {}", long_name(), message);
  if (sh.status == sol::SOLVED)
    w.write("; dominance {}", dominance_ub);
  w.write("\n{} iteration(s)", iteration);
  outer_sh.HandleSolution(sh.status, w.c_str(), &sh.solution[0], 0, 0);
}

SolverPtr create_ssdsolver(const char *) { return SolverPtr(new SSDSolver()); }
}

/*
 IBM/ILOG CP solver for AMPL.

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

 Author: Victor Zverovich (based on the older version by Robert Fourer)

 October 2000: Linear/Nonlinear version (Robert Fourer)
 June 2012:    Updated to Concert 12.4 (Victor Zverovich)

 Possible improvements: Some sort of variable preference mechanism.

 Reference: "Extending an Algebraic Modeling Language to
 Support Constraint Programming" by Robert Fourer and David M. Gay,
 INFORMS Journal on Computing, Fall 2002, vol. 14, no. 4, 322-344
 (http://joc.journal.informs.org/content/14/4/322).
 */

#include "gurobi.h"

#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>

#include <iostream>

using std::strcmp;
using std::vector;

namespace {


const mp::OptionValueInfo OPTIMIZERS[] = {
  {
    "auto",
    "CP Optimizer if the problem has nonlinear objective/constraints "
    "or logical constraints, CPLEX otherwise", 0
  },
  {
    "cp",
    "CP Optimizer", 0
  },
  {
    "cplex",
    "CPLEX Optimizer", 0
  }
};


mp::OptionError GetOptionValueError(
    const mp::SolverOption &opt, fmt::StringRef message) {
  throw mp::OptionError(fmt::format(
      "Can't get value of option {}: {}", opt.name(), message));
}


bool HasNonlinearObj(const mp::Problem &p) {
  if (p.num_objs() == 0)
    return false;
  mp::NumericExpr expr = p.obj(0).nonlinear_expr();
  return expr && !mp::Cast<mp::NumericConstant>(expr);
}

std::string ConvertSolutionStatus(
    GRBmodel* model, const mp::Interrupter &interrupter, int &solve_code) {
  namespace sol = mp::sol;
  int optimstatus;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus) );
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter.Stop()) {
      solve_code = 600;
      return "interrupted";
    }
    int solcount;
    GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solcount) );
    if (solcount>0) {
      solve_code = sol::UNCERTAIN;
      return "feasible solution";
    }
    solve_code = sol::FAILURE + 1;
    return "unknown solution status";
  case GRB_OPTIMAL:
    solve_code = sol::SOLVED;
    return "optimal solution";
  case GRB_INFEASIBLE:
    solve_code = sol::INFEASIBLE;
    return "infeasible problem";
  case GRB_UNBOUNDED:
    solve_code = sol::UNBOUNDED;
    return "unbounded problem";
  case GRB_INF_OR_UNBD:
    solve_code = sol::INFEASIBLE + 1;
    return "infeasible or unbounded problem";
  case GRB_NUMERIC:
    solve_code = sol::FAILURE;
    return "error";
  }
}

bool InterruptGurobi(void *model) {
  GRBterminate( static_cast<GRBmodel*>(model) );
  return true;
}

}  // namespace

namespace mp {

GurobiBackend::GurobiBackend() :
   SolverImpl<Problem>("gurobi", 0, 0, MULTIPLE_SOL | MULTIPLE_OBJ)
   {
  InitBackend();

  options_[DEBUGEXPR] = false;
  options_[USENUMBEROF] = true;
  options_[SOLUTION_LIMIT] = -1;

  int a,b,c;
  GRBversion(&a, &b, &c);
  set_long_name(fmt::format("Gurobi {}.{}.{}", a, b, c));
  set_version(fmt::format("AMPL/Gurobi Optimizer [{}.{}.{}]", a,b,c));

  AddSuffix("priority", 0, suf::VAR);

  set_option_header(
      "Gurobi Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``ilogcp_options``. For example::\n"
      "\n"
      "  ampl: option ilogcp_options 'optimalitytolerance=1e-6 "
      "searchtype=restart';\n");

  AddStrOption("optimizer",
      "Specifies which optimizer to use. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``auto``.",
      &GurobiBackend::GetOptimizer, &GurobiBackend::SetOptimizer, OPTIMIZERS);

}

GurobiBackend::~GurobiBackend() {
  CloseBackend();
}

void GurobiBackend::InitBackend() {
  GRB_CALL( GRBemptyenv(&env) );

  GRB_CALL( GRBsetstrparam(env, "LogFile", "gurobi.log") );

  GRB_CALL( GRBstartenv(env) );

  /* Create an empty model */
  GRB_CALL( GRBnewmodel(env, &model, "amplgurobimodel", 0, NULL, NULL, NULL, NULL, NULL) );

}

void GurobiBackend::CloseBackend() {
  /* Free model */
  GRBfreemodel(model);

  /* Free environment */
  GRBfreeenv(env);
}

bool GurobiBackend::IsMIP() const {
  int isMIP;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_IS_MIP, &isMIP) );
  return isMIP;
}

bool GurobiBackend::IsQCP() const {
  int isQCP;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_IS_QCP, &isQCP) );
  return isQCP;
}

int GurobiBackend::NumberOfConstraints() const {
  int nc;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &nc) );
  return nc;
}

int GurobiBackend::NumberOfVariables() const {
  int nv;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &nv) );
  return nv;
}

int GurobiBackend::NumberOfObjectives() const {
  int no;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_NUMOBJ, &no) );
  return no;
}

void GurobiBackend::PrimalSolution(std::vector<double> &x) {
  int num_vars = NumberOfVariables();
  x.resize(num_vars);
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, num_vars, x.data());
}

void GurobiBackend::DualSolution(std::vector<double> &pi) {
  int num_cons = NumberOfConstraints();
  pi.resize(num_cons);
  GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, num_cons, pi.data());
}

double GurobiBackend::ObjectiveValue() const {
  double objval;
  GRB_CALL( GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval) );
  return objval;
}

double GurobiBackend::NodeCount() const {
  double ndcount;
  GRB_CALL( GRBgetdblattr(model, GRB_DBL_ATTR_NODECOUNT, &ndcount) );
  return ndcount;
}

double GurobiBackend::Niterations() const {
  double ni;
  GRB_CALL( GRBgetdblattr(model, GRB_DBL_ATTR_ITERCOUNT, &ni) );
  return ni;
}

std::string GurobiBackend::GetOptimizer(const SolverOption &) const {
  switch (optimizer_) {
  default:
    assert(false);
    // Fall through.
  case AUTO:  return "auto";
  case CP:    return "cp";
  case CPLEX: return "cplex";
  }
}

void GurobiBackend::SetOptimizer(const SolverOption &opt, fmt::StringRef value) {
  if (value == "auto")
    optimizer_ = AUTO;
  else if (value == "cp")
    optimizer_ = CP;
  else if (value == "cplex")
    optimizer_ = CPLEX;
  else
    throw InvalidOptionValue(opt, value);
}

void GurobiBackend::SetBoolOption(
    const SolverOption &opt, int value, Option id) {
  if (value != 0 && value != 1)
    throw InvalidOptionValue(opt, value);
  options_[id] = value;
}

void GurobiBackend::DoSetIntOption(
    const SolverOption &opt, int value, Option id) {
  if (value < 0)
    throw InvalidOptionValue(opt, value);
  options_[id] = value;
}


void GurobiBackend::SolveWithGurobi(
    Problem &p,
    Stats &stats, SolutionHandler &sh) {
  interrupter()->SetHandler(InterruptGurobi, model);

//  GRBwrite(model, "model.lp");

  stats.setup_time = GetTimeAndReset(stats.time);
  GRB_CALL( GRBoptimize(model) );
  stats.solution_time = GetTimeAndReset(stats.time);

  // Convert solution status.
  int solve_code = 0;
  std::string status =
      ConvertSolutionStatus(model, *interrupter(), solve_code);

  fmt::MemoryWriter writer;
  writer.write("{}: {}\n", long_name(), status);
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  vector<double> solution, dual_solution;
  if (solve_code < sol::INFEASIBLE) {
    PrimalSolution(solution);

    if (IsMIP()) {
      writer << NodeCount() << " nodes, ";
    } else {                                    // Also for QCP
      DualSolution(dual_solution);
    }
    writer << Niterations() << " iterations";

    if (NumberOfObjectives() > 0) {
      writer.write(", objective {}", FormatObjValue(ObjectiveValue()));
    }
  }
  sh.HandleSolution(solve_code, writer.c_str(),
      solution.empty() ? 0 : solution.data(),
      dual_solution.empty() ? 0 : dual_solution.data(), obj_value);
}

void GurobiBackend::Solve(Problem &p, SolutionHandler &sh) {
  Resolve(p, sh);
}

void GurobiBackend::InitProblemModificationPhase(const Problem &p) {
  stats.time = steady_clock::now();

  optimizer = optimizer_;
  if (optimizer == AUTO) {
    if (p.num_logical_cons() != 0 || p.has_nonlinear_cons() ||
        HasNonlinearObj(p)) {
      optimizer = CP;
    } else {
      optimizer = CPLEX;
    }
  }

}

void GurobiBackend::AddVariables(int n, double *lbs, double *ubs, var::Type *types) {
  std::vector<char> vtypes(n, GRB_CONTINUOUS);
  for (int var = 0; var < n; ++var) {
    if (types[var]!=var::Type::CONTINUOUS)
      vtypes[var] = GRB_INTEGER;
  }
  GRB_CALL( GRBaddvars(model, n, 0,
                       &n, &n, lbs, lbs,                  // placeholders, no matrix here
                       lbs, ubs, vtypes.data(), NULL) );
}
void GurobiBackend::AddLinearObjective( obj::Type sense, int nnz,
                         const double* c, const int* v) {
  if (1>=NumberOfObjectives()) {
    GRB_CALL( GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE,
                          obj::Type::MAX==sense ? GRB_MAXIMIZE : GRB_MINIMIZE) );
    for (int i = 0; i < nnz; ++i) {
      GRB_CALL( GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, v[i], c[i]) );
    }
  } else {
//    TODO
//    GRB_CALL( GRBsetobjectiven(model, 0, 1, 0.0, 0.0, 0.0, "primary",
//                               0.0, nnz, (int*)v, (double*)c) );
  }
}
void GurobiBackend::AddLinearConstraint(int nnz, const double* c, const int* v,
                         double lb, double ub) {
  /// TODO Separate abstraction
  if (lb==ub)
    GRB_CALL( GRBaddconstr(model, nnz, (int*)v, (double*)c, GRB_EQUAL, lb, NULL) );
  else {            // Let solver deal with lb>~ub etc.
    if (lb>MinusInfinity()) {
      GRB_CALL( GRBaddconstr(model, nnz, (int*)v, (double*)c, GRB_GREATER_EQUAL, lb, NULL) );
    }
    if (ub<Infinity()) {
      GRB_CALL( GRBaddconstr(model, nnz, (int*)v, (double*)c, GRB_LESS_EQUAL, ub, NULL) );
    }
  }
}
void GurobiBackend::FinishProblemModificationPhase() {
}

void GurobiBackend::Convert(Problem &p) {
  InitProblemModificationPhase(p);
}

void GurobiBackend::Resolve(Problem &p, SolutionHandler &sh) {

  SolveWithGurobi(p, stats, sh);
  double output_time = GetTimeAndReset(stats.time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          stats.setup_time, stats.solution_time, output_time);
  }
}


SolverPtr create_ilogcp(const char *) { return SolverPtr(new GurobiBackend()); }
}

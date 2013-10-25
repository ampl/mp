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

#include "solvers/ilogcp/ilogcp.h"

#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>

#include "solvers/ilogcp/concert.h"
#include "solvers/ilogcp/ilogcp_date.h"

using std::strcmp;
using std::vector;

#ifndef ILOGCP_NO_VERS
static char xxxvers[] = "ilogcp_options\0\n"
  "AMPL/IBM ILOG CP Optimizer Driver Version " qYYYYMMDD "\n";
#endif

namespace {

const char *InferenceLevels[] = {
  "default",
  "low",
  "basic",
  "medium",
  "extended",
  0
};

const char *Flags[] = {
  "off",
  "on",
  0
};

const char *SearchTypes[] = {
  "depthfirst",
  "restart",
  "multipoint",
  0
};

const char *Verbosities[] = {
  "quiet",
  "terse",
  "normal",
  "verbose",
  0
};

const char *TimeModes[] = {
  "cputime",
  "elapsedtime",
  0
};

ampl::OptionError GetOptionValueError(
    fmt::StringRef name, fmt::StringRef message) {
  throw ampl::OptionError(fmt::Format(
      "Can't get value of option {}: {}") << name.c_str() << message.c_str());
}

// An integer option.
class IntOption : public ampl::TypedSolverOption<int> {
 private:
  IloCP cp_;
  IloCP::IntParam param_;

 public:
  IntOption(const char *name, const char *description,
      IloCP cp, IloCP::IntParam p)
  : TypedSolverOption<int>(name, description), cp_(cp), param_(p) {}

  int GetValue() const;
  void SetValue(int value);
};

int IntOption::GetValue() const {
  try {
    return static_cast<int>(cp_.getParameter(param_));
  } catch (const IloException &e) {
    throw GetOptionValueError(name(), e.getMessage());
  }
}

void IntOption::SetValue(int value) {
  try {
    cp_.setParameter(param_, value);
  } catch (const IloException &) {
    throw ampl::InvalidOptionValue(name(), value);
  }
}

// An enumerated option.
class EnumOption : public ampl::TypedSolverOption<std::string> {
 private:
  IloCP cp_;
  IloCP::IntParam param_;
  int start_;           // start value for the enumerated options
  const char **values_; // string values for enum options
  bool accepts_auto_;   // true if the option accepts IloCP::Auto value

 public:
  EnumOption(const char *name, const char *description,
      IloCP cp, IloCP::IntParam p, int start,
      const char **values, bool accepts_auto = false)
  : ampl::TypedSolverOption<std::string>(name, description),
    cp_(cp), param_(p), start_(start), values_(values),
    accepts_auto_(accepts_auto) {
  }

  std::string GetValue() const;
  void SetValue(const char *value);
};

std::string EnumOption::GetValue() const {
  IloInt value = 0;
  try {
    value = cp_.getParameter(param_);
  } catch (const IloException &e) {
    throw GetOptionValueError(name(), e.getMessage());
  }
  if (value == IloCP::Auto && accepts_auto_)
    return "auto";
  if (values_) {
    for (int i = 0; values_[i]; ++i) {
      if (i + start_ == value)
        return values_[i];
    }
  }
  return str(fmt::Format("{}") << value);
}

void EnumOption::SetValue(const char *value) {
  try {
    char *end = 0;
    long intval = std::strtol(value, &end, 0);
    if (!*end) {
      if (intval != -1 || !accepts_auto_)
        intval += start_;
      cp_.setParameter(param_, intval);
      return;
    }
    if (values_) {
      // Search for a value in the list of known values.
      // Use linear search since the number of values is small.
      for (int i = 0; values_[i]; ++i) {
        if (strcmp(value, values_[i]) == 0) {
          cp_.setParameter(param_, i + start_);
          return;
        }
      }
    }
    if (accepts_auto_ && strcmp(value, "auto") == 0) {
      cp_.setParameter(param_, IloCP::Auto);
      return;
    }
  } catch (const IloException &) {}
  throw ampl::InvalidOptionValue(name(), value);
}

void GetSolution(IloCP cp, IloNumVarArray vars, vector<double> &solution) {
  for (int j = 0, n = vars.getSize(); j < n; ++j) {
    IloNumVar &v = vars[j];
    solution[j] = cp.isExtracted(v) ? cp.getValue(v) : v.getLB();
  }
}

bool HasNonlinearObj(const ampl::Problem &p) {
  if (p.num_objs() == 0)
    return false;
  ampl::NumericExpr expr = p.nonlinear_obj_expr(0);
  return expr && !ampl::Cast<ampl::NumericConstant>(expr);
}

std::string ConvertSolutionStatus(
    IloAlgorithm alg, const ampl::SignalHandler &sh, int &solve_code) {
  switch (alg.getStatus()) {
  default:
    // Fall through.
  case IloAlgorithm::Unknown:
    if (sh.stop()) {
      solve_code = 600;
      return "interrupted";
    }
    solve_code = 501;
    return "unknown solution status";
  case IloAlgorithm::Feasible:
    if (sh.stop()) {
      solve_code = 600;
      return "interrupted";
    }
    solve_code = 100;
    return "feasible solution";
  case IloAlgorithm::Optimal:
    solve_code = 0;
    return "optimal solution";
  case IloAlgorithm::Infeasible:
    solve_code = 200;
    return "infeasible problem";
  case IloAlgorithm::Unbounded:
    solve_code = 300;
    return "unbounded problem";
  case IloAlgorithm::InfeasibleOrUnbounded:
    solve_code = 201;
    return "infeasible or unbounded problem";
  case IloAlgorithm::Error:
    solve_code = 500;
    return "error";
  }
}
}

namespace ampl {

std::string IlogCPSolver::GetOptionHeader() {
  return "IlogCP Directives for AMPL\n"
      "--------------------------\n"
      "\n"
      "To set these directives, assign a string specifying their values to the AMPL "
      "option ilogcp_options.  For example:\n"
      "\n"
      "  ampl: option ilogcp_options 'optimalitytolerance=1e-6 searchtype=restart';\n"
      "\n"
      "Where both a number and a keyword are given, either may be used to specify "
      "the option setting.\n";
}

IlogCPSolver::IlogCPSolver() :
   Solver("ilogcp", 0, YYYYMMDD, MULTIPLE_SOL), cp_(env_), cplex_(env_),
   optimizer_(AUTO) {
  cp_.setIntParameter(IloCP::LogVerbosity, IloCP::Quiet);
  cplex_.setParam(IloCplex::MIPDisplay, 0);

  options_[DEBUGEXPR] = false;
  options_[USENUMBEROF] = true;
  options_[SOLUTION_LIMIT] = -1;

  set_long_name(fmt::Format("ilogcp {}.{}.{}")
      << IloConcertVersion::_ILO_MAJOR_VERSION
      << IloConcertVersion::_ILO_MINOR_VERSION
      << IloConcertVersion::_ILO_TECH_VERSION);
  set_version(fmt::Format("AMPL/IBM ILOG CP Optimizer [{} {}.{}.{}]")
      << IloConcertVersion::_ILO_NAME << IloConcertVersion::_ILO_MAJOR_VERSION
      << IloConcertVersion::_ILO_MINOR_VERSION
      << IloConcertVersion::_ILO_TECH_VERSION);

  AddSuffix("priority", 0, ASL_Sufkind_var);

  AddStrOption("optimizer",
      "Specifies which optimizer to use.  Possible values:\n"
      "      auto  = CP Optimizer if the problem has\n"
      "              nonlinear objective/constraints\n"
      "              or logical constraints, CPLEX\n"
      "              otherwise (default)\n"
      "      cp    = CP Optimizer\n"
      "      cplex = CPLEX Optimizer\n",
      &IlogCPSolver::GetOptimizer, &IlogCPSolver::SetOptimizer);

  // The following options are not implemented because corresponding
  // constraints are never generated by the driver:
  // - CountInferenceLevel
  // - SequenceInferenceLevel
  // - AllMinDistanceInferenceLevel
  // - ElementInferenceLevel
  // - PrecedenceInferenceLevel
  // - IntervalSequenceInferenceLevel
  // - NoOverlapInferenceLevel
  // - CumulFunctionInferenceLevel
  // - StateFunctionInferenceLevel

  AddOption(OptionPtr(new EnumOption("alldiffinferencelevel",
      "Inference level for 'alldiff' constraints.  Possible values:\n"
      "      0 = default\n"
      "      1 = low\n"
      "      2 = basic\n"
      "      3 = medium\n"
      "      4 = extended\n",
      cp_, IloCP::AllDiffInferenceLevel, IloCP::Default, InferenceLevels)));

  AddOption(OptionPtr(new IntOption("branchlimit",
      "Limit on the number of branches made before "
      "terminating a search.  Default = no limit.",
      cp_, IloCP::BranchLimit)));

  AddOption(OptionPtr(new IntOption("choicepointlimit",
      "Limit on the number of choice points created"
      "before terminating a search.  Default = no limit.",
      cp_, IloCP::ChoicePointLimit)));

  AddOption(OptionPtr(new EnumOption("constraintaggregation",
      "0 or 1 (default 1):  Whether to aggregate basic constraints.",
      cp_, IloCP::ConstraintAggregation, IloCP::Off, Flags)));

  AddIntOption("debugexpr",
      "0 or 1 (default 0):  Whether to print debugging "
      "information for expression trees.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::SetBoolOption, DEBUGEXPR);

  AddOption(OptionPtr(new EnumOption("defaultinferencelevel",
      "Default inference level for constraints.  Possible values:\n"
      "      1 = low\n"
      "      2 = basic\n"
      "      3 = medium\n"
      "      4 = extended\n",
      cp_, IloCP::DefaultInferenceLevel, IloCP::Default, InferenceLevels)));

  AddOption(OptionPtr(new EnumOption("distributeinferencelevel",
      "Inference level for 'distribute' constraints.  Possible values:\n"
      "      0 = default\n"
      "      1 = low\n"
      "      2 = basic\n"
      "      3 = medium\n"
      "      4 = extended\n",
      cp_, IloCP::DistributeInferenceLevel, IloCP::Default, InferenceLevels)));

  AddOption(OptionPtr(new EnumOption("dynamicprobing",
      "Use probing during search.  Possible values:\n"
      "     -1 = auto (default)\n"
      "      0 = off\n"
      "      1 = on\n",
      cp_, IloCP::DynamicProbing, IloCP::Off, Flags, true)));

  AddDblOption("dynamicprobingstrength",
      "Effort dedicated to dynamic probing as a factor "
      "of the total search effort.  Default = 0.03.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::DynamicProbingStrength);

  AddOption(OptionPtr(new IntOption("faillimit",
      "Limit on the number of failures allowed before terminating a search.  "
      "Default = no limit",
      cp_, IloCP::FailLimit)));

  AddOption(OptionPtr(new IntOption("logperiod",
      "Specifies how often the information in the search log is displayed.",
      cp_, IloCP::LogPeriod)));

  AddOption(OptionPtr(new EnumOption("logverbosity",
      "Verbosity of the search log.  Possible values:\n"
      "      0 = quiet (default)\n"
      "      1 = terse\n"
      "      2 = normal\n"
      "      3 = verbose\n",
      cp_, IloCP::LogVerbosity, IloCP::Quiet, Verbosities)));

  AddOption(OptionPtr(new IntOption("multipointnumberofsearchpoints",
      "Number of solutions for the multi-point search "
      "algorithm.  Default = 30.",
      cp_, IloCP::MultiPointNumberOfSearchPoints)));

  AddDblOption("optimalitytolerance",
      "Absolute tolerance on the objective value.  Default = 0.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::OptimalityTolerance);

  AddOption(OptionPtr(new EnumOption("outlev",
      "Synonym for \"logverbosity\".",
      cp_, IloCP::LogVerbosity, IloCP::Quiet, Verbosities)));

  AddOption(OptionPtr(new EnumOption("propagationlog",
      "Level of propagation trace reporting.  Possible values:\n"
      "      0 = quiet (default)\n"
      "      1 = terse\n"
      "      2 = normal\n"
      "      3 = verbose\n",
      cp_, IloCP::PropagationLog, IloCP::Quiet, Verbosities)));

  AddOption(OptionPtr(new IntOption("randomseed",
      "Seed for the random number generator.  Default = 0.",
      cp_, IloCP::RandomSeed)));

  AddDblOption("relativeoptimalitytolerance",
      "Relative tolerance on the objective value.  Default = 1e-4.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::RelativeOptimalityTolerance);

  AddOption(OptionPtr(new IntOption("restartfaillimit",
      "Number of failures allowed before restarting  search.  Default = 100.",
      cp_, IloCP::RestartFailLimit)));

  AddDblOption("restartgrowthfactor",
      "Increase of the number of allowed failures "
      "before restarting search.  Default = 1.05.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::RestartGrowthFactor);

  AddOption(OptionPtr(new EnumOption("searchtype",
      "Type of search used for solving a problem.  Possible values:\n"
      "      0 = depthfirst\n"
      "      1 = restart (default)\n"
      "      2 = multipoint\n",
      cp_, IloCP::SearchType, IloCP::DepthFirst, SearchTypes, true)));

  AddIntOption("solutionlimit",
      "Limit on the number of feasible solutions found before terminating "
      "a search.  Leaving the solution limit unspecified will make the "
      "optimizer search for an optimal solution if there is an objective "
      "function or for a feasible solution otherwise.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::DoSetIntOption,
      IlogCPSolver::SOLUTION_LIMIT);

  AddOption(OptionPtr(new EnumOption("temporalrelaxation",
      "0 or 1 (default 1):  Whether to use temporal relaxation.",
      cp_, IloCP::TemporalRelaxation, IloCP::Off, Flags)));

  AddDblOption("timelimit",
      "Limit on the CPU time spent solving before "
      "terminating a search.  Default = no limit.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::TimeLimit);

  AddOption(OptionPtr(new EnumOption("timemode",
      "Specifies how the time is measured in CP Optimizer.  Possible values:\n"
      "      0 = cputime (default)\n"
      "      1 = elapsedtime\n",
      cp_, IloCP::TimeMode, IloCP::CPUTime, TimeModes)));

  AddIntOption("usenumberof",
      "0 or 1 (default 1):  Whether to consolidate 'numberof' expressions "
      "by use of IloDistribute constraints.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::SetBoolOption,
      IlogCPSolver::USENUMBEROF);

  AddOption(OptionPtr(new EnumOption("workers",
      "Number of workers to run in parallel to solve a problem.  "
      "In addition to numeric values this option accepts the value "
      "\"auto\" since CP Optimizer version 12.3.  Default = auto.",
      cp_, IloCP::Workers, 0, 0, true)));

  AddIntOption<IlogCPSolver, int>("mipdisplay",
      "Frequency of displaying branch-and-bound information "
      "(for optimizing integer variables):\n"
      "      0 (default) = never\n"
      "      1 = each integer feasible solution\n"
      "      2 = every \"mipinterval\" nodes\n"
      "      3 = every \"mipinterval\" nodes plus\n"
      "          information on LP relaxations\n"
      "          (as controlled by \"display\")\n"
      "      4 = same as 2, plus LP relaxation info.\n"
      "      5 = same as 2, plus LP subproblem info.\n",
      &IlogCPSolver::GetCPLEXIntOption, &IlogCPSolver::SetCPLEXIntOption,
      IloCplex::MIPDisplay);

  AddIntOption<IlogCPSolver, int>("mipinterval",
      "Frequency of node logging for mipdisplay 2 or 3. Default = 0.",
      &IlogCPSolver::GetCPLEXIntOption, &IlogCPSolver::SetCPLEXIntOption,
      IloCplex::MIPInterval);
}

IlogCPSolver::~IlogCPSolver() {
  env_.end();
}

std::string IlogCPSolver::GetOptimizer(const char *) const {
  switch (optimizer_) {
  default:
    assert(false);
    // Fall through.
  case AUTO:  return "auto";
  case CP:    return "cp";
  case CPLEX: return "cplex";
  }
}

void IlogCPSolver::SetOptimizer(const char *name, const char *value) {
  if (strcmp(value, "auto") == 0)
    optimizer_ = AUTO;
  else if (strcmp(value, "cp") == 0)
    optimizer_ = CP;
  else if (strcmp(value, "cplex") == 0)
    optimizer_ = CPLEX;
  else
    throw InvalidOptionValue(name, value);
}

void IlogCPSolver::SetBoolOption(const char *name, int value, Option opt) {
  if (value != 0 && value != 1)
    throw InvalidOptionValue(name, value);
  options_[opt] = value;
}

void IlogCPSolver::DoSetIntOption(const char *name, int value, Option opt) {
  if (value < 0)
    throw InvalidOptionValue(name, value);
  options_[opt] = value;
}

double IlogCPSolver::GetCPDblOption(
    const char *name, IloCP::NumParam param) const {
  try {
    return cp_.getParameter(param);
  } catch (const IloException &e) {
    throw GetOptionValueError(name, e.getMessage());
  }
}

void IlogCPSolver::SetCPDblOption(
    const char *name, double value, IloCP::NumParam param) {
  try {
    cp_.setParameter(param, value);
  } catch (const IloException &) {
    throw InvalidOptionValue(name, value);
  }
}

int IlogCPSolver::GetCPLEXIntOption(const char *name, int param) const {
  // Use CPXgetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  int value = 0;
  int result = CPXgetintparam(cplex_.getImpl()->getCplexEnv(), param, &value);
  if (result != 0)
    throw GetOptionValueError(name, fmt::Format("CPLEX error = {}") << result);
  return value;
}

void IlogCPSolver::SetCPLEXIntOption(const char *name, int value, int param) {
  // Use CPXsetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  if (CPXsetintparam(cplex_.getImpl()->getCplexEnv(), param, value) != 0)
    throw InvalidOptionValue(name, value);
}

void IlogCPSolver::SolveWithCP(
    Problem &p, const NLToConcertConverter &converter, Stats &stats) {
  IloNumVarArray vars = converter.vars();
  IloIntVarArray priority_vars(env_);
  Suffix priority_suffix = p.suffix("priority", ASL_Sufkind_var);
  if (priority_suffix && priority_suffix.has_values()) {
    for (int i = 0, n = p.num_vars(); i < n; ++i) {
      if (priority_suffix.int_value(i) > 0)
        priority_vars.add(vars[i]);
    }
  }

  class InteruptCP : public Interruptible {
  private:
    IloCP &cp_;
  public:
    explicit InteruptCP(IloCP &cp) : cp_(cp) {}
    void Interrupt() { cp_.abortSearch(); }
  } interrupt(cp_);
  SignalHandler sh(*this, &interrupt);

  int num_objs = p.num_objs();
  vector<double> solution(p.num_vars());
  std::set< vector<double> > solutions;
  std::string feasible_sol_message =
      str(fmt::Format("{}: feasible solution") << long_name());
  bool multiple_sols = need_multiple_solutions();
  int solution_limit = GetOption(SOLUTION_LIMIT);
  if (solution_limit == -1)
    solution_limit = p.num_objs() > 0 ? INT_MAX : 1;
  stats.setup_time = GetTimeAndReset(stats.time);
  if (priority_vars.getSize() != 0) {
    IloSearchPhase phase(env_, priority_vars);
    cp_.startNewSearch(phase);
  } else {
    cp_.startNewSearch();
  }
  int num_solutions = 0;
  int solve_code = -1;
  std::string status;
  while (cp_.next()) {
    if (solution_limit == INT_MAX && !multiple_sols)
      continue;
    IloAlgorithm::Status s = cp_.getStatus();
    if (s != IloAlgorithm::Feasible && s != IloAlgorithm::Optimal)
      continue;
    GetSolution(cp_, vars, solution);
    if (num_objs == 0 && !solutions.insert(solution).second)
      continue;
    if (multiple_sols) {
      double obj_value = num_objs > 0 ?
          cp_.getObjValue() : std::numeric_limits<double>::quiet_NaN();
      HandleFeasibleSolution(p, feasible_sol_message,
          solution.data(), 0, obj_value);
    }
    if (++num_solutions >= solution_limit) {
      if (p.num_objs() > 0) {
        solve_code = 403;
        status = "solution limit";
      }
      break;
    }
  }
  cp_.endSearch();
  stats.solution_time = GetTimeAndReset(stats.time);

  // Convert solution status.
  if (solve_code == -1) {
    status = ConvertSolutionStatus(cp_, sh, solve_code);
    if (p.num_objs() > 0) {
      if (cp_.getInfo(IloCP::FailStatus) == IloCP::SearchStoppedByLimit) {
        solve_code = 400;
        status = "limit";
      }
    } else if (solve_code == 100)
      solve_code = 0;
  }
  p.set_solve_code(solve_code);

  fmt::Writer writer;
  writer.Format("{}: {}\n") << long_name() << status;
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  writer.Format("{} choice points, {} fails")
      << cp_.getInfo(IloCP::NumberOfChoicePoints)
      << cp_.getInfo(IloCP::NumberOfFails);
  IloAlgorithm::Status cpstatus = cp_.getStatus();
  if (cpstatus == IloAlgorithm::Feasible || cpstatus == IloAlgorithm::Optimal) {
    GetSolution(cp_, vars, solution);
    if (num_objs > 0) {
      obj_value = cp_.getObjValue();
      writer.Format(", objective {}") << ObjPrec(obj_value);
    }
  } else {
    solution.clear();
  }
  HandleSolution(p, writer.c_str(),
      solution.empty() ? 0 : &solution[0], 0, obj_value);
}

void IlogCPSolver::SolveWithCPLEX(
    Problem &p, const NLToConcertConverter &converter, Stats &stats) {
  IloCplex::Aborter aborter(env_);
  struct InteruptCPLEX : public Interruptible {
  private:
    IloCplex::Aborter &aborter_;
  public:
    explicit InteruptCPLEX(IloCplex::Aborter &aborter) : aborter_(aborter) {}
    void Interrupt() { aborter_.abort(); }
  } interrupt(aborter);
  cplex_.use(aborter);
  SignalHandler sh(*this, &interrupt);

  stats.setup_time = GetTimeAndReset(stats.time);
  cplex_.solve();
  stats.solution_time = GetTimeAndReset(stats.time);

  // Convert solution status.
  int solve_code = 0;
  std::string status = ConvertSolutionStatus(cplex_, sh, solve_code);
  p.set_solve_code(solve_code);

  fmt::Writer writer;
  writer.Format("{}: {}\n") << long_name() << status;
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  vector<double> solution, dual_solution;
  if (solve_code < 200) {
    int num_vars = p.num_vars();
    solution.resize(num_vars);
    IloNumVarArray vars = converter.vars();
    for (int j = 0; j < num_vars; ++j) {
      IloNumVar &v = vars[j];
      solution[j] = cplex_.isExtracted(v) ? cplex_.getValue(v) : v.getLB();
    }

    if (cplex_.isMIP()) {
      writer << cplex_.getNnodes() << " nodes, ";
    } else {
      IloRangeArray cons = converter.cons();
      IloInt num_cons = cons.getSize();
      dual_solution.resize(num_cons);
      for (IloInt i = 0; i < num_cons; ++i)
        dual_solution[i] = cplex_.getDual(cons[i]);
    }
    writer << cplex_.getNiterations() << " iterations";

    if (p.num_objs() > 0) {
      obj_value = cplex_.getObjValue();
      writer.Format(", objective {}") << ObjPrec(obj_value);
    }
  }
  HandleSolution(p, writer.c_str(),
      solution.empty() ? 0 : solution.data(),
      dual_solution.empty() ? 0 : dual_solution.data(), obj_value);
}

void IlogCPSolver::DoSolve(Problem &p) {
  Stats stats = {steady_clock::now()};

  Optimizer optimizer = optimizer_;
  if (optimizer == AUTO) {
    if (p.num_logical_cons() != 0 || p.num_nonlinear_cons() != 0 ||
        HasNonlinearObj(p)) {
      if (p.num_continuous_vars() != 0)
        throw Error("CP Optimizer doesn't support continuous variables");
      optimizer = CP;
    } else {
      optimizer = CPLEX;
    }
  }

  unsigned flags = 0;
  if (GetOption(USENUMBEROF) != 0)
    flags |= NLToConcertConverter::USENUMBEROF;
  if (GetOption(DEBUGEXPR) != 0)
    flags |= NLToConcertConverter::DEBUG;
  if (optimizer == CP)
    flags |= NLToConcertConverter::MULTIOBJ;
  NLToConcertConverter converter(env_, flags);
  converter.Convert(p);

  try {
    IloAlgorithm &cp_alg = cp_;
    (optimizer == CP ? cp_alg : cplex_).extract(converter.model());
  } catch (IloAlgorithm::CannotExtractException &e) {
    const IloExtractableArray &extractables = e.getExtractables();
    if (extractables.getSize() == 0)
      throw;
    throw UnsupportedExprError::CreateFromExprString(
        str(fmt::Format("{}") << extractables[0]));
  }

  if (optimizer == CP)
    SolveWithCP(p, converter, stats);
  else
    SolveWithCPLEX(p, converter, stats);
  double output_time = GetTimeAndReset(stats.time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n")
            << stats.setup_time << stats.solution_time << output_time;
  }
}
}

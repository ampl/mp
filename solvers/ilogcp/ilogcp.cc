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
 March 2020:   Updated to Concert 12.10 (Christian Valente)

 Possible improvements: Some sort of variable preference mechanism.

 Reference: "Extending an Algebraic Modeling Language to
 Support Constraint Programming" by Robert Fourer and David M. Gay,
 INFORMS Journal on Computing, Fall 2002, vol. 14, no. 4, 322-344
 (http://joc.journal.informs.org/content/14/4/322).
 */

#include "ilogcp.h"

#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>

#include "concert.h"

using std::strcmp;
using std::vector;

#ifndef ILOGCP_NO_VERS
static char xxxvers[] = "ilogcp_options\0\n"
  "AMPL/IBM ILOG CP Optimizer Driver Version " qYYYYMMDD "\n";
#endif

namespace {

const mp::OptionValueInfo INFERENCE_LEVELS[] = {
  // Default must be the first item.
  {"default",  0, IloCP::Default },
  {"low",      0, IloCP::Low     },
  {"basic",    0, IloCP::Basic   },
  {"medium",   0, IloCP::Medium  },
  {"extended", 0, IloCP::Extended}
};

const mp::OptionValueInfo FLAGS[] = {
  // Auto must be the first item.
  {"auto", 0, IloCP::Auto},
  {"off",  0, IloCP::Off },
  {"on",   0, IloCP::On  }
};

const mp::OptionValueInfo SEARCH_TYPES[] = {
  {"auto",       0, IloCP::Auto      },
  {"depthfirst", 0, IloCP::DepthFirst},
  {"restart",    0, IloCP::Restart   },
  {"multipoint", 0, IloCP::MultiPoint}
  #if CPX_VERSION >= 12090000
  ,{"iterativediving", 0, IloCP::IterativeDiving}
  #endif
  #if CPX_VERSION >= 22010000
  ,{"neighborhood", "(beta in CP Optimizer 22.1)", IloCP::Neighborhood}
  #endif
};

#if CPX_VERSION >= 12090000
const mp::OptionValueInfo KPI_DISPLAYS[] = {
  {"singleline",    "all KPIs on one line (default)", IloCP::SingleLine},
  {"multiplelines", "one line per KPI", IloCP::MultipleLines}
};
#endif

const mp::OptionValueInfo VERBOSITIES[] = {
  {"quiet",   0, IloCP::Quiet  },
  {"terse",   0, IloCP::Terse  },
  {"normal",  0, IloCP::Normal },
  {"verbose", 0, IloCP::Verbose}
};

const mp::OptionValueInfo TIME_MODES[] = {
  {"cputime",     0, IloCP::CPUTime    },
  {"elapsedtime", 0, IloCP::ElapsedTime}
};

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

const mp::OptionValueInfo AUTO_VALUE[] = {
  {"auto", 0, IloCP::Auto}
};

mp::OptionError GetOptionValueError(
    const mp::SolverOption &opt, fmt::StringRef message) {
  throw mp::OptionError(fmt::format(
      "Can't get value of option {}: {}", opt.name(), message));
}

// An integer option (64-bit, as CP Optimizer's IloInt parameters).
class IntOption : public mp::TypedSolverOption<fmt::LongLong> {
 private:
  IloCP cp_;
  IloCP::IntParam param_;

 public:
  IntOption(const char *name, const char *description,
      IloCP cp, IloCP::IntParam p)
  : mp::TypedSolverOption<fmt::LongLong>(name, description),
    cp_(cp), param_(p) {}

  void GetValue(fmt::LongLong &value) const;
  void SetValue(fmt::LongLong value);
};

void IntOption::GetValue(fmt::LongLong &value) const {
  try {
    value = cp_.getParameter(param_);   // IloInt is 64-bit: no truncation
  } catch (const IloException &e) {
    throw GetOptionValueError(*this, e.getMessage());
  }
}

void IntOption::SetValue(fmt::LongLong value) {
  try {
    cp_.setParameter(param_, static_cast<IloInt>(value));
  } catch (const IloException &) {
    throw mp::InvalidOptionValue(name(), value);
  }
}

// An enumerated option.
class EnumOption : public mp::TypedSolverOption<std::string> {
 private:
  IloCP cp_;
  IloCP::IntParam param_;
  int start_;           // start value for the enumerated options
  bool accepts_auto_;   // true if the option accepts IloCP::Auto value

 public:
  EnumOption(const char *name, const char *description,
      IloCP cp, IloCP::IntParam p, int start,
      mp::ValueArrayRef values, bool accepts_auto = false)
  : mp::TypedSolverOption<std::string>(name, description, values),
    cp_(cp), param_(p), start_(start), accepts_auto_(accepts_auto) {
  }

  void GetValue(std::string &value) const;
  void SetValue(fmt::StringRef value);
};

void EnumOption::GetValue(std::string &value) const {
  IloInt int_value = 0;
  try {
    int_value = cp_.getParameter(param_);
  } catch (const IloException &e) {
    throw GetOptionValueError(*this, e.getMessage());
  }
  for (mp::ValueArrayRef::iterator
      i = values().begin(), e = values().end(); i != e; ++i) {
    if (std::get<intptr_t>(i->data) == (intptr_t)int_value) {
      value = i->value;
      return;
    }
  }
  value = fmt::format("{}", int_value);
}

void EnumOption::SetValue(fmt::StringRef value) {
  try {
    // Copy the value adding a terminating null for strtol.
    fmt::MemoryWriter writer;
    writer << value;
    const char *str = writer.c_str();
    char *end = 0;
    long intval = std::strtol(str, &end, 0);
    if (!*end) {
      if (intval != -1 || !accepts_auto_)
        intval += start_;
      // With a full value table, accept only the listed values
      // (the library may silently accept undocumented ones).
      if (values().size() > 1) {
        bool known = accepts_auto_ && intval == IloCP::Auto;
        for (mp::ValueArrayRef::iterator
            i = values().begin(), e = values().end(); !known && i != e; ++i)
          known = std::get<intptr_t>(i->data) == (intptr_t)intval;
        if (!known)
          throw mp::InvalidOptionValue(name(), value);
      }
      cp_.setParameter(param_, intval);
      return;
    }
    // Search for a value in the list of known values.
    // Use linear search since the number of values is small.
    for (mp::ValueArrayRef::iterator
        i = values().begin(), e = values().end(); i != e; ++i) {
      if (strcmp(str, i->value) == 0) {
        cp_.setParameter(param_, std::get<intptr_t>(i->data));
        return;
      }
    }
  } catch (const IloException &) {}
  throw mp::InvalidOptionValue(name(), value);
}

void GetSolution(IloCP cp, IloNumVarArray vars, vector<double> &solution) {
  for (IloInt j = 0, n = vars.getSize(); j < n; ++j) {
    IloNumVar &v = vars[j];
    solution[j] = cp.isExtracted(v) ? cp.getValue(v) : v.getLB();
  }
}

bool HasNonlinearObj(const mp::Problem &p) {
  if (p.num_objs() == 0)
    return false;
  mp::NumericExpr expr = p.obj(0).nonlinear_expr();
  return expr && !mp::Cast<mp::NumericConstant>(expr);
}

// Driver-specific solve result codes (see mp::sol for the standard ones).
enum IlogCPResultCode {
  LIMIT_FEAS_FAIL        = mp::sol::LIMIT_FEAS_NEW,        // 430
  LIMIT_FEAS_BRANCH      = mp::sol::LIMIT_FEAS_NEW + 1,    // 431
  LIMIT_FEAS_CHOICEPOINT = mp::sol::LIMIT_FEAS_NEW + 2,    // 432
  LIMIT_NO_FEAS_FAIL        = mp::sol::LIMIT_NO_FEAS_NEW,      // 490
  LIMIT_NO_FEAS_BRANCH      = mp::sol::LIMIT_NO_FEAS_NEW + 1,  // 491
  LIMIT_NO_FEAS_CHOICEPOINT = mp::sol::LIMIT_NO_FEAS_NEW + 2   // 492
};

struct LimitKind {
  int code_feas, code_nofeas;
  const char *name;
};

// Determine which CP Optimizer limit stopped the search.
LimitKind DetectCPLimit(IloCP cp) {
  namespace sol = mp::sol;
  if (cp.getInfo(IloCP::NumberOfFails) >= cp.getParameter(IloCP::FailLimit))
    return {LIMIT_FEAS_FAIL, LIMIT_NO_FEAS_FAIL, "fail limit"};
  if (cp.getInfo(IloCP::NumberOfBranches) >= cp.getParameter(IloCP::BranchLimit))
    return {LIMIT_FEAS_BRANCH, LIMIT_NO_FEAS_BRANCH, "branch limit"};
  if (cp.getInfo(IloCP::NumberOfChoicePoints)
      >= cp.getParameter(IloCP::ChoicePointLimit))
    return {LIMIT_FEAS_CHOICEPOINT, LIMIT_NO_FEAS_CHOICEPOINT,
            "choice point limit"};
  if (cp.getParameter(IloCP::TimeLimit) < IloInfinity)
    return {sol::LIMIT_FEAS_TIME, sol::LIMIT_NO_FEAS_TIME, "time limit"};
  return {sol::LIMIT_FEAS, sol::LIMIT_NO_FEAS, "limit"};
}

std::string LimitMessage(const char *name, bool has_sol) {
  return fmt::format("{}, {}feasible solution", name, has_sol ? "" : "no ");
}

// Convert the CP Optimizer status to an MP solve result code.
std::string ConvertCPStatus(IloCP cp, const mp::Interrupter &interrupter,
                            bool has_obj, bool has_sol, int &solve_code) {
  namespace sol = mp::sol;
  switch (cp.getStatus()) {
  case IloAlgorithm::Optimal:
    solve_code = sol::SOLVED;
    return "optimal solution";
  case IloAlgorithm::Feasible:
    if (!has_obj) {   // satisfaction problem: a solution is the answer
      solve_code = sol::SOLVED;
      return "feasible solution";
    }
    break;
  case IloAlgorithm::Infeasible:
    solve_code = sol::INFEASIBLE;
    return "infeasible problem";
  case IloAlgorithm::Unbounded:
    solve_code = has_sol ? sol::UNBOUNDED_FEAS : sol::UNBOUNDED_NO_FEAS;
    return "unbounded problem";
  case IloAlgorithm::InfeasibleOrUnbounded:
    solve_code = sol::LIMIT_INF_UNB;
    return "infeasible or unbounded problem";
  case IloAlgorithm::Error:
    solve_code = sol::FAILURE;
    return "error";
  default:
    break;
  }
  // Feasible with an objective (optimality not proven) or Unknown.
  IloInt fail_status = cp.getInfo(IloCP::FailStatus);
  if (interrupter.Stop() || fail_status == IloCP::SearchStoppedByAbort) {
    solve_code = has_sol ? sol::LIMIT_FEAS_INTERRUPT : sol::LIMIT_NO_FEAS_INTERRUPT;
    return LimitMessage("interrupted", has_sol);
  }
  if (fail_status == IloCP::SearchStoppedByLimit) {
    LimitKind lim = DetectCPLimit(cp);
    solve_code = has_sol ? lim.code_feas : lim.code_nofeas;
    return LimitMessage(lim.name, has_sol);
  }
  if (has_sol) {
    solve_code = sol::UNCERTAIN;
    return "feasible solution, optimality not proven";
  }
  solve_code = sol::FAILURE + 1;
  return "unknown solution status";
}

// Convert the CPLEX status to an MP solve result code.
std::string ConvertCPLEXStatus(IloCplex cplex, const mp::Interrupter &interrupter,
                               bool has_sol, int &solve_code) {
  namespace sol = mp::sol;
  auto limit = [&](int code_feas, int code_nofeas, const char *name) {
    solve_code = has_sol ? code_feas : code_nofeas;
    return LimitMessage(name, has_sol);
  };
  switch (cplex.getCplexStatus()) {
  case IloCplex::Optimal:
  case IloCplex::OptimalTol:
  case IloCplex::OptimalPopulated:
  case IloCplex::OptimalPopulatedTol:
  case IloCplex::MultiObjOptimal:
    solve_code = sol::SOLVED;
    return "optimal solution";
  case IloCplex::OptimalInfeas:
    solve_code = sol::UNCERTAIN + 10;
    return "optimal with unscaled infeasibilities";
  case IloCplex::FirstOrder:
    solve_code = sol::OPTIMAL_LOCALLY;
    return "first-order optimal solution";
  case IloCplex::NumBest:
    solve_code = sol::UNCERTAIN;
    return "feasible solution, numerical difficulties";
  case IloCplex::Feasible:
  case IloCplex::MultiObjNonOptimal:
    solve_code = sol::UNCERTAIN;
    return "feasible solution, optimality not proven";
  case IloCplex::Infeasible:
  case IloCplex::MultiObjInfeasible:
    solve_code = sol::INFEASIBLE;
    return "infeasible problem";
  case IloCplex::Unbounded:
  case IloCplex::MultiObjUnbounded:
  case IloCplex::RelaxationUnbounded:
    solve_code = has_sol ? sol::UNBOUNDED_FEAS : sol::UNBOUNDED_NO_FEAS;
    return "unbounded problem";
  case IloCplex::InfOrUnbd:
  case IloCplex::MultiObjInfOrUnbd:
    solve_code = sol::LIMIT_INF_UNB;
    return "infeasible or unbounded problem";
  case IloCplex::AbortTimeLim:
  case IloCplex::AbortDetTimeLim:
    return limit(sol::LIMIT_FEAS_TIME, sol::LIMIT_NO_FEAS_TIME, "time limit");
  case IloCplex::AbortItLim:
    return limit(sol::LIMIT_FEAS_ITER, sol::LIMIT_NO_FEAS_ITER, "iteration limit");
  case IloCplex::NodeLimFeas:
  case IloCplex::NodeLimInfeas:
    return limit(sol::LIMIT_FEAS_NODES, sol::LIMIT_NO_FEAS_NODES, "node limit");
  case IloCplex::SolLim:
  case IloCplex::PopulateSolLim:
    solve_code = sol::LIMIT_FEAS_NUMSOLS;
    return "solution limit";
  case IloCplex::AbortObjLim:
  case IloCplex::AbortPrimObjLim:
  case IloCplex::AbortDualObjLim:
    return limit(sol::LIMIT_FEAS_BESTOBJ, sol::LIMIT_NO_FEAS_CUTOFF, "objective limit");
  case IloCplex::MemLimFeas:
  case IloCplex::MemLimInfeas:
  case IloCplex::FailFeasNoTree:
  case IloCplex::FailInfeasNoTree:
    return limit(sol::LIMIT_FEAS_SOFTMEM, sol::LIMIT_NO_FEAS_SOFTMEM, "memory limit");
  case IloCplex::AbortUser:
  case IloCplex::MultiObjStopped:
    if (interrupter.Stop())
      return limit(sol::LIMIT_FEAS_INTERRUPT, sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted");
    return limit(sol::LIMIT_FEAS, sol::LIMIT_NO_FEAS, "stopped");
  case IloCplex::FailFeas:
    solve_code = sol::LIMIT_FEAS_FAILURE;
    return "failure, feasible solution";
  case IloCplex::FailInfeas:
    solve_code = sol::FAILURE;
    return "failure";
  default:
    break;
  }
  if (interrupter.Stop())
    return limit(sol::LIMIT_FEAS_INTERRUPT, sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted");
  if (has_sol) {
    solve_code = sol::UNCERTAIN;
    return "feasible solution, unknown status";
  }
  solve_code = sol::FAILURE + 1;
  return "unknown solution status";
}

// IIS suffix values (same table as the MP backends).
enum IISValue {
  IIS_NON = 0, IIS_LOW = 1, IIS_FIX = 2, IIS_UPP = 3, IIS_MEM = 4, IIS_PMEM = 5
};

const char IIS_TABLE[] =
    "\n"
    "0\tnon\tnot in the iis\n"
    "1\tlow\tlower bound in the iis\n"
    "2\tfix\tboth bounds in the iis\n"
    "3\tupp\tupper bound in the iis\n"
    "4\tmem\tmember\n"
    "5\tpmem\tpossible member\n"
    "6\tplow\tpossibly lower bound\n"
    "7\tpupp\tpossibly upper bound\n";

template <typename Status>
int IISValueFromConflict(Status s, Status member, Status possible) {
  return s == member ? IIS_MEM : s == possible ? IIS_PMEM : IIS_NON;
}

bool InterruptCP(void *cp) {
  static_cast<IloCP*>(cp)->abortSearch();
  return true;
}

bool InterruptCPLEX(void *aborter) {
  static_cast<IloCplex::Aborter*>(aborter)->abort();
  return true;
}
}  // namespace

namespace mp {

#ifndef YYYYMMDD
  #define YYYYMMDD 20160218             // Update
#endif

IlogCPSolver::IlogCPSolver() :
   SolverImpl<Problem>("ilogcp", 0, YYYYMMDD, MULTIPLE_SOL | MULTIPLE_OBJ),
   cp_(env_), cplex_(env_), optimizer_(AUTO), optimizer(AUTO) {
  cp_.setIntParameter(IloCP::LogVerbosity, IloCP::Quiet);
  cplex_.setParam(IloCplex::MIPDisplay, 0);

  options_[DEBUGEXPR] = false;
  options_[USENUMBEROF] = true;
  options_[SOLUTION_LIMIT] = -1;
  options_[IIS_FIND] = 0;

  cp_defaults_.opttol = cp_.getParameter(IloCP::OptimalityTolerance);
  cp_defaults_.relopttol = cp_.getParameter(IloCP::RelativeOptimalityTolerance);
  cp_defaults_.seed = cp_.getParameter(IloCP::RandomSeed);

  AddSolveResults({
    { LIMIT_FEAS_FAIL, "fail limit, feasible solution" },
    { LIMIT_FEAS_BRANCH, "branch limit, feasible solution" },
    { LIMIT_FEAS_CHOICEPOINT, "choice point limit, feasible solution" },
    { LIMIT_NO_FEAS_FAIL, "fail limit, no feasible solution" },
    { LIMIT_NO_FEAS_BRANCH, "branch limit, no feasible solution" },
    { LIMIT_NO_FEAS_CHOICEPOINT, "choice point limit, no feasible solution" }
  });

  set_long_name(fmt::format("ilogcp {}.{}.{}",
      IloConcertVersion::_ILO_MAJOR_VERSION,
      IloConcertVersion::_ILO_MINOR_VERSION,
      IloConcertVersion::_ILO_TECH_VERSION));
  set_version(fmt::format("AMPL/IBM ILOG CP Optimizer [{} {}.{}.{}]",
      IloConcertVersion::_ILO_NAME, IloConcertVersion::_ILO_MAJOR_VERSION,
      IloConcertVersion::_ILO_MINOR_VERSION,
      IloConcertVersion::_ILO_TECH_VERSION));

  AddSuffix("priority", 0, suf::VAR);

  set_option_header(
      "IBM ILOG CPLEX CP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``ilogcp_options``. For example::\n"
      "\n"
      "  ampl: option ilogcp_options 'optimalitytolerance=1e-6 "
      "searchtype=restart';\n");

  AddStrOption("tech:optimizer optimizer",
      "Specifies which optimizer to use. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``auto``.",
      &IlogCPSolver::GetOptimizer, &IlogCPSolver::SetOptimizer, OPTIMIZERS);

  // CP options:

  // The following options are not implemented because corresponding
  // constraints are never generated by the driver:
  // - CountInferenceLevel, CountDifferentInferenceLevel
  // - SequenceInferenceLevel
  // - AllMinDistanceInferenceLevel
  // - PrecedenceInferenceLevel
  // - IntervalSequenceInferenceLevel
  // - NoOverlapInferenceLevel
  // - CumulFunctionInferenceLevel
  // - StateFunctionInferenceLevel

  AddOption(OptionPtr(new EnumOption("cp:alldiffinferencelevel alldiffinferencelevel",
      "Inference level for ``alldiff`` constraints. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``default``, which allows the inference "
      "strength of all ``alldiff`` constraints to be controlled via "
      "``cp:defaultinferencelevel``.",
      cp_, IloCP::AllDiffInferenceLevel, IloCP::Default, INFERENCE_LEVELS)));

  AddOption(OptionPtr(new IntOption("lim:branch branchlimit",
      "Limit on the number of branches made before "
      "terminating a search. Default = no limit.",
      cp_, IloCP::BranchLimit)));

  AddOption(OptionPtr(new IntOption("lim:choicepoint choicepointlimit",
      "Limit on the number of choice points created "
      "before terminating a search. Default = no limit.",
      cp_, IloCP::ChoicePointLimit)));

#if CPX_VERSION >= 12080000
  // Changed behaviour:
  // https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.8.0/ilog.odms.studio.help/CP_Optimizer/Release_notes/topics/relnotes_V1280_changes_api.html
  AddOption(OptionPtr(new EnumOption("pre:solve cppresolve presolve",
      "0 or 1 (default 1):  Whether to activate presolve on the CP optimizer.",
      cp_, IloCP::IntParam::Presolve, IloCP::Off,
      ValueArrayRef(FLAGS, 1))));
#else
  AddOption(OptionPtr(new EnumOption("pre:aggregate constraintaggregation",
      "0 or 1 (default 1):  Whether to aggregate basic constraints.",
      cp_, IloCP::ConstraintAggregation, IloCP::Off,
      ValueArrayRef(FLAGS, 1))));
#endif

  AddIntOption("tech:debugexpr debugexpr",
      "0 or 1 (default 0):  Whether to print debugging "
      "information for expression trees.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::SetBoolOption, DEBUGEXPR);

  AddOption(OptionPtr(new EnumOption("cp:defaultinferencelevel defaultinferencelevel",
      "Inference level for constraints that have inference level set to "
      "``default``. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``basic``.",
      cp_, IloCP::DefaultInferenceLevel, IloCP::Default,
      ValueArrayRef(INFERENCE_LEVELS, 1))));

  AddOption(OptionPtr(new EnumOption("cp:distributeinferencelevel distributeinferencelevel",
      "Inference level for aggregated ``numberof`` (``IloDistribute``) "
      "constraints. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``default``, which allows the inference "
      "strength of all aggregated ``numberof`` constraints to be controlled "
      "via ``cp:defaultinferencelevel``.",
      cp_, IloCP::DistributeInferenceLevel, IloCP::Default, INFERENCE_LEVELS)));

  AddOption(OptionPtr(new EnumOption("cp:dynamicprobing dynamicprobing",
      "Use probing during search. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``auto``.",
      cp_, IloCP::DynamicProbing, IloCP::Off, FLAGS, true)));

  AddDblOption("cp:dynamicprobingstrength dynamicprobingstrength",
      "Effort dedicated to dynamic probing as a factor "
      "of the total search effort. Default = 0.03; "
      "must be between 0.001 and 1000.",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::DynamicProbingStrength);

  AddOption(OptionPtr(new EnumOption("cp:elementinferencelevel elementinferencelevel",
      "Inference level for ``element`` (``IloElement``) constraints. "
      "Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``default``, which allows the inference "
      "strength of all ``element`` constraints to be controlled via "
      "``cp:defaultinferencelevel``.",
      cp_, IloCP::ElementInferenceLevel, IloCP::Default, INFERENCE_LEVELS)));

  AddOption(OptionPtr(new IntOption("lim:fail faillimit",
      "Limit on the number of failures allowed before terminating a search. "
      "Default = no limit.",
      cp_, IloCP::FailLimit)));

  AddOption(OptionPtr(new IntOption("tech:logperiod logperiod",
      "Specifies how often (in branches) the information in the search log "
      "is displayed. Default = 1000.",
      cp_, IloCP::LogPeriod)));

  AddOption(OptionPtr(new EnumOption("tech:outlev outlev logverbosity",
      "Verbosity of the search log. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``quiet``. With the CPLEX optimizer, values "
      "other than ``quiet`` turn on the CPLEX log (see tech:mipdisplay).",
      cp_, IloCP::LogVerbosity, IloCP::Quiet, VERBOSITIES)));

  AddOption(OptionPtr(new IntOption(
      "cp:multipointnumberofsearchpoints multipointnumberofsearchpoints",
      "Number of solutions for the multi-point search "
      "algorithm. Default = 30.",
      cp_, IloCP::MultiPointNumberOfSearchPoints)));

  AddDblOption("mip:gapabs optimalitytolerance absmipgap",
      "Absolute tolerance on the objective value. Default = 1e-9 "
      "(CP Optimizer); with the CPLEX optimizer, a non-default value "
      "is applied as the absolute MIP gap.",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::OptimalityTolerance);

#if CPX_VERSION < 12080000
  // Not available in 12.8, documented only in the Java docs
  // https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.8.0/ilog.odms.studio.help/CP_Optimizer/Release_notes/topics/relnotes_V1280_changes_java.html
  AddOption(OptionPtr(new EnumOption("tech:propagationlog propagationlog",
      "Level of propagation trace reporting. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``quiet``.",
      cp_, IloCP::PropagationLog, IloCP::Quiet, VERBOSITIES)));
#endif

  AddOption(OptionPtr(new IntOption("tech:seed randomseed seed",
      "Seed for the random number generator. Default = 0. "
      "With the CPLEX optimizer, a non-default value is applied as the "
      "CPLEX random seed.",
      cp_, IloCP::RandomSeed)));

  AddDblOption("mip:gap relativeoptimalitytolerance mipgap",
      "Relative tolerance on the objective value. Default = 1e-4; "
      "with the CPLEX optimizer, a non-default value is applied as the "
      "relative MIP gap.",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::RelativeOptimalityTolerance);

  AddOption(OptionPtr(new IntOption("cp:restartfaillimit restartfaillimit",
      "Number of failures allowed before restarting search "
      "(``cp:searchtype=restart``). Default = 100.",
      cp_, IloCP::RestartFailLimit)));

  AddDblOption("cp:restartgrowthfactor restartgrowthfactor",
      "Increase of the number of allowed failures "
      "before restarting search. Default = 1.15; must be >= 1.",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::RestartGrowthFactor);

  AddOption(OptionPtr(new EnumOption("cp:searchtype searchtype",
      "Type of search used for solving a problem. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``auto``.",
      cp_, IloCP::SearchType, IloCP::DepthFirst, SEARCH_TYPES, true)));

  AddIntOption("lim:sol solutionlimit sollimit",
      "Limit on the number of feasible solutions found before terminating "
      "a search. Leaving the solution limit unspecified will make the "
      "optimizer search for an optimal solution if there is an objective "
      "function or for a feasible solution otherwise. "
      "Also applied to the CPLEX optimizer.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::DoSetIntOption,
      IlogCPSolver::SOLUTION_LIMIT);

  AddOption(OptionPtr(new EnumOption("cp:temporalrelaxation temporalrelaxation",
      "Whether to use the internal temporal relaxation. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``auto``.",
      cp_, IloCP::TemporalRelaxation, IloCP::Off, FLAGS, true)));

  AddDblOption("lim:time timelimit timelim",
      "Limit on the time (measured as specified by tech:timemode, "
      "elapsed time by default) spent solving before terminating a search. "
      "Default = no limit. Also applied to the CPLEX optimizer.",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::TimeLimit);

  AddOption(OptionPtr(new EnumOption("tech:timemode timemode",
      "Specifies how the time is measured in CP Optimizer. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``elapsedtime``.",
      cp_, IloCP::TimeMode, IloCP::CPUTime, TIME_MODES)));

  AddIntOption("cvt:usenumberof usenumberof",
      "0 or 1 (default 1):  Whether to aggregate ``numberof`` expressions "
      "by use of ``IloDistribute`` constraints.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::SetBoolOption,
      IlogCPSolver::USENUMBEROF);

  AddOption(OptionPtr(new EnumOption("tech:threads workers threads",
      "Number of workers to run in parallel to solve a problem. "
      "In addition to numeric values this option accepts the value "
      "``auto``. Default = ``auto`` (one worker per CPU core). "
      "Also applied as the number of threads of the CPLEX optimizer.",
      cp_, IloCP::Workers, 0, AUTO_VALUE, true)));

#if CPX_VERSION >= 12060100
  AddDblOption("cp:fdsemphasis failuredirectedsearchemphasis",
      "Specifies the number of workers that use failure-directed search once "
      "it has started. The value does not have to be integer. For example, "
      "value 1.5 means that first worker spends 100% of the time by "
      "failure-directed search, second worker 50% and remaining workers 0%. "
      "Default = auto (-1; depends on actual performance of the "
      "failure-directed search).",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::FailureDirectedSearchEmphasis);

  AddStrOption("tech:dumpfile dumpfile",
      "Specifies the name of a file where to dump the model before "
      "solving it. This file name must have extension ``.cpo``. "
      "Default = \"\" (don't dump the model).",
      &IlogCPSolver::GetFile, &IlogCPSolver::SetFile, DUMP_FILE);

  AddStrOption("tech:writemodel exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. With CP Optimizer this file name must have extension "
      "``.cpo``; with the CPLEX optimizer the usual CPLEX formats "
      "(``.lp``, ``.mps``, ``.sav``) can be used. "
      "Default = \"\" (don't export the model).",
      &IlogCPSolver::GetFile, &IlogCPSolver::SetFile, EXPORT_FILE);
#endif

#if CPX_VERSION >= 12060200
  AddIntOption("tech:warninglevel warninglevel",
      "Specifies the highest warning level to be displayed, all warnings "
      "higher than this level are masked. CP Optimizer warning levels run "
      "from 1 to 3, so setting this option to 0 turns off all warnings. "
      "Warnings issued may indicate potential errors or inefficiencies in "
      "your model. Default = 2.",
      &IlogCPSolver::GetCPOption<int>, &IlogCPSolver::SetCPOption<int>,
      IloCP::WarningLevel);
#endif

#if CPX_VERSION >= 12090000
  AddOption(OptionPtr(new EnumOption("tech:kpidisplay kpidisplay",
      "How KPIs are displayed in the CP Optimizer log. Possible values:\n"
      "\n"
      ".. value-table::\n",
      cp_, IloCP::KPIDisplay, IloCP::SingleLine, KPI_DISPLAYS)));

  AddOption(OptionPtr(new EnumOption("tech:printmodeldetails printmodeldetailsinmessages",
      "0 or 1 (default 1):  Whether CP Optimizer prints relevant model "
      "fragments in error and warning messages.",
      cp_, IloCP::PrintModelDetailsInMessages, IloCP::Off,
      ValueArrayRef(FLAGS, 1))));

  AddOption(OptionPtr(new EnumOption("tech:modelanonymizer modelanonymizer",
      "0 or 1 (default 0):  Whether to anonymize names in models written "
      "by tech:dumpfile.",
      cp_, IloCP::ModelAnonymizer, IloCP::Off, ValueArrayRef(FLAGS, 1))));

  AddOption(OptionPtr(new EnumOption("tech:usefilelocations usefilelocations",
      "0 or 1 (default 1):  Whether source file locations are processed "
      "and exported by CP Optimizer.",
      cp_, IloCP::UseFileLocations, IloCP::Off, ValueArrayRef(FLAGS, 1))));

  AddOption(OptionPtr(new EnumOption("tech:logsearchtags logsearchtags",
      "0 or 1 (default 0):  Whether to display failure tags "
      "in the CP Optimizer log.",
      cp_, IloCP::LogSearchTags, IloCP::Off, ValueArrayRef(FLAGS, 1))));

  AddOption(OptionPtr(new EnumOption("cp:fds failuredirectedsearch",
      "0 or 1 (default 1):  Whether to use failure-directed search.",
      cp_, IloCP::FailureDirectedSearch, IloCP::Off, ValueArrayRef(FLAGS, 1))));

  AddOption(OptionPtr(new IntOption("cp:fdsmaxmemory failuredirectedsearchmaxmemory",
      "Maximum memory in bytes available to failure-directed search. "
      "Default = 104857600 (100 MB).",
      cp_, IloCP::FailureDirectedSearchMaxMemory)));

  // Conflict refiner (IIS):

  AddIntOption("iis:find iisfind iis",
      "0 or 1 (default 0):  Whether to run the conflict refiner when the "
      "problem is infeasible and return an IIS in suffix .iis "
      "(``mem`` = member, ``pmem`` = possible member, ``non`` = not in the "
      "IIS) on constraints and, with CP Optimizer, on variables. "
      "Members among logical constraints are counted in the solve message.",
      &IlogCPSolver::DoGetIntOption, &IlogCPSolver::SetBoolOption,
      IlogCPSolver::IIS_FIND);

  AddDblOption("lim:conflicttime conflictrefinertimelimit",
      "Time limit in seconds for the CP Optimizer conflict refiner. "
      "Default = no limit.",
      &IlogCPSolver::GetCPOption<double>, &IlogCPSolver::SetCPOption<double>,
      IloCP::ConflictRefinerTimeLimit);

  AddOption(OptionPtr(new IntOption("lim:conflictiter conflictrefineriterationlimit",
      "Limit on the number of iterations of the CP Optimizer conflict "
      "refiner. Default = no limit.",
      cp_, IloCP::ConflictRefinerIterationLimit)));

  AddOption(OptionPtr(new IntOption("lim:conflictbranch conflictrefinerbranchlimit",
      "Limit on the total number of branches of the CP Optimizer conflict "
      "refiner. Default = no limit.",
      cp_, IloCP::ConflictRefinerBranchLimit)));

  AddOption(OptionPtr(new IntOption("lim:conflictfail conflictrefinerfaillimit",
      "Limit on the total number of failures of the CP Optimizer conflict "
      "refiner. Default = no limit.",
      cp_, IloCP::ConflictRefinerFailLimit)));

  AddOption(OptionPtr(new EnumOption("iis:onvariables conflictrefineronvariables",
      "0 or 1 (default 0):  Whether the CP Optimizer conflict refiner "
      "also refines variable domains.",
      cp_, IloCP::ConflictRefinerOnVariables, IloCP::Off,
      ValueArrayRef(FLAGS, 1))));
#endif

  // CPLEX options:

  AddIntOption<IlogCPSolver, int>("tech:mipdisplay mipdisplay",
      "Frequency of displaying branch-and-bound information "
      "(for optimizing integer variables with the CPLEX optimizer):\n"
      "\n"
      "| 0 (default) - never\n"
      "| 1 - each integer feasible solution\n"
      "| 2 - every ``tech:mipinterval`` nodes\n"
      "| 3 - every ``tech:mipinterval`` nodes plus cutting planes info\n"
      "| 4 - same as 3, plus LP relaxation info\n"
      "| 5 - same as 4, plus LP subproblem info.\n",
      &IlogCPSolver::GetCPLEXIntOption, &IlogCPSolver::SetCPLEXIntOption,
      IloCplex::MIPDisplay);

  AddIntOption<IlogCPSolver, int>("tech:mipinterval mipinterval",
      "Frequency of node logging for tech:mipdisplay >= 2. Default = 0 "
      "(automatic).",
      &IlogCPSolver::GetCPLEXIntOption, &IlogCPSolver::SetCPLEXIntOption,
      IloCplex::MIPInterval);
}

IlogCPSolver::~IlogCPSolver() {
  env_.end();
}

std::string IlogCPSolver::GetOptimizer(const SolverOption &) const {
  switch (optimizer_) {
  default:
    assert(false);
    // Fall through.
  case AUTO:  return "auto";
  case CP:    return "cp";
  case CPLEX: return "cplex";
  }
}

void IlogCPSolver::SetOptimizer(const SolverOption &opt, fmt::StringRef value) {
  if (value == "auto")
    optimizer_ = AUTO;
  else if (value == "cp")
    optimizer_ = CP;
  else if (value == "cplex")
    optimizer_ = CPLEX;
  else
    throw InvalidOptionValue(opt, value);
}

void IlogCPSolver::SetBoolOption(
    const SolverOption &opt, int value, Option id) {
  if (value != 0 && value != 1)
    throw InvalidOptionValue(opt, value);
  options_[id] = value;
}

void IlogCPSolver::DoSetIntOption(
    const SolverOption &opt, int value, Option id) {
  if (value < 0)
    throw InvalidOptionValue(opt, value);
  options_[id] = value;
}

template <typename T>
T IlogCPSolver::GetCPOption(
    const SolverOption &opt, typename ParamTraits<T>::Type param) const {
  try {
    return cp_.getParameter(param);
  } catch (const IloException &e) {
    throw GetOptionValueError(opt, e.getMessage());
  }
}

template <typename T>
void IlogCPSolver::SetCPOption(
    const SolverOption &opt, T value, typename ParamTraits<T>::Type param) {
  try {
    cp_.setParameter(param, value);
  } catch (const IloException &) {
    throw InvalidOptionValue(opt, value);
  }
}

int IlogCPSolver::GetCPLEXIntOption(const SolverOption &opt, int param) const {
  // Use CPXgetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  int value = 0;
  int result = CPXgetintparam(cplex_.getImpl()->getCplexEnv(), param, &value);
  if (result != 0)
    throw GetOptionValueError(opt, fmt::format("CPLEX error = {}", result));
  return value;
}

void IlogCPSolver::SetCPLEXIntOption(
    const SolverOption &opt, int value, int param) {
  // Use CPXsetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  if (CPXsetintparam(cplex_.getImpl()->getCplexEnv(), param, value) != 0)
    throw InvalidOptionValue(opt, value);
}

void ReportVariableConflict(Problem &p, IloCP cp, IloNumVarArray vars,
                            fmt::MemoryWriter &writer);

void IlogCPSolver::SolveWithCP(
    Problem &p, const MPToConcertConverter &converter,
    Stats &stats, SolutionHandler &sh) {
	
#if CPX_VERSION >= 12060100
  if (!filenames_[DUMP_FILE].empty())
    cp_.dumpModel(filenames_[DUMP_FILE].c_str());
  if (!filenames_[EXPORT_FILE].empty())
    cp_.exportModel(filenames_[EXPORT_FILE].c_str());
#endif

  IloNumVarArray vars = converter.vars();
  IloIntVarArray priority_vars(env_);
  IntSuffix priority_suffix = p.suffixes(suf::VAR).Find<int>("priority");
  if (priority_suffix) {
    for (int i = 0, n = p.num_vars(); i < n; ++i) {
      if (priority_suffix.value(i) > 0)
        priority_vars.add(vars[i]);
    }
  }

  interrupter()->SetHandler(InterruptCP, &cp_);

  int num_objs = p.num_objs();
  vector<double> solution(p.num_vars());
  std::set< vector<double> > solutions;
  std::string feasible_sol_message =
      fmt::format("{}: feasible solution", long_name());
  bool multiple_sols = need_multiple_solutions();
  int solution_limit = GetOption(SOLUTION_LIMIT);
  if (solution_limit == -1)
    solution_limit = p.num_objs() > 0 ? INT_MAX : 1;
  stats.setup_time = GetTimeAndReset(stats.time);
  if (priority_vars.getSize() != 0) {
    #if CPX_VERSION >= 12080000
    IloGoal goal(env_, priority_vars);
    #else
    IloSearchPhase goal(env_, priority_vars);
    #endif  
    cp_.startNewSearch(goal);
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
      sh.HandleFeasibleSolution(feasible_sol_message,
          solution.data(), 0, obj_value);
    }
    if (++num_solutions >= solution_limit) {
      if (num_objs > 0) {
        solve_code = sol::LIMIT_FEAS_NUMSOLS;
        status = "solution limit";
      }
      break;
    }
  }
  cp_.endSearch();
  stats.solution_time = GetTimeAndReset(stats.time);

  IloAlgorithm::Status cpstatus = cp_.getStatus();
  bool has_sol = cpstatus == IloAlgorithm::Feasible
      || cpstatus == IloAlgorithm::Optimal;

  // Convert solution status.
  if (solve_code == -1)
    status = ConvertCPStatus(cp_, *interrupter(), num_objs > 0, has_sol,
                             solve_code);

  fmt::MemoryWriter writer;
  writer.write("{}: {}\n", long_name(), status);
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  writer.write("{} choice points, {} fails",
      cp_.getInfo(IloCP::NumberOfChoicePoints),
      cp_.getInfo(IloCP::NumberOfFails));
  if (has_sol) {
    GetSolution(cp_, vars, solution);
    if (num_objs > 0) {
      obj_value = cp_.getObjValue();
      writer.write(", objective {}", FormatObjValue(obj_value));
    }
  } else {
    solution.clear();
  }
  if (GetOption(IIS_FIND) != 0 &&
      (solve_code == sol::INFEASIBLE || solve_code == sol::LIMIT_INF_UNB)) {
    ReportConflict(p, cp_, converter, writer, false);
    ReportVariableConflict(p, cp_, vars, writer);
  }
  sh.HandleSolution(solve_code, writer.c_str(),
      solution.empty() ? 0 : &solution[0], 0, obj_value);
}

// Run the conflict refiner of the optimizer and report the result.
template <typename Alg>
void IlogCPSolver::ReportConflict(
    Problem &p, Alg alg, const MPToConcertConverter &converter,
    fmt::MemoryWriter &writer, bool with_var_bounds) {
  typedef typename Alg::ConflictStatus Status;
  IloRangeArray cons = converter.cons();
  IloConstraintArray lcons = converter.lcons();
  IloNumVarArray vars = converter.vars();
  IloConstraintArray all(env_);
  for (IloInt i = 0, n = cons.getSize(); i < n; ++i)
    all.add(cons[i]);
  for (IloInt i = 0, n = lcons.getSize(); i < n; ++i)
    all.add(lcons[i]);
  // Variable bounds as explicit constraints (used with CPLEX, whose
  // conflict refiner reports constraints only).
  int num_vars = p.num_vars();
  std::vector<IloInt> lb_index(num_vars, -1), ub_index(num_vars, -1);
  IloConstraintArray bound_cons(env_);
  IloModel model = converter.model();
  if (with_var_bounds) {
    for (int i = 0; i < num_vars; ++i) {
      if (!alg.isExtracted(vars[i])) continue;
      double lb = vars[i].getLB(), ub = vars[i].getUB();
      if (lb > -IloInfinity) {
        lb_index[i] = all.getSize();
        IloConstraint c = (vars[i] >= lb);
        bound_cons.add(c);
        all.add(c);
      }
      if (ub < IloInfinity) {
        ub_index[i] = all.getSize();
        IloConstraint c = (vars[i] <= ub);
        bound_cons.add(c);
        all.add(c);
      }
    }
    // The refiner only considers extracted constraints:
    // add the bound constraints to the model temporarily.
    model.add(bound_cons);
  }
  IloNumArray prefs(env_, all.getSize());
  for (IloInt i = 0, n = all.getSize(); i < n; ++i)
    prefs[i] = 1.0;
  bool found = false;
  try {
    found = alg.refineConflict(all, prefs) != IloFalse;
  } catch (const IloException &e) {
    writer.write("\nConflict refiner failed: {}", e.getMessage());
    found = false;
  }
  if (!found) {
    writer.write("\nConflict refiner did not find a conflict");
    if (with_var_bounds) model.remove(bound_cons);
    bound_cons.endElements();
    bound_cons.end();
    all.end();
    prefs.end();
    return;
  }
  // Suffix .iis on constraints covers algebraic then logical constraints.
  int num_con = p.num_algebraic_cons(), num_lcon = p.num_logical_cons();
  int num_lcon_mem = 0, num_con_mem = 0;
  std::vector<int> con_iis(num_con + num_lcon, IIS_NON);
  for (int i = 0; i < num_con && i < cons.getSize(); ++i) {
    con_iis[i] = IISValueFromConflict<Status>(alg.getConflict(cons[i]),
        Alg::ConflictMember, Alg::ConflictPossibleMember);
    if (con_iis[i] != IIS_NON) ++num_con_mem;
  }
  for (int i = 0; i < num_lcon && i < lcons.getSize(); ++i) {
    con_iis[num_con + i] = IISValueFromConflict<Status>(
        alg.getConflict(lcons[i]), Alg::ConflictMember, Alg::ConflictPossibleMember);
    if (con_iis[num_con + i] != IIS_NON) ++num_lcon_mem;
  }
  writer.write("\nIIS: {} algebraic and {} logical constraints in conflict",
               num_con_mem, num_lcon_mem);
  try {
    MutIntSuffix suf = p.suffixes(suf::CON).template Add<int>(
          "iis", suf::CON | suf::OUTPUT, num_con + num_lcon, IIS_TABLE);
    for (int i = 0; i < num_con + num_lcon; ++i)
      suf.set_value(i, con_iis[i]);
  } catch (const Error &e) {
    writer.write("\nCannot return suffix .iis on constraints: {}", e.what());
  }
  if (with_var_bounds) {
    int num_var_mem = 0;
    std::vector<int> var_iis(num_vars, IIS_NON);
    for (int i = 0; i < num_vars; ++i) {
      bool lb_in = lb_index[i] >= 0
          && alg.getConflict(all[lb_index[i]]) != Alg::ConflictExcluded;
      bool ub_in = ub_index[i] >= 0
          && alg.getConflict(all[ub_index[i]]) != Alg::ConflictExcluded;
      var_iis[i] = lb_in && ub_in ? IIS_FIX : lb_in ? IIS_LOW
                                           : ub_in ? IIS_UPP : IIS_NON;
      if (var_iis[i] != IIS_NON) ++num_var_mem;
    }
    writer.write(", {} variable bounds", num_var_mem);
    try {
      MutIntSuffix suf = p.suffixes(suf::VAR).template Add<int>(
            "iis", suf::VAR | suf::OUTPUT, num_vars, IIS_TABLE);
      for (int i = 0; i < num_vars; ++i)
        suf.set_value(i, var_iis[i]);
    } catch (const Error &e) {
      writer.write("\nCannot return suffix .iis on variables: {}", e.what());
    }
  }
  if (with_var_bounds) model.remove(bound_cons);
  bound_cons.endElements();
  bound_cons.end();
  all.end();
  prefs.end();
}

// Variable conflicts are available from CP Optimizer only.
void ReportVariableConflict(Problem &p, IloCP cp, IloNumVarArray vars,
                            fmt::MemoryWriter &writer) {
  int num_vars = p.num_vars(), num_mem = 0;
  std::vector<int> var_iis(num_vars, IIS_NON);
  for (int i = 0; i < num_vars; ++i) {
    if (!cp.isExtracted(vars[i])) continue;
    var_iis[i] = IISValueFromConflict<IloCP::ConflictStatus>(
        cp.getConflict(vars[i]), IloCP::ConflictMember, IloCP::ConflictPossibleMember);
    if (var_iis[i] != IIS_NON) ++num_mem;
  }
  writer.write(", {} variables", num_mem);
  try {
    MutIntSuffix suf = p.suffixes(suf::VAR).Add<int>(
          "iis", suf::VAR | suf::OUTPUT, num_vars, IIS_TABLE);
    for (int i = 0; i < num_vars; ++i)
      suf.set_value(i, var_iis[i]);
  } catch (const Error &e) {
    writer.write("\nCannot return suffix .iis on variables: {}", e.what());
  }
}

void IlogCPSolver::SolveWithCPLEX(
    Problem &p, const MPToConcertConverter &converter,
    Stats &stats, SolutionHandler &sh) {
  IloCplex::Aborter aborter(env_);
  cplex_.use(aborter);
  interrupter()->SetHandler(InterruptCPLEX, &aborter);

  if (!filenames_[EXPORT_FILE].empty())
    cplex_.exportModel(filenames_[EXPORT_FILE].c_str());

  ApplyCommonOptionsToCPLEX();

  stats.setup_time = GetTimeAndReset(stats.time);
  cplex_.solve();
  stats.solution_time = GetTimeAndReset(stats.time);

  // Convert solution status.
  IloAlgorithm::Status cpxstatus = cplex_.getStatus();
  bool has_sol = cpxstatus == IloAlgorithm::Feasible
      || cpxstatus == IloAlgorithm::Optimal;
  int solve_code = 0;
  std::string status =
      ConvertCPLEXStatus(cplex_, *interrupter(), has_sol, solve_code);

  fmt::MemoryWriter writer;
  writer.write("{}: {}\n", long_name(), status);
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  vector<double> solution, dual_solution;
  if (has_sol) {
    int num_vars = p.num_vars();
    solution.resize(num_vars);
    IloNumVarArray vars = converter.vars();
    for (int j = 0; j < num_vars; ++j) {
      IloNumVar &v = vars[j];
      solution[j] = cplex_.isExtracted(v) ? cplex_.getValue(v) : v.getLB();
    }

    if (cplex_.isMIP()) {
      writer << cplex_.getNnodes() << " nodes, ";
    } else if (!cplex_.isQC()) {
      IloRangeArray cons = converter.cons();
      IloInt num_cons = cons.getSize();
      dual_solution.resize(num_cons);
      for (IloInt i = 0; i < num_cons; ++i)
        dual_solution[i] = cplex_.getDual(cons[i]);
    }
    writer << cplex_.getNiterations() << " iterations";

    if (p.num_objs() > 0) {
      obj_value = cplex_.getObjValue();
      writer.write(", objective {}", FormatObjValue(obj_value));
    }
  }
  if (GetOption(IIS_FIND) != 0 &&
      (solve_code == sol::INFEASIBLE || solve_code == sol::LIMIT_INF_UNB))
    ReportConflict(p, cplex_, converter, writer, true);
  sh.HandleSolution(solve_code, writer.c_str(),
      solution.empty() ? 0 : solution.data(),
      dual_solution.empty() ? 0 : dual_solution.data(), obj_value);
}

void IlogCPSolver::ApplyCommonOptionsToCPLEX() {
  CPXENVptr env = cplex_.getImpl()->getCplexEnv();
  IloNum time_limit = cp_.getParameter(IloCP::TimeLimit);
  if (time_limit < IloInfinity)
    CPXsetdblparam(env, CPXPARAM_TimeLimit, time_limit);
  int solution_limit = GetOption(SOLUTION_LIMIT);
  if (solution_limit > 0)
    CPXsetintparam(env, CPXPARAM_MIP_Limits_Solutions, solution_limit);
  IloInt workers = cp_.getParameter(IloCP::Workers);
  if (workers > 0)
    CPXsetintparam(env, CPXPARAM_Threads, static_cast<int>(workers));
  IloInt seed = cp_.getParameter(IloCP::RandomSeed);
  if (seed != cp_defaults_.seed)
    CPXsetintparam(env, CPXPARAM_RandomSeed, static_cast<int>(seed));
  IloNum opttol = cp_.getParameter(IloCP::OptimalityTolerance);
  if (opttol != cp_defaults_.opttol)
    CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_AbsMIPGap, opttol);
  IloNum relopttol = cp_.getParameter(IloCP::RelativeOptimalityTolerance);
  if (relopttol != cp_defaults_.relopttol)
    CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, relopttol);
  if (cp_.getParameter(IloCP::LogVerbosity) != IloCP::Quiet) {
    int level = cp_.getParameter(IloCP::LogVerbosity) >= IloCP::Normal ? 2 : 1;
    int mip_display = 0;
    CPXgetintparam(env, CPXPARAM_MIP_Display, &mip_display);
    if (mip_display == 0)
      CPXsetintparam(env, CPXPARAM_MIP_Display, level);
    CPXsetintparam(env, CPXPARAM_Simplex_Display, level);
    CPXsetintparam(env, CPXPARAM_Barrier_Display, level);
  }
}

void IlogCPSolver::Solve(Problem &p, SolutionHandler &sh) {
  Convert(p);
  Resolve(p, sh);
}

void IlogCPSolver::InitProblemModificationPhase(const Problem &p) {
  stats.time = std::chrono::steady_clock::now();

  optimizer = optimizer_;
  if (optimizer == AUTO) {
    if (p.num_logical_cons() != 0 || p.has_nonlinear_cons() ||
        HasNonlinearObj(p)) {
      optimizer = CP;
    } else {
      optimizer = CPLEX;
    }
  }
  if (optimizer == CP) {
    for (int i = 0, num_vars = p.num_vars(); i < num_vars; ++i) {
      if (p.var(i).type() == mp::var::CONTINUOUS)
        throw Error("CP Optimizer doesn't support continuous variables");
    }
  }

  if (GetOption(USENUMBEROF) != 0)
    converter_flags_ |= MPToConcertConverter::USENUMBEROF;
  if (GetOption(DEBUGEXPR) != 0)
    converter_flags_ |= MPToConcertConverter::DEBUG;
  converter_.reset(new MPToConcertConverter(env_, converter_flags_));
}

void IlogCPSolver::AddVariable(Problem::Variable var) {
  converter_->AddVariable(var);
}
void IlogCPSolver::AddCommonExpression(Problem::CommonExpr cexpr) {
  converter_->AddCommonExpression(cexpr);
}
void IlogCPSolver::AddObjective(Problem::Objective obj) {
  converter_->AddObjective(obj);
}
void IlogCPSolver::AddAlgebraicConstraint(Problem::AlgebraicCon con) {
  converter_->AddAlgebraicConstraint(con);
}
void IlogCPSolver::AddLogicalConstraint(Problem::LogicalCon lcon) {
  converter_->AddLogicalConstraint(lcon);
}
void IlogCPSolver::FinishProblemModificationPhase() {
  converter_->FinishProblemModificationPhase();
}

void IlogCPSolver::Convert(Problem &p) {
  InitProblemModificationPhase(p);
  converter_->Convert(p);
}

void IlogCPSolver::Resolve(Problem &p, SolutionHandler &sh) {
  try {
    IloAlgorithm &cp_alg = cp_;
    (optimizer == CP ? cp_alg : cplex_).extract(converter_->model());
  } catch (IloAlgorithm::CannotExtractException &e) {
    const IloExtractableArray &extractables = e.getExtractables();
    if (extractables.getSize() == 0)
      throw;
    throw mp::MakeUnsupportedError("{}", extractables[0].getObject());
  }

  if (optimizer == CP)
    SolveWithCP(p, *converter_, stats, sh);
  else
    SolveWithCPLEX(p, *converter_, stats, sh);
  double output_time = GetTimeAndReset(stats.time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          stats.setup_time, stats.solution_time, output_time);
  }
}


SolverPtr create_ilogcp(const char *) { return SolverPtr(new IlogCPSolver()); }
}

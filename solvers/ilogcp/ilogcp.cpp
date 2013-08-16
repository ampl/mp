/*-------------------------------------------------------------------------*/
/* AMPL/IBM ILOG CP Optimizer driver                         Robert Fourer */
/*                                                                         */
/* Name           : ilogcp.cpp                                             */
/* Title          : AMPL/IBM ILOG CP Optimizer driver                      */
/* By             : Robert Fourer                                          */
/* Date           : October 2000                                           */
/*                                                                         */
/* A driver to link AMPL linear integer programs with ILOG Concert 1.0     */
/* October 2000: Linear/Nonlinear version                                  */
/* June 2012:    Updated to Concert 12.4 (Victor Zverovich)                */
/*-------------------------------------------------------------------------*/
//
// Possible improvements: Some sort of variable preference mechanism.
//
// Reference: "Extending an Algebraic Modeling Language to
// Support Constraint Programming" by Robert Fourer and David M. Gay,
// INFORMS Journal on Computing, Fall 2002, vol. 14, no. 4, 322-344
// (http://joc.journal.informs.org/content/14/4/322).

#include "ilogcp.h"
#include "ilogcp_date.h"

#include <cctype>
#include <cstdlib>

using std::vector;

#ifndef ILOGCP_NO_VERS
static char xxxvers[] = "ilogcp_options\0\n"
  "AMPL/IBM ILOG CP Optimizer Driver Version " qYYYYMMDD "\n";
#endif

#ifndef M_PI
# define M_PI 3.14159265358979323846
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

void RequireNonzeroConstRHS(ampl::BinaryExpr e, const std::string &func_name) {
  ampl::NumericConstant num = ampl::Cast<ampl::NumericConstant>(e.rhs());
  if (!num || num.value() != 0) {
    throw ampl::UnsupportedExprError::CreateFromExprString(
        func_name + " with nonzero second parameter");
  }
}

inline IloInt CastToInt(double value) {
  IloInt int_value = static_cast<int>(value);
  if (int_value != value) {
    throw ampl::Error(str(
        fmt::Format("value {} can't be represented as int") << value));
  }
  return int_value;
}

ampl::OptionError GetOptionValueError(
    fmt::StringRef name, fmt::StringRef message) {
  throw ampl::OptionError(fmt::Format(
      "Can't get value of option {}: {}") << name.c_str() << message.c_str());
}
}

namespace ampl {

Optimizer::Optimizer(IloEnv env, const Problem &p) : cons_(env, p.num_cons()) {}

Optimizer::~Optimizer() {}

CPLEXOptimizer::CPLEXOptimizer(IloEnv env, const Problem &p)
: Optimizer(env, p), cplex_(env), aborter_(env), started_(false) {
  cplex_.setParam(IloCplex::MIPDisplay, 0);
  cplex_.use(aborter_);
}

bool CPLEXOptimizer::FindNextSolution() {
  if (!started_)
    return false;
  started_ = false;
  cplex_.solve();
  return true;
}

void CPLEXOptimizer::GetSolutionInfo(
    fmt::Formatter &format_message, vector<double> &dual_values) const {
  if (cplex_.isMIP()) {
    format_message("{} nodes, ") << cplex_.getNnodes();
  } else {
    IloRangeArray cons = Optimizer::cons();
    IloInt num_cons = cons.getSize();
    dual_values.resize(num_cons);
    for (IloInt i = 0; i < num_cons; ++i)
      dual_values[i] = cplex_.getDual(cons[i]);
  }
  format_message("{} iterations") << cplex_.getNiterations();
}

CPOptimizer::CPOptimizer(IloEnv env, const Problem &p) :
    Optimizer(env, p), solver_(env) {
  solver_.setIntParameter(IloCP::LogVerbosity, IloCP::Quiet);
  if (p.num_objs() == 0)
    solver_.setIntParameter(IloCP::SolutionLimit, 1);
}

void CPOptimizer::GetSolutionInfo(
    fmt::Formatter &format_message, vector<double> &) const {
  format_message("{} choice points, {} fails")
      << solver_.getInfo(IloCP::NumberOfChoicePoints)
      << solver_.getInfo(IloCP::NumberOfFails);
}

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

IloExpr NLToConcertConverter::VisitPLTerm(PiecewiseLinearTerm t) {
  IloNumArray slopes(env_), breakpoints(env_);
  int num_breakpoints = t.num_breakpoints();
  for (int i = 0; i < num_breakpoints; ++i) {
    slopes.add(t.slope(i));
    breakpoints.add(t.breakpoint(i));
  }
  slopes.add(t.slope(num_breakpoints));
  return IloPiecewiseLinear(vars_[t.var_index()], breakpoints, slopes, 0, 0);
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
   Solver<IlogCPSolver>("ilogcp", 0, YYYYMMDD), gotopttype_(false) {
  options_[DEBUGEXPR] = 0;
  options_[OPTIMIZER] = AUTO;
  options_[TIMING] = 0;
  options_[USENUMBEROF] = 1;

  set_long_name(fmt::Format("ilogcp {}.{}.{}")
      << IloConcertVersion::_ILO_MAJOR_VERSION
      << IloConcertVersion::_ILO_MINOR_VERSION
      << IloConcertVersion::_ILO_TECH_VERSION);
  set_version(fmt::Format("AMPL/IBM ILOG CP Optimizer [{} {}.{}.{}]")
      << IloConcertVersion::_ILO_NAME << IloConcertVersion::_ILO_MAJOR_VERSION
      << IloConcertVersion::_ILO_MINOR_VERSION
      << IloConcertVersion::_ILO_TECH_VERSION);

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

  AddOption(OptionPtr(new EnumCPOption("alldiffinferencelevel",
      "Inference level for 'alldiff' constraints.  Possible values:\n"
      "      0 = default\n"
      "      1 = low\n"
      "      2 = basic\n"
      "      3 = medium\n"
      "      4 = extended\n",
      this, IloCP::AllDiffInferenceLevel, IloCP::Default, InferenceLevels)));

  AddOption(OptionPtr(new IntCPOption("branchlimit",
      "Limit on the number of branches made before "
      "terminating a search.  Default = no limit.",
      this, IloCP::BranchLimit)));

  AddOption(OptionPtr(new IntCPOption("choicepointlimit",
      "Limit on the number of choice points created"
      "before terminating a search.  Default = no limit.",
      this, IloCP::ChoicePointLimit)));

  AddOption(OptionPtr(new EnumCPOption("constraintaggregation",
      "0 or 1 (default 1):  Whether to aggregate basic constraints.",
      this, IloCP::ConstraintAggregation, IloCP::Off, Flags)));

  AddIntOption("debugexpr",
      "0 or 1 (default 0):  Whether to print debugging "
      "information for expression trees.",
      &IlogCPSolver::GetBoolOption, &IlogCPSolver::SetBoolOption, DEBUGEXPR);

  AddOption(OptionPtr(new EnumCPOption("defaultinferencelevel",
      "Default inference level for constraints.  Possible values:\n"
      "      1 = low\n"
      "      2 = basic\n"
      "      3 = medium\n"
      "      4 = extended\n",
      this, IloCP::DefaultInferenceLevel, IloCP::Default, InferenceLevels)));

  AddOption(OptionPtr(new EnumCPOption("distributeinferencelevel",
      "Inference level for 'distribute' constraints.  Possible values:\n"
      "      0 = default\n"
      "      1 = low\n"
      "      2 = basic\n"
      "      3 = medium\n"
      "      4 = extended\n",
      this, IloCP::DistributeInferenceLevel, IloCP::Default, InferenceLevels)));

  AddOption(OptionPtr(new EnumCPOption("dynamicprobing",
      "Use probing during search.  Possible values:\n"
      "     -1 = auto (default)\n"
      "      0 = off\n"
      "      1 = on\n",
      this, IloCP::DynamicProbing, IloCP::Off, Flags, true)));

  AddDblOption("dynamicprobingstrength",
      "Effort dedicated to dynamic probing as a factor "
      "of the total search effort.  Default = 0.03.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::DynamicProbingStrength);

  AddOption(OptionPtr(new IntCPOption("faillimit",
      "Limit on the number of failures allowed before terminating a search.  "
      "Default = no limit",
      this, IloCP::FailLimit)));

  AddOption(OptionPtr(new IntCPOption("logperiod",
      "Specifies how often the information in the search log is displayed.",
      this, IloCP::LogPeriod)));

  AddOption(OptionPtr(new EnumCPOption("logverbosity",
      "Verbosity of the search log.  Possible values:\n"
      "      0 = quiet (default)\n"
      "      1 = terse\n"
      "      2 = normal\n"
      "      3 = verbose\n",
      this, IloCP::LogVerbosity, IloCP::Quiet, Verbosities)));

  AddIntOption<int>("mipdisplay",
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

  AddIntOption<int>("mipinterval",
      "Frequency of node logging for mipdisplay 2 or 3. Default = 1.",
      &IlogCPSolver::GetCPLEXIntOption, &IlogCPSolver::SetCPLEXIntOption,
      IloCplex::MIPInterval);

  AddOption(OptionPtr(new IntCPOption("multipointnumberofsearchpoints",
      "Number of solutions for the multi-point search "
      "algorithm.  Default = 30.",
      this, IloCP::MultiPointNumberOfSearchPoints)));

  AddDblOption("optimalitytolerance",
      "Absolute tolerance on the objective value.  Default = 0.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::OptimalityTolerance);

  AddStrOption("optimizer",
      "Specifies which optimizer to use.  Possible values:\n"
      "      auto  = CP Optimizer if the problem has\n"
      "              nonlinear objective/constraints\n"
      "              or logical constraints, CPLEX\n"
      "              otherwise (default)\n"
      "      cp    = CP Optimizer\n"
      "      cplex = CPLEX Optimizer\n",
      &IlogCPSolver::GetOptimizer, &IlogCPSolver::SetOptimizer);

  AddOption(OptionPtr(new EnumCPOption("outlev",
      "Synonym for \"logverbosity\".",
      this, IloCP::LogVerbosity, IloCP::Quiet, Verbosities)));

  AddOption(OptionPtr(new EnumCPOption("propagationlog",
      "Level of propagation trace reporting.  Possible values:\n"
      "      0 = quiet (default)\n"
      "      1 = terse\n"
      "      2 = normal\n"
      "      3 = verbose\n",
      this, IloCP::PropagationLog, IloCP::Quiet, Verbosities)));

  AddOption(OptionPtr(new IntCPOption("randomseed",
      "Seed for the random number generator.  Default = 0.",
      this, IloCP::RandomSeed)));

  AddDblOption("relativeoptimalitytolerance",
      "Relative tolerance on the objective value.  Default = 1e-4.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::RelativeOptimalityTolerance);

  AddOption(OptionPtr(new IntCPOption("restartfaillimit",
      "Number of failures allowed before restarting  search.  Default = 100.",
      this, IloCP::RestartFailLimit)));

  AddDblOption("restartgrowthfactor",
      "Increase of the number of allowed failures "
      "before restarting search.  Default = 1.05.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::RestartGrowthFactor);

  AddOption(OptionPtr(new EnumCPOption("searchtype",
      "Type of search used for solving a problem.  Possible values:\n"
      "      0 = depthfirst\n"
      "      1 = restart (default)\n"
      "      2 = multipoint\n",
      this, IloCP::SearchType, IloCP::DepthFirst, SearchTypes, true)));

  AddOption(OptionPtr(new IntCPOption("solutionlimit",
      "Limit on the number of feasible solutions found before terminating "
      "a search.  Leaving the solution limit unspecified will make the "
      "optimizer search for an optimal solution if there is an objective "
      "function or for a feasible solution otherwise.",
      this, IloCP::SolutionLimit)));

  AddOption(OptionPtr(new EnumCPOption("temporalrelaxation",
      "0 or 1 (default 1):  Whether to use temporal relaxation.",
      this, IloCP::TemporalRelaxation, IloCP::Off, Flags)));

  AddDblOption("timelimit",
      "Limit on the CPU time spent solving before "
      "terminating a search.  Default = no limit.",
      &IlogCPSolver::GetCPDblOption, &IlogCPSolver::SetCPDblOption,
      IloCP::TimeLimit);

  AddOption(OptionPtr(new EnumCPOption("timemode",
      "Specifies how the time is measured in CP Optimizer.  Possible values:\n"
      "      0 = cputime (default)\n"
      "      1 = elapsedtime\n",
      this, IloCP::TimeMode, IloCP::CPUTime, TimeModes)));

  AddIntOption("timing",
      "0 or 1 (default 0):  Whether to display timings for the run.\n",
      &IlogCPSolver::GetBoolOption, &IlogCPSolver::SetBoolOption,
      IlogCPSolver::TIMING);

  AddIntOption("usenumberof",
      "0 or 1 (default 1):  Whether to consolidate 'numberof' expressions "
      "by use of IloDistribute constraints.",
      &IlogCPSolver::GetBoolOption, &IlogCPSolver::SetBoolOption,
      IlogCPSolver::USENUMBEROF);

  AddOption(OptionPtr(new EnumCPOption("workers",
      "Number of workers to run in parallel to solve a problem.  "
      "In addition to numeric values this option accepts the value "
      "\"auto\" since CP Optimizer version 12.3.  Default = 1.",
      this, IloCP::Workers, 0, 0, true)));
}

IlogCPSolver::~IlogCPSolver() {
  env_.end();
}

CPOptimizer *IlogCPSolver::GetCPForOption(fmt::StringRef option_name) const {
  CPOptimizer *cp = dynamic_cast<CPOptimizer*>(optimizer_.get());
  if (!cp) {
    throw OptionError(
        fmt::Format("Invalid option {} for CPLEX optimizer")
                << option_name.c_str());
  }
  return cp;
}

CPLEXOptimizer *IlogCPSolver::GetCPLEXForOption(
    fmt::StringRef option_name) const {
  CPLEXOptimizer *cplex = dynamic_cast<CPLEXOptimizer*>(optimizer_.get());
  if (!cplex) {
    throw OptionError(
        fmt::Format("Invalid option {} for CP optimizer")
                << option_name.c_str());
  }
  return cplex;
}

std::string IlogCPSolver::GetOptimizer(const char *) const {
  switch (options_[OPTIMIZER]) {
  default:
    assert(false);
    // Fall through.
  case AUTO:  return "auto";
  case CP:    return "cp";
  case CPLEX: return "cplex";
  }
}

void IlogCPSolver::SetOptimizer(const char *name, const char *value) {
  int opt = 0;
  if (strcmp(value, "auto") == 0)
    opt = AUTO;
  else if (strcmp(value, "cp") == 0)
    opt = CP;
  else if (strcmp(value, "cplex") == 0)
    opt = CPLEX;
  else
    throw InvalidOptionValue(name, value);
  if (!gotopttype_)
    options_[OPTIMIZER] = opt;
}

void IlogCPSolver::SetBoolOption(const char *name, int value, Option opt) {
  if (!gotopttype_)
    return;
  if (value != 0 && value != 1)
    throw InvalidOptionValue(name, value);
  options_[opt] = value;
}

std::string IlogCPSolver::EnumCPOption::GetValue() const {
  int value = 0;
  try {
    value = solver_.GetCPForOption(name())->solver().getParameter(param_);
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

void IlogCPSolver::EnumCPOption::SetValue(const char *value) {
  if (!solver_.gotopttype_)
    return;
  CPOptimizer *cp = solver_.GetCPForOption(name());
  try {
    char *end = 0;
    long intval = std::strtol(value, &end, 0);
    if (!*end) {
      if (intval != -1 || !accepts_auto_)
        intval += start_;
      cp->solver().setParameter(param_, intval);
      return;
    }
    if (values_) {
      // Search for a value in the list of known values.
      // Use linear search since the number of values is small.
      for (int i = 0; values_[i]; ++i) {
        if (strcmp(value, values_[i]) == 0) {
          cp->solver().setParameter(param_, i + start_);
          return;
        }
      }
    }
    if (accepts_auto_ && strcmp(value, "auto") == 0) {
      cp->solver().setParameter(param_, IloCP::Auto);
      return;
    }
  } catch (const IloException &) {}
  throw InvalidOptionValue(name(), value);
}

int IlogCPSolver::IntCPOption::GetValue() const {
  try {
    return solver_.GetCPForOption(name())->solver().getParameter(param_);
  } catch (const IloException &e) {
    throw GetOptionValueError(name(), e.getMessage());
  }
}

void IlogCPSolver::IntCPOption::SetValue(int value) {
  if (!solver_.gotopttype_)
    return;
  try {
    solver_.GetCPForOption(name())->solver().setParameter(param_, value);
  } catch (const IloException &) {
    throw InvalidOptionValue(name(), value);
  }
}

double IlogCPSolver::GetCPDblOption(
    const char *name, IloCP::NumParam param) const {
  try {
    return GetCPForOption(name)->solver().getParameter(param);
  } catch (const IloException &e) {
    throw GetOptionValueError(name, e.getMessage());
  }
}

void IlogCPSolver::SetCPDblOption(
    const char *name, double value, IloCP::NumParam param) {
  if (!gotopttype_)
    return;
  try {
    GetCPForOption(name)->solver().setParameter(param, value);
  } catch (const IloException &) {
    throw InvalidOptionValue(name, value);
  }
}

int IlogCPSolver::GetCPLEXIntOption(const char *name, int param) const {
  // Use CPXgetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  cpxenv *env = GetCPLEXForOption(name)->cplex().getImpl()->getCplexEnv();
  int value = 0;
  int result = CPXgetintparam(env, param, &value);
  if (result != 0)
    throw GetOptionValueError(name, fmt::Format("CPLEX error = {}") << result);
  return value;
}

void IlogCPSolver::SetCPLEXIntOption(const char *name, int value, int param) {
  if (!gotopttype_)
    return;
  // Use CPXsetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  cpxenv *env = GetCPLEXForOption(name)->cplex().getImpl()->getCplexEnv();
  if (CPXsetintparam(env, param, value) != 0)
    throw InvalidOptionValue(name, value);
}

void IlogCPSolver::CreateOptimizer(const Problem &p) {
  int &opt = options_[OPTIMIZER];
  if (opt == AUTO) {
    opt = CPLEX;
    if (p.num_nonlinear_objs() + p.num_nonlinear_cons() +
        p.num_logical_cons() != 0) {
      opt = CP;
    }
  }
  if (opt == CPLEX)
    optimizer_.reset(new CPLEXOptimizer(env_, p));
  else
    optimizer_.reset(new CPOptimizer(env_, p));
}

bool IlogCPSolver::ParseOptions(char **argv, unsigned flags) {
  // Get optimizer type.
  gotopttype_ = false;
  if (!BasicSolver::ParseOptions(argv, BasicSolver::NO_OPTION_ECHO))
    return false;
  CreateOptimizer(problem());

  // Parse remaining options.
  gotopttype_ = true;
  return BasicSolver::ParseOptions(argv, flags);
}

void IlogCPSolver::Solve(Problem &p) {
  double start_time = xectim_();

  // Set up optimization problem in ILOG Concert.

  int num_continuous_vars = p.num_continuous_vars();
  if (num_continuous_vars != 0 && GetOption(OPTIMIZER) == CP)
    throw Error("CP Optimizer doesn't support continuous variables");

  if (!optimizer_.get())
    CreateOptimizer(p);

  int num_vars = p.num_vars();
  IloNumVarArray vars(env_, num_vars);
  for (int j = 0; j < num_continuous_vars; j++)
    vars[j] = IloNumVar(env_, p.var_lb(j), p.var_ub(j), ILOFLOAT);
  for (int j = num_continuous_vars; j < num_vars; j++)
    vars[j] = IloNumVar(env_, p.var_lb(j), p.var_ub(j), ILOINT);

  NLToConcertConverter converter(env_, vars,
      GetOption(USENUMBEROF) != 0, GetOption(DEBUGEXPR) != 0);
  IloModel model = converter.model();
  int num_objs = p.num_objs();
  if (num_objs > 0) {
    NumericExpr expr(p.nonlinear_obj_expr(0));
    NumericConstant constant(Cast<NumericConstant>(expr));
    IloExpr ilo_expr(env_, constant ? constant.value() : 0);
    if (p.num_nonlinear_objs() > 0 && !constant)
      ilo_expr += converter.Visit(expr);
    LinearObjExpr linear = p.linear_obj_expr(0);
    for (LinearObjExpr::iterator
        i = linear.begin(), end = linear.end(); i != end; ++i) {
      ilo_expr += i->coef() * vars[i->var_index()];
    }
    IloObjective obj(env_, ilo_expr,
        p.obj_type(0) == MIN ? IloObjective::Minimize : IloObjective::Maximize);
    IloAdd(model, obj);
  }

  if (int n_cons = p.num_cons()) {
    IloRangeArray cons(optimizer_->cons());
    for (int i = 0; i < n_cons; ++i) {
      IloExpr expr(env_);
      LinearConExpr linear = p.linear_con_expr(i);
      for (LinearConExpr::iterator
          j = linear.begin(), end = linear.end(); j != end; ++j) {
        expr += j->coef() * vars[j->var_index()];
      }
      if (i < p.num_nonlinear_cons())
        expr += converter.Visit(p.nonlinear_con_expr(i));
      cons[i] = (p.con_lb(i) <= expr <= p.con_ub(i));
    }
    model.add(cons);
  }

  if (int n_lcons = p.num_logical_cons()) {
    IloConstraintArray cons(env_, n_lcons);
    for (int i = 0; i < n_lcons; ++i)
      cons[i] = converter.Visit(p.logical_con_expr(i));
    model.add(cons);
  }

  converter.FinishBuildingNumberOf();

  double define_time = xectim_() - start_time;
  start_time = xectim_();

  // Solve the problem.
  IloAlgorithm alg(optimizer_->algorithm());
  try {
    alg.extract(model);
  } catch (IloAlgorithm::CannotExtractException &e) {
    const IloExtractableArray &extractables = e.getExtractables();
    if (extractables.getSize() == 0)
      throw;
    throw UnsupportedExprError::CreateFromExprString(
        str(fmt::Format("{}") << extractables[0]));
  }
  SignalHandler sh(*this, optimizer_.get());
  vector<double> solution, dual_solution;
  double setup_time = xectim_() - start_time;
  double solve_time = 0, output_time = 0;
  start_time = xectim_();
  optimizer_->StartSearch();
  solve_time += xectim_() - start_time;
  start_time = xectim_();
  for (bool succeeded = true, first = true; succeeded; first = false) {
    start_time = xectim_();
    succeeded = optimizer_->FindNextSolution();
    solve_time += xectim_() - start_time;
    start_time = xectim_();
    if (!succeeded && !first)
      continue;
    // Convert solution status.
    bool has_solution = false;
    int solve_code = 0;
    const char *status;
    switch (alg.getStatus()) {
    default:
      // Fall through.
    case IloAlgorithm::Unknown:
      if (sh.stop()) {
        solve_code = 600;
        status = "interrupted";
      } else {
        solve_code = 501;
        status = "unknown solution status";
      }
      break;
    case IloAlgorithm::Feasible:
      has_solution = true;
      if (sh.stop()) {
        solve_code = 600;
        status = "interrupted";
      } else {
        solve_code = 100;
        status = "feasible solution";
      }
      break;
    case IloAlgorithm::Optimal:
      has_solution = true;
      solve_code = 0;
      status = "optimal solution";
      break;
    case IloAlgorithm::Infeasible:
      solve_code = 200;
      status = "infeasible problem";
      break;
    case IloAlgorithm::Unbounded:
      solve_code = 300;
      status = "unbounded problem";
      break;
    case IloAlgorithm::InfeasibleOrUnbounded:
      solve_code = 201;
      status = "infeasible or unbounded problem";
      break;
    case IloAlgorithm::Error:
      solve_code = 500;
      status = "error";
      break;
    }
    p.set_solve_code(solve_code);

    fmt::Formatter format_message;
    format_message("{}: {}\n") << long_name() << status;
    double obj_value = std::numeric_limits<double>::quiet_NaN();
    solution.clear();
    dual_solution.clear();
    if (has_solution) {
      solution.resize(num_vars);
      for (int j = 0, n = p.num_vars(); j < n; ++j) {
        IloNumVar &v = vars[j];
        solution[j] = alg.isExtracted(v) ? alg.getValue(v) : v.getLB();
      }
      optimizer_->GetSolutionInfo(format_message, dual_solution);
      if (num_objs > 0) {
        obj_value = alg.getObjValue();
        format_message(", objective {}") << ObjPrec(obj_value);
      }
    }
    HandleSolution(format_message.c_str(), solution.empty() ? 0 : &solution[0],
        dual_solution.empty() ? 0 : &dual_solution[0], obj_value);
    output_time += xectim_() - start_time;
  }
  start_time = xectim_();
  optimizer_->EndSearch();
  solve_time += xectim_() - start_time;

  if (GetOption(TIMING)) {
    std::cerr << "\n"
        << "Define = " << define_time + read_time() << "\n"
        << "Setup =  " << setup_time << "\n"
        << "Solve =  " << solve_time << "\n"
        << "Output = " << output_time << "\n";
  }
}
}

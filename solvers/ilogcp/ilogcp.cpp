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
#include <algorithm>
#include <iostream>

#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include "getstub.h"
#include "nlp.h"
#include "opcode.hd"

using namespace std;

using ampl::NumberOfExpr;

static char xxxvers[] = "ilogcp_options\0\n"
  "AMPL/IBM ILOG CP Optimizer Driver Version " qYYYYMMDD "\n";

// for suppressing "String literal to char*" warnings
#define CSTR(s) const_cast<char*>(s)

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

namespace {
struct DriverOptionInfo : Option_Info {
  ampl::IlogCPDriver *driver;
};

char *skip_nonspace(char *s) {
  while (*s && !isspace(*s))
    ++s;
  return s;
}

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

// Information about a constraint programming solver option.
struct CPOptionInfo {
  IloCP::IntParam param;
  int start;           // start value for the enumerated options
  bool accepts_auto;   // true if the option accepts IloCP::Auto value
  const char **values; // string values for enum options
};

#define CP_INT_OPTION_FULL(name, start, accepts_auto, values) \
  const CPOptionInfo name = {IloCP::name, start, accepts_auto, values};
#define CP_INT_OPTION(name) CP_INT_OPTION_FULL(name, 0, false, 0)
#define CP_ENUM_OPTION(name, start, values) \
  CP_INT_OPTION_FULL(name, start, false, values)

CP_ENUM_OPTION(AllDiffInferenceLevel, IloCP::Default, InferenceLevels)
CP_INT_OPTION(BranchLimit)
CP_INT_OPTION(ChoicePointLimit)
CP_ENUM_OPTION(ConstraintAggregation, IloCP::Off, Flags)
CP_ENUM_OPTION(DefaultInferenceLevel, IloCP::Default, InferenceLevels)
CP_ENUM_OPTION(DistributeInferenceLevel, IloCP::Default, InferenceLevels)
CP_INT_OPTION_FULL(DynamicProbing, IloCP::Off, true, Flags)
CP_INT_OPTION(FailLimit)
CP_INT_OPTION(LogPeriod)
CP_ENUM_OPTION(LogVerbosity, IloCP::Quiet, Verbosities)
CP_INT_OPTION(MultiPointNumberOfSearchPoints)
CP_ENUM_OPTION(PropagationLog, IloCP::Quiet, Verbosities)
CP_INT_OPTION(RandomSeed)
CP_INT_OPTION(RestartFailLimit)
CP_INT_OPTION_FULL(SearchType, IloCP::DepthFirst, true, SearchTypes)
CP_INT_OPTION(SolutionLimit)
CP_ENUM_OPTION(TemporalRelaxation, IloCP::Off, Flags)
CP_ENUM_OPTION(TimeMode, IloCP::CPUTime, TimeModes)
CP_INT_OPTION_FULL(Workers, 0, true, 0)

class SameExpr {
private:
  NumberOfExpr expr_;
  unsigned num_args_;
  
public:
  SameExpr(NumberOfExpr e) :
  expr_(e), num_args_(std::distance(e.begin(), e.end())) {}
  
  // Returns true if the stored expression is the same as the argument's
  // expression.
  bool operator()(const ampl::NumberOf &nof) const;
};

bool SameExpr::operator()(const ampl::NumberOf &nof) const {
  if (nof.num_vars() != num_args_)
    return false;
  
  for (NumberOfExpr::iterator i = expr_.begin(), end = expr_.end(),
       j = nof.expr().begin(); i != end; ++i, ++j) {
    if (!Equal(*i, *j))
      return false;
  }
  return true;
}
}

namespace ampl {
  
IloIntVar NumberOf::Add(real value, IloEnv env) {
  for (int i = 0, n = values_.getSize(); i < n; ++i)
    if (values_[i] == value)
      return cards_[i];
  values_.add(value);
  IloIntVar cardVar(env, IloIntMin, IloIntMax);
  cards_.add(cardVar);
  return cardVar;
}

Optimizer::Optimizer(IloEnv env, ASL_fg *asl) :
  vars_(env, n_var), cons_(env, n_con) {}

Optimizer::~Optimizer() {}

void CPLEXOptimizer::set_option(const void *key, int value) {
  // Use CPXsetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  if (CPXsetintparam(cplex_.getImpl()->getCplexEnv(),
      reinterpret_cast<size_t>(key), value) != 0) {
    throw IloWrongUsage();
  }
}

void CPLEXOptimizer::set_option(const void *key, double value) {
  cplex_.setParam(
      static_cast<IloCplex::NumParam>(reinterpret_cast<size_t>(key)), value);
}

void CPLEXOptimizer::get_solution(ASL_fg *asl, char *message,
    std::vector<double> &primal, std::vector<double> &dual) const {
  IloNum objValue = cplex_.getObjValue();
  primal.resize(n_var);
  IloNumVarArray vars = Optimizer::vars();
  for (int j = 0; j < n_var; j++) primal[j] = cplex_.getValue(vars[j]);
  if (cplex_.isMIP()) {
    message += g_fmtop(message, cplex_.getNnodes());
    message += Sprintf(message, " nodes, ");
    message += g_fmtop(message, cplex_.getNiterations());
    message += Sprintf(message, " iterations, objective ");
    g_fmtop(message, objValue);
  } else {
    message += g_fmtop(message, cplex_.getNiterations());
    message += Sprintf(message, " iterations, objective ");
    g_fmtop(message, objValue);
    dual.resize(n_con);
    IloRangeArray cons = Optimizer::cons();
    for (int i = 0; i < n_con; i++)
      dual[i] = cplex_.getDual(cons[i]);
  }
}

void CPOptimizer::set_option(const void *key, int value) {
  const CPOptionInfo *info = static_cast<const CPOptionInfo*>(key);
  if (value != -1 || !info->accepts_auto)
    value += info->start;
  solver_.setParameter(info->param, value);
}

void CPOptimizer::set_option(const void *key, double value) {
  solver_.setParameter(
      static_cast<IloCP::NumParam>(reinterpret_cast<size_t>(key)), value);
}

void CPOptimizer::get_solution(ASL_fg *asl, char *message,
    std::vector<double> &primal, std::vector<double> &) const {
  message += g_fmtop(message, solver_.getNumberOfChoicePoints());
  message += Sprintf(message, " choice points, ");
  message += g_fmtop(message, solver_.getNumberOfFails());
  message += Sprintf(message, " fails");
  if (n_obj > 0) {
    message += Sprintf(message, ", objective ");
    g_fmtop(message, solver_.getValue(obj()));
  }
  primal.resize(n_var);
  IloNumVarArray vars = Optimizer::vars();
  for (int j = 0; j < n_var; j++)
    primal[j] = solver_.getValue(vars[j]);
}

#define SPACE "                                "

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
keyword IlogCPDriver::keywords_[] = {
  // The options must be in alphabetical order.

  KW(CSTR("alldiffinferencelevel"),
      IlogCPDriver::set_cp_int_option, &AllDiffInferenceLevel,
      CSTR("Inference level for 'alldiff' constraints.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = default\n"
          SPACE "      1 = low\n"
          SPACE "      2 = basic\n"
          SPACE "      3 = medium\n"
          SPACE "      4 = extended\n")),

  KW(CSTR("branchlimit"),
      IlogCPDriver::set_cp_int_option, &BranchLimit,
      CSTR("Limit on the number of branches made before\n"
          SPACE "terminating a search.  Default = no limit.\n")),

  KW(CSTR("choicepointlimit"),
      IlogCPDriver::set_cp_int_option, &ChoicePointLimit,
      CSTR("Limit on the number of choice points created\n"
          SPACE "before terminating a search. Default = no limit.\n")),

  KW(CSTR("constraintaggregation"),
      IlogCPDriver::set_cp_int_option, &ConstraintAggregation,
      CSTR("0 or 1 (default 1):  Whether to aggregate basic\n"
          SPACE "constraints.\n")),

  KW(CSTR("debugexpr"), IlogCPDriver::set_int_option, IlogCPDriver::DEBUGEXPR,
      CSTR("0 or 1 (default 0):  Whether to print debugging\n"
          SPACE "information for expression trees.\n")),

  KW(CSTR("defaultinferencelevel"),
      IlogCPDriver::set_cp_int_option, &DefaultInferenceLevel,
      CSTR("Default inference level for constraints.\n"
          SPACE "Possible values:\n"
          SPACE "      1 = low\n"
          SPACE "      2 = basic\n"
          SPACE "      3 = medium\n"
          SPACE "      4 = extended\n")),

  KW(CSTR("distributeinferencelevel"),
      IlogCPDriver::set_cp_int_option, &DistributeInferenceLevel,
      CSTR("Inference level for 'distribute' constraints.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = default\n"
          SPACE "      1 = low\n"
          SPACE "      2 = basic\n"
          SPACE "      3 = medium\n"
          SPACE "      4 = extended\n")),

  KW(CSTR("dynamicprobing"),
      IlogCPDriver::set_cp_int_option, &DynamicProbing,
      CSTR("Use probing during search.  Possible values:\n"
          SPACE "     -1 = auto (default)\n"
          SPACE "      0 = off\n"
          SPACE "      1 = on\n")),

  KW(CSTR("dynamicprobingstrength"),
      IlogCPDriver::set_cp_dbl_option, IloCP::DynamicProbingStrength,
      CSTR("Effort dedicated to dynamic probing as a factor\n"
          SPACE "of the total search effort.  Default = 0.03.\n")),

  KW(CSTR("faillimit"),
      IlogCPDriver::set_cp_int_option, &FailLimit,
      CSTR("Limit on the number of failures allowed before\n"
          SPACE "terminating a search.  Default = no limit.\n")),

  KW(CSTR("logperiod"),
      IlogCPDriver::set_cp_int_option, &LogPeriod,
      CSTR("Specifies how often the information in the\n"
          SPACE "search log is displayed.\n")),

  KW(CSTR("logverbosity"), IlogCPDriver::set_cp_int_option, &LogVerbosity,
      CSTR("Verbosity of the search log.  Possible values:\n"
          SPACE "      0 = quiet (default)\n"
          SPACE "      1 = terse\n"
          SPACE "      2 = normal\n"
          SPACE "      3 = verbose\n")),

  KW(CSTR("mipdisplay"), IlogCPDriver::set_cplex_int_option, IloCplex::MIPDisplay,
      CSTR("Frequency of displaying branch-and-bound\n"
          SPACE "information (for optimizing integer variables):\n"
          SPACE "      0 (default) = never\n"
          SPACE "      1 = each integer feasible solution\n"
          SPACE "      2 = every \"mipinterval\" nodes\n"
          SPACE "      3 = every \"mipinterval\" nodes plus\n"
          SPACE "          information on LP relaxations\n"
          SPACE "          (as controlled by \"display\")\n"
          SPACE "      4 = same as 2, plus LP relaxation info.\n"
          SPACE "      5 = same as 2, plus LP subproblem info.\n")),

  KW(CSTR("mipinterval"), IlogCPDriver::set_cplex_int_option, IloCplex::MIPInterval,
      CSTR("Frequency of node logging for mipdisplay 2 or 3.\n"
          SPACE "Default = 1.\n")),

  KW(CSTR("multipointnumberofsearchpoints"),
      IlogCPDriver::set_cp_int_option, &MultiPointNumberOfSearchPoints,
      CSTR("Number of solutions for the multi-point search\n"
          SPACE "algorithm.  Default = 30.\n")),

  KW(CSTR("optimalitytolerance"),
      IlogCPDriver::set_cp_dbl_option, IloCP::OptimalityTolerance,
      CSTR("Absolute tolerance on the objective value.\n"
          SPACE "Default = 0.\n")),

  KW(CSTR("optimizer"), IlogCPDriver::set_optimizer, 0,
      CSTR("Specifies which optimizer to use.\n"
          SPACE "Possible values:\n"
          SPACE "      auto  = CP Optimizer if the problem has\n"
          SPACE "              nonlinear objective/constraints\n"
          SPACE "              or logical constraints, CPLEX\n"
          SPACE "              otherwise (default)\n"
          SPACE "      cp    = CP Optimizer\n"
          SPACE "      cplex = CPLEX Optimizer\n")),

  KW(CSTR("outlev"), IlogCPDriver::set_cp_int_option,
      &LogVerbosity, CSTR("Synonym for \"logverbosity\".\n")),

  KW(CSTR("propagationlog"), IlogCPDriver::set_cp_int_option, &PropagationLog,
      CSTR("Level of propagation trace reporting.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = quiet (default)\n"
          SPACE "      1 = terse\n"
          SPACE "      2 = normal\n"
          SPACE "      3 = verbose\n")),

  KW(CSTR("randomseed"), IlogCPDriver::set_cp_int_option, &RandomSeed,
      CSTR("Seed for the random number generator.\n"
          SPACE "Default = 0.\n")),

  KW(CSTR("relativeoptimalitytolerance"),
      IlogCPDriver::set_cp_dbl_option, IloCP::RelativeOptimalityTolerance,
      CSTR("Relative tolerance on the objective value.\n"
          SPACE "Default = 1e-4.\n")),

  KW(CSTR("restartfaillimit"),
      IlogCPDriver::set_cp_int_option, &RestartFailLimit,
      CSTR("Number of failures allowed before restarting\n"
          SPACE "search.  Default = 100.\n")),

  KW(CSTR("restartgrowthfactor"),
      IlogCPDriver::set_cp_dbl_option, IloCP::RestartGrowthFactor,
      CSTR("Increase of the number of allowed failures\n"
          SPACE "before restarting search.  Default = 1.05.\n")),

  KW(CSTR("searchtype"), IlogCPDriver::set_cp_int_option, &SearchType,
      CSTR("Type of search used for solving a problem.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = depthfirst\n"
          SPACE "      1 = restart (default)\n"
          SPACE "      2 = multipoint\n")),

  KW(CSTR("solutionlimit"),
      IlogCPDriver::set_cp_int_option, &SolutionLimit,
      CSTR("Limit on the number of feasible solutions found\n"
          SPACE "before terminating a search. Default = no limit.\n")),

  KW(CSTR("temporalrelaxation"),
      IlogCPDriver::set_cp_int_option, &TemporalRelaxation,
      CSTR("0 or 1 (default 1):  Whether to use temporal\n"
          SPACE "relaxation.\n")),

  KW(CSTR("timelimit"),
      IlogCPDriver::set_cp_dbl_option, IloCP::TimeLimit,
      CSTR("Limit on the CPU time spent solving before\n"
          SPACE "terminating a search.  Default = no limit.\n")),

  KW(CSTR("timemode"),
      IlogCPDriver::set_cp_int_option, &TimeMode,
      CSTR("Specifies how the time is measured in CP\n"
          SPACE "Optimizer.  Possible values:\n"
          SPACE "      0 = cputime (default)\n"
          SPACE "      1 = elapsedtime\n")),

  KW(CSTR("timing"), IlogCPDriver::set_bool_option, IlogCPDriver::TIMING,
      CSTR("0 or 1 (default 0):  Whether to display timings\n"
          SPACE "for the run.\n")),

  KW(CSTR("usenumberof"), IlogCPDriver::set_bool_option, IlogCPDriver::USENUMBEROF,
      CSTR("0 or 1 (default 1):  Whether to consolidate\n"
          SPACE "'numberof' expressions by use of IloDistribute\n"
          SPACE "constraints.\n")),

  KW(CSTR("version"), Ver_val, 0,
      CSTR("Single-word phrase:  report version details\n"
          SPACE "before solving the problem.\n")),

  KW(CSTR("wantsol"), WS_val, 0,
      CSTR("In a stand-alone invocation (no -AMPL on the\n"
          SPACE "command line), what solution information to\n"
          SPACE "write.  Sum of\n"
          SPACE "      1 = write .sol file\n"
          SPACE "      2 = primal variables to stdout\n"
          SPACE "      4 = dual variables to stdout\n"
          SPACE "      8 = suppress solution message\n")),

  KW(CSTR("workers"),
      IlogCPDriver::set_cp_int_option, &Workers,
      CSTR("Number of workers to run in parallel to solve a\n"
          SPACE "problem.  In addition to numeric values this\n"
          SPACE "option accepts the value \"auto\" since CP\n"
          SPACE "Optimizer version 12.3.  Default = 1.\n"))
};

IlogCPDriver::IlogCPDriver() :
   mod_(env_), asl(Driver::asl()), gotopttype(false), n_badvals(0) {
  char *s;
  int n;
  size_t L;

  options_[DEBUGEXPR] = 0;
  options_[OPTIMIZER] = AUTO;
  options_[TIMING] = 0;
  options_[USENUMBEROF] = 1;

  version_.resize(L = strlen(IloConcertVersion::_ILO_NAME) + 100);
  n = snprintf(s = &version_[0], L, "AMPL/IBM ILOG CP Optimizer [%s %d.%d.%d]",
      IloConcertVersion::_ILO_NAME, IloConcertVersion::_ILO_MAJOR_VERSION,
      IloConcertVersion::_ILO_MINOR_VERSION,
      IloConcertVersion::_ILO_TECH_VERSION);
  DriverOptionInfo *doi = 0;
  oinfo_.reset(doi = new DriverOptionInfo());
  oinfo_->sname = CSTR("ilogcp");
  snprintf(oinfo_->bsname = s + n + 1, L - n, "ilogcp %d.%d.%d",
      IloConcertVersion::_ILO_MAJOR_VERSION,
      IloConcertVersion::_ILO_MINOR_VERSION,
      IloConcertVersion::_ILO_TECH_VERSION);
  oinfo_->opname = CSTR("ilogcp_options");
  oinfo_->keywds = keywords_;
  oinfo_->n_keywds = sizeof(keywords_) / sizeof(*keywords_);
  oinfo_->version = &version_[0];
  oinfo_->driver_date = YYYYMMDD;
  doi->driver = this;
}

IlogCPDriver::~IlogCPDriver() {
  env_.end();
}

char *IlogCPDriver::set_optimizer(Option_Info *oi, keyword *kw, char *value) {
  IlogCPDriver *d = static_cast<DriverOptionInfo*>(oi)->driver;
  int opt = 0;
  char *end = skip_nonspace(value);
  size_t length = end - value;
  if (strncmp(value, "auto", length) == 0) {
    opt = AUTO;
  } else if (strncmp(value, "cp", length) == 0) {
    opt = CP;
  } else if (strncmp(value, "cplex", length) == 0) {
    opt = CPLEX;
  } else {
    ++d->n_badvals;
    cerr << "Invalid value " << string(value, end)
             << " for option " << kw->name << endl;
  }
  if (!d->gotopttype)
    d->options_[OPTIMIZER] = opt;
  return end;
}

char *IlogCPDriver::set_int_option(Option_Info *oi, keyword *kw, char *value) {
  IlogCPDriver *d = static_cast<DriverOptionInfo*>(oi)->driver;
  if (!d->gotopttype)
    return skip_nonspace(value);
  keyword thiskw(*kw);
  thiskw.info = d->options_ + reinterpret_cast<size_t>(kw->info);
  return I_val(oi, &thiskw, value);
}

char *IlogCPDriver::set_bool_option(Option_Info *oi, keyword *kw, char *value) {
  IlogCPDriver *d = static_cast<DriverOptionInfo*>(oi)->driver;
  if (!d->gotopttype)
    return skip_nonspace(value);
  keyword thiskw(*kw);
  int intval = 0;
  thiskw.info = &intval;
  char *result = I_val(oi, &thiskw, value);
  if (intval != 0 && intval != 1) {
    ++d->n_badvals;
    cerr << "Invalid value " << value
        << " for option " << kw->name << endl;
  } else d->options_[reinterpret_cast<size_t>(kw->info)] = intval;
  return result;
}

void IlogCPDriver::set_option(keyword *kw, int value) {
  try {
    optimizer_->set_option(kw->info, value);
  } catch (const IloException &) {
    cerr << "Invalid value " << value << " for option " << kw->name << endl;
    ++n_badvals;
  }
}

char *IlogCPDriver::set_cp_int_option(Option_Info *oi, keyword *kw, char *value) {
  IlogCPDriver *d = static_cast<DriverOptionInfo*>(oi)->driver;
  if (!d->gotopttype)
    return skip_nonspace(value);
  if (d->get_option(OPTIMIZER) != CP) {
    ++d->n_badvals;
    cerr << "Invalid option " << kw->name << " for CPLEX optimizer" << endl;
    return skip_nonspace(value);
  }
  const CPOptionInfo *info = static_cast<const CPOptionInfo*>(kw->info);
  char c = *value;
  if ((info->values || info->accepts_auto) &&
      !isdigit(c) && c != '+' && c != '-') {
    char *end = skip_nonspace(value);
    if (info->values) {
      // Search for a value in the list of known values.
      // Use linear search since the number of values is small.
      for (int i = 0; info->values[i]; ++i) {
        if (strncmp(value, info->values[i], end - value) == 0) {
          d->set_option(kw, i);
          return end;
        }
      }
    }
    if (info->accepts_auto && strncmp(value, "auto", end - value) == 0) {
      d->set_option(kw, IloCP::Auto);
      return end;
    }
    cerr << "Invalid value " << string(value, end)
             << " for option " << kw->name << endl;
    ++d->n_badvals;
    return end;
  }
  keyword thiskw(*kw);
  int intval = 0;
  thiskw.info = &intval;
  char *result = I_val(oi, &thiskw, value);
  d->set_option(kw, intval);
  return result;
}

char *IlogCPDriver::set_cp_dbl_option(Option_Info *oi, keyword *kw, char *value) {
  IlogCPDriver *d = static_cast<DriverOptionInfo*>(oi)->driver;
  if (!d->gotopttype)
    return skip_nonspace(value);
  if (d->get_option(OPTIMIZER) != CP) {
    ++d->n_badvals;
    cerr << "Invalid option " << kw->name << " for CPLEX optimizer" << endl;
    return skip_nonspace(value);
  }
  keyword thiskw(*kw);
  double dblval = 0;
  thiskw.info = &dblval;
  char *result = D_val(oi, &thiskw, value);
  try {
    d->optimizer_->set_option(kw->info, dblval);
  } catch (const IloException &) {
    cerr << "Invalid value " << dblval << " for option " << kw->name << endl;
    ++d->n_badvals;
  }
  return result;
}

char *IlogCPDriver::set_cplex_int_option(Option_Info *oi, keyword *kw, char *value) {
  IlogCPDriver *d = static_cast<DriverOptionInfo*>(oi)->driver;
  if (!d->gotopttype)
    return skip_nonspace(value);
  if (d->get_option(OPTIMIZER) != CPLEX) {
    ++d->n_badvals;
    cerr << "Invalid option " << kw->name << " for CP optimizer" << endl;
    return skip_nonspace(value);
  }
  keyword thiskw(*kw);
  int intval = 0;
  thiskw.info = &intval;
  char *result = I_val(oi, &thiskw, value);
  d->set_option(kw, intval);
  return result;
}

bool IlogCPDriver::parse_options(char **argv) {
  // Get optimizer type.
  gotopttype = false;
  oinfo_->option_echo &= ~ASL_OI_echo;
  if (getopts(argv, oinfo_.get()))
    return false;

  int &opt = options_[OPTIMIZER];
  if (opt == AUTO)
    opt = nlo + nlc + n_lcon == 0 ? CPLEX : CP;
  if (opt == CPLEX)
    optimizer_.reset(new CPLEXOptimizer(env_, asl));
  else optimizer_.reset(new CPOptimizer(env_, asl));

  // Parse remaining options.
  gotopttype = true;
  n_badvals = 0;
  oinfo_->option_echo |= ASL_OI_echo;
  if (getopts(argv, oinfo_.get()) || n_badvals != 0)
    return false;

  debug_ = get_option(DEBUGEXPR);
  return true;
}

bool IlogCPDriver::show_version() const {
  return (oinfo_->flags & ASL_OI_show_version) != 0;
}

int IlogCPDriver::wantsol() const {
  return oinfo_->wantsol;
}

IloNumExprArray IlogCPDriver::ConvertArgs(VarArgExpr e) {
  IloNumExprArray args(env_);
  for (VarArgExpr::iterator i = e.begin(); NumericExpr arg = *i; ++i)
    args.add(Visit(arg));
  return args;
}

IloExpr IlogCPDriver::VisitIf(IfExpr e) {
  IloConstraint condition(Visit(e.condition()));
  IloNumVar var(env_, -IloInfinity, IloInfinity);
  mod_.add(IloIfThen(env_, condition, var == Visit(e.true_expr())));
  mod_.add(IloIfThen(env_, !condition, var == Visit(e.false_expr())));
  return var;
}

IloExpr IlogCPDriver::VisitAtan2(BinaryExpr e) {
  IloNumExpr y(Visit(e.lhs())), x(Visit(e.rhs()));
  IloNumExpr atan(IloArcTan(y / x));
  IloNumVar result(env_, -IloInfinity, IloInfinity);
  mod_.add(IloIfThen(env_, x >= 0, result == atan));
  mod_.add(IloIfThen(env_, x <= 0 && y >= 0, result == atan + M_PI));
  mod_.add(IloIfThen(env_, x <= 0 && y <= 0, result == atan - M_PI));
  return result;
}

IloExpr IlogCPDriver::VisitSum(SumExpr e) {
  IloExpr sum(env_);
  for (SumExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += Visit(*i);
  return sum;
}

IloExpr IlogCPDriver::VisitRound(BinaryExpr e) {
  NumericConstant num = Cast<NumericConstant>(e.rhs());
  if (!num || num.value() != 0)
    throw UnsupportedExprError("round with nonzero second parameter");
  // Note that IloOplRound rounds half up.
  return IloOplRound(Visit(e.lhs()));
}

IloExpr IlogCPDriver::VisitTrunc(BinaryExpr e) {
  NumericConstant num = Cast<NumericConstant>(e.rhs());
  if (!num || num.value() != 0)
    throw UnsupportedExprError("trunc with nonzero second parameter");
  return IloTrunc(Visit(e.lhs()));
}

IloExpr IlogCPDriver::VisitCount(CountExpr e) {
  IloExpr sum(env_);
  for (CountExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    sum += Visit(*i);
  return sum;
}

IloExpr IlogCPDriver::VisitNumberOf(NumberOfExpr e) {
  NumericExpr target = e.target();
  NumericConstant num = Cast<NumericConstant>(target);
  if (!num || !get_option(USENUMBEROF)) {
    IloExpr sum(env_);
    IloExpr concert_target(Visit(target));
    for (NumberOfExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      sum += (Visit(*i) == concert_target);
    return sum;
  }

  // If the first operand is constant, add it to the driver's data structure
  // that collects these operators.

  // Did we previously see a number-of operator
  // having the same expression-list?
  vector<NumberOf>::reverse_iterator np =
      find_if(numberofs_.rbegin(), numberofs_.rend(), SameExpr(e));

  // New expression-list:
  // Build a new numberof structure.
  if (np == numberofs_.rend()) {
    IloIntArray values(env_);
    values.add(num.value());

    IloIntVarArray vars(env_);
    for (NumberOfExpr::iterator i = e.begin(), end = e.end(); i != end; ++i) {
      IloIntVar var(env_, IloIntMin, IloIntMax);
      vars.add(var);
      mod_.add(var == Visit(*i));
    }

    IloIntVar cardVar(env_, IloIntMin, IloIntMax);
    IloIntVarArray cards(env_);
    cards.add(cardVar);
    numberofs_.push_back(NumberOf(cards, values, vars, e));
    return cardVar;
  }

  // Previously seen expression-list:
  // Add to its numberof structure.
  return np->Add(num.value(), env_);
}

IloExpr IlogCPDriver::VisitPLTerm(PiecewiseLinearTerm t) {
  IloNumArray slopes(env_), breakpoints(env_);
  int num_breakpoints = t.num_breakpoints();
  for (int i = 0; i < num_breakpoints; ++i) {
    slopes.add(t.slope(i));
    breakpoints.add(t.breakpoint(i));
  }
  slopes.add(t.slope(num_breakpoints));
  return IloPiecewiseLinear(vars_[t.var_index()], breakpoints, slopes, 0, 0);
}

IloConstraint IlogCPDriver::VisitExists(IteratedLogicalExpr e) {
  IloOr disjunction(env_);
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i) {
    disjunction.add(Visit(*i));
  }
  return disjunction;
}

IloConstraint IlogCPDriver::VisitForAll(IteratedLogicalExpr e) {
  IloAnd conjunction(env_);
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i) {
    conjunction.add(Visit(*i));
  }
  return conjunction;
}

IloConstraint IlogCPDriver::VisitImplication(ImplicationExpr e) {
  IloConstraint condition(Visit(e.condition()));
  return IloIfThen(env_,  condition, Visit(e.true_expr())) &&
      IloIfThen(env_, !condition, Visit(e.false_expr()));
}

IloConstraint IlogCPDriver::VisitAllDiff(AllDiffExpr e) {
  IloIntVarArray vars(env_);
  for (AllDiffExpr::iterator i = e.begin(), end = e.end(); i != end; ++i) {
    if (Variable var = Cast<Variable>(*i)) {
      vars.add(vars_[var.index()]);
      continue;
    }
    IloIntVar var(env_, IloIntMin, IloIntMax);
    mod_.add(var == Visit(*i));
    vars.add(var);
  }
  return IloAllDiff(env_, vars);
}

void IlogCPDriver::FinishBuildingNumberOf() {
  for (vector<NumberOf>::const_iterator
      i = numberofs_.begin(), end = numberofs_.end(); i != end; ++i) {
    mod_.add(i->Convert(env_));
  }
  numberofs_.clear();
}

/*----------------------------------------------------------------------

  Main Program

----------------------------------------------------------------------*/

int IlogCPDriver::run(char **argv) {
  /*** Initialize timers ***/

  double Times[5];
  Times[0] = xectim_();

  /*** Get name of .nl file; read problem sizes ***/

  char *stub = getstub(&argv, oinfo_.get());
  if (!stub) {
    usage_noexit_ASL(oinfo_.get(), 1);
    return 1;
  }
  FILE *nl = jac0dim(stub, strlen(stub));

  /*** Read coefficients & bounds & expression tree from .nl file ***/

  Uvx = static_cast<real*>(Malloc(n_var * sizeof(real)));
  Urhsx = static_cast<real*>(Malloc(n_con * sizeof(real)));

  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; i++)
    r_ops_int[i] = reinterpret_cast<efunc*>(i);
  asl->I.r_ops_ = r_ops_int;
  want_derivs = 0;
  fg_read(nl, ASL_allow_CLP);
  asl->I.r_ops_ = 0;

  if (!parse_options(argv))
    return 1;

  /*-------------------------------------------------------------------

     Set up optimization problem in ILOG Concert

   -------------------------------------------------------------------*/

  vars_ = optimizer_->vars();

  int n_var_int = nbv + niv + nlvbi + nlvci + nlvoi;
  int n_var_cont = n_var - n_var_int;
  if (n_var_cont != 0 && get_option(OPTIMIZER) == CP) {
    cerr << "CP Optimizer doesn't support continuous variables" << endl;
    return 1;
  }
  for (int j = 0; j < n_var_cont; j++)
    vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOFLOAT);
  for (int j = n_var_cont; j < n_var; j++)
    vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOINT);

  if (n_obj > 0) {
    NumericExpr expr(GetNonlinearObjExpr(0));
    NumericConstant constant(Cast<NumericConstant>(expr));
    IloExpr ilo_expr(env_, constant ? constant.value() : 0);
    if (0 < nlo)
      ilo_expr += Visit(expr);
    for (ograd *og = Ograd[0]; og; og = og->next)
      ilo_expr += og->coef * vars_[og->varno];
    IloObjective MinOrMax(env_, ilo_expr,
        objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
    optimizer_->set_obj(MinOrMax);
    IloAdd(mod_, MinOrMax);
  }

  IloRangeArray Con(optimizer_->cons());

  for (int i = 0; i < n_con; i++) {
    IloExpr conExpr(env_);
    for (cgrad *cg = Cgrad[i]; cg; cg = cg->next)
      conExpr += (cg -> coef) * vars_[cg -> varno];
    if (i < nlc)
      conExpr += Visit(GetNonlinearConExpr(i));
    Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
  }

  IloConstraintArray LCon(env_, n_lcon);

  for (int i = 0; i < n_lcon; i++)
    LCon[i] = Visit(GetLogicalConExpr(i));

  if (n_con > 0) mod_.add(Con);
  if (n_lcon > 0) mod_.add(LCon);

  FinishBuildingNumberOf();

  int timing = get_option(TIMING);
  Times[1] = timing ? xectim_() : 0;

  // Solve the problem.
  IloAlgorithm alg(optimizer_->algorithm());
  alg.extract (mod_);
  Times[2] = timing ? xectim_() : 0;
  IloBool successful = alg.solve();
  Times[3] = timing ? xectim_() : 0;

  // Convert solution status.
  const char *message;
  switch (alg.getStatus()) {
  default:
    // Fall through.
  case IloAlgorithm::Unknown:
    solve_result_num = 501;
    message = "unknown solution status";
    break;
  case IloAlgorithm::Feasible:
    solve_result_num = 100;
    message = "feasible solution";
    break;
  case IloAlgorithm::Optimal:
    solve_result_num = 0;
    message = "optimal solution";
    break;
  case IloAlgorithm::Infeasible:
    solve_result_num = 200;
    message = "infeasible problem";
    break;
  case IloAlgorithm::Unbounded:
    solve_result_num = 300;
    message = "unbounded problem";
    break;
  case IloAlgorithm::InfeasibleOrUnbounded:
    solve_result_num = 201;
    message = "infeasible or unbounded problem";
    break;
  case IloAlgorithm::Error:
    solve_result_num = 500;
    message = "error";
    break;
  }

  char sMsg[256];
  int sSoFar = Sprintf(sMsg, "%s: %s\n", oinfo_->bsname, message);
  vector<real> primal, dual;
  if (successful)
    optimizer_->get_solution(asl, sMsg + sSoFar, primal, dual);
  write_sol(sMsg, primal.empty() ? 0 : &primal[0],
      dual.empty() ? 0 : &dual[0], oinfo_.get());

  if (timing) {
    Times[4] = xectim_();
    cerr << endl
        << "Define = " << Times[1] - Times[0] << endl
        << "Setup =  " << Times[2] - Times[1] << endl
        << "Solve =  " << Times[3] - Times[2] << endl
        << "Output = " << Times[4] - Times[3] << endl;
  }
  return 0;
}
}

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

static char xxxvers[] = "ilogcp_options\0\n"
	"AMPL/IBM ILOG CP Optimizer Driver Version " qYYYYMMDD "\n";

// for suppressing "String literal to char*" warnings
#define CSTR(s) const_cast<char*>(s)

namespace {
struct DriverOptionInfo : Option_Info {
  Driver *driver;
};

// Returns the constant term in the first objective.
real objconst0(ASL_fg *a) {
  expr *e = a->I.obj_de_->e;
  return reinterpret_cast<size_t>(e->op) == OPNUM ?
      reinterpret_cast<expr_n*>(e)->v : 0;
}

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
}

Optimizer::Optimizer(IloEnv env, ASL_fg *asl) :
  vars_(env, n_var), cons_(env, n_con) {}

Optimizer::~Optimizer() {}

void CPLEXOptimizer::set_option(const void *key, int value) {
  cplex_.setParam(
      static_cast<IloCplex::IntParam>(reinterpret_cast<size_t>(key)), value);
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
keyword Driver::keywords_[] = {
  // The options must be in alphabetical order.

  KW(CSTR("alldiffinferencelevel"),
      Driver::set_cp_int_option, &AllDiffInferenceLevel,
      CSTR("Inference level for 'alldiff' constraints.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = default\n"
          SPACE "      1 = low\n"
          SPACE "      2 = basic\n"
          SPACE "      3 = medium\n"
          SPACE "      4 = extended\n")),

  KW(CSTR("branchlimit"),
      Driver::set_cp_int_option, &BranchLimit,
      CSTR("Limit on the number of branches made before\n"
          SPACE "terminating a search.  Default = no limit.\n")),

  KW(CSTR("choicepointlimit"),
      Driver::set_cp_int_option, &ChoicePointLimit,
      CSTR("Limit on the number of choice points created\n"
          SPACE "before terminating a search. Default = no limit.\n")),

  KW(CSTR("constraintaggregation"),
      Driver::set_cp_int_option, &ConstraintAggregation,
      CSTR("0 or 1 (default 1):  Whether to aggregate basic\n"
          SPACE "constraints.\n")),

  KW(CSTR("debugexpr"), Driver::set_int_option, Driver::DEBUGEXPR,
      CSTR("0 or 1 (default 0):  Whether to print debugging\n"
          SPACE "information for expression trees.\n")),

  KW(CSTR("defaultinferencelevel"),
      Driver::set_cp_int_option, &DefaultInferenceLevel,
      CSTR("Default inference level for constraints.\n"
          SPACE "Possible values:\n"
          SPACE "      1 = low\n"
          SPACE "      2 = basic\n"
          SPACE "      3 = medium\n"
          SPACE "      4 = extended\n")),

  KW(CSTR("distributeinferencelevel"),
      Driver::set_cp_int_option, &DistributeInferenceLevel,
      CSTR("Inference level for 'distribute' constraints.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = default\n"
          SPACE "      1 = low\n"
          SPACE "      2 = basic\n"
          SPACE "      3 = medium\n"
          SPACE "      4 = extended\n")),

  KW(CSTR("dynamicprobing"),
      Driver::set_cp_int_option, &DynamicProbing,
      CSTR("Use probing during search.  Possible values:\n"
          SPACE "     -1 = auto (default)\n"
          SPACE "      0 = off\n"
          SPACE "      1 = on\n")),

  KW(CSTR("dynamicprobingstrength"),
      Driver::set_cp_dbl_option, IloCP::DynamicProbingStrength,
      CSTR("Effort dedicated to dynamic probing as a factor\n"
          SPACE "of the total search effort.  Default = 0.03.\n")),

  KW(CSTR("faillimit"),
      Driver::set_cp_int_option, &FailLimit,
      CSTR("Limit on the number of failures allowed before\n"
          SPACE "terminating a search.  Default = no limit.\n")),

  KW(CSTR("ilogcplex"), Driver::use_cplex, Driver::ILOGOPTTYPE,
      CSTR("Single-word phrase:  use IBM ILOG CPLEX\n"
          SPACE "optimizer.\n")),

  KW(CSTR("ilogsolver"), Driver::use_cpoptimizer, Driver::ILOGOPTTYPE,
      CSTR("Single-word phrase:  use IBM ILOG Constraint\n"
          SPACE "Programming optimizer (default).\n")),

  KW(CSTR("logperiod"),
      Driver::set_cp_int_option, &LogPeriod,
      CSTR("Specifies how often the information in the\n"
          SPACE "search log is displayed.\n")),

  KW(CSTR("logverbosity"), Driver::set_cp_int_option, &LogVerbosity,
      CSTR("Verbosity of the search log.  Possible values:\n"
          SPACE "      0 = quiet\n"
          SPACE "      1 = terse\n"
          SPACE "      2 = normal (default)\n"
          SPACE "      3 = verbose\n")),

  KW(CSTR("multipointnumberofsearchpoints"),
      Driver::set_cp_int_option, &MultiPointNumberOfSearchPoints,
      CSTR("Number of solutions for the multi-point search\n"
          SPACE "algorithm.  Default = 30.\n")),

  KW(CSTR("optimalitytolerance"),
      Driver::set_cp_dbl_option, IloCP::OptimalityTolerance,
      CSTR("Absolute tolerance on the objective value.\n"
          SPACE "Default = 0.\n")),

  KW(CSTR("outlev"), Driver::set_cp_int_option,
      &LogVerbosity, CSTR("Synonym for \"logverbosity\".\n")),

  KW(CSTR("propagationlog"), Driver::set_cp_int_option, &PropagationLog,
      CSTR("Level of propagation trace reporting.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = quiet (default)\n"
          SPACE "      1 = terse\n"
          SPACE "      2 = normal\n"
          SPACE "      3 = verbose\n")),

  KW(CSTR("randomseed"), Driver::set_cp_int_option, &RandomSeed,
      CSTR("Seed for the random number generator.\n"
          SPACE "Default = 0.\n")),

  KW(CSTR("relativeoptimalitytolerance"),
      Driver::set_cp_dbl_option, IloCP::RelativeOptimalityTolerance,
      CSTR("Relative tolerance on the objective value.\n"
          SPACE "Default = 1e-4.\n")),

  KW(CSTR("restartfaillimit"),
      Driver::set_cp_int_option, &RestartFailLimit,
      CSTR("Number of failures allowed before restarting\n"
          SPACE "search.  Default = 100.\n")),

  KW(CSTR("restartgrowthfactor"),
      Driver::set_cp_dbl_option, IloCP::RestartGrowthFactor,
      CSTR("Increase of the number of allowed failures\n"
          SPACE "before restarting search.  Default = 1.05.\n")),

  KW(CSTR("searchtype"), Driver::set_cp_int_option, &SearchType,
      CSTR("Type of search used for solving a problem.\n"
          SPACE "Possible values:\n"
          SPACE "      0 = depthfirst\n"
          SPACE "      1 = restart (default)\n"
          SPACE "      2 = multipoint\n")),

  KW(CSTR("solutionlimit"),
      Driver::set_cp_int_option, &SolutionLimit,
      CSTR("Limit on the number of feasible solutions found\n"
          SPACE "before terminating a search. Default = no limit.\n")),

  KW(CSTR("temporalrelaxation"),
      Driver::set_cp_int_option, &TemporalRelaxation,
      CSTR("0 or 1 (default 1):  Whether to use temporal\n"
          SPACE "relaxation.\n")),

  KW(CSTR("timelimit"),
      Driver::set_cp_dbl_option, IloCP::TimeLimit,
      CSTR("Limit on the CPU time spent solving before\n"
          SPACE "terminating a search.  Default = no limit.\n")),

  KW(CSTR("timemode"),
      Driver::set_cp_int_option, &TimeMode,
      CSTR("Specifies how the time is measured in CP\n"
          SPACE "Optimizer.  Possible values:\n"
          SPACE "      0 = cputime (default)\n"
          SPACE "      1 = elapsedtime\n")),

  KW(CSTR("timing"), Driver::set_bool_option, Driver::TIMING,
      CSTR("0 or 1 (default 0):  Whether to display timings\n"
          SPACE "for the run.\n")),

  KW(CSTR("usenumberof"), Driver::set_bool_option, Driver::USENUMBEROF,
      CSTR("0 or 1 (default 1):  Whether to consolidate\n"
          SPACE "'numberof' expressions by use of IloDistribute\n"
          SPACE "constraints.\n")),

  KW(CSTR("wantsol"), WS_val, 0,
      CSTR("In a stand-alone invocation (no -AMPL on the\n"
          SPACE "command line), what solution information to\n"
          SPACE "write.  Sum of\n"
          SPACE "      1 = write .sol file\n"
          SPACE "      2 = primal variables to stdout\n"
          SPACE "      4 = dual variables to stdout\n"
          SPACE "      8 = suppress solution message\n")),

  KW(CSTR("workers"),
      Driver::set_cp_int_option, &Workers,
      CSTR("Number of workers to run in parallel to solve a\n"
          SPACE "problem.  In addition to numeric values this\n"
          SPACE "option accepts the value \"auto\" since CP\n"
          SPACE "Optimizer version 12.3.  Default = 1.\n"))
};

Driver::Driver() :
   mod_(env_), asl(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))),
   gotopttype(false), n_badvals(0) {
	char *s;
	int n;
	size_t L;

	options_[DEBUGEXPR] = 0;
	options_[ILOGOPTTYPE] = DEFAULT_OPT;
	options_[TIMING] = 0;
	options_[USENUMBEROF] = 1;

	version_.resize(L = strlen(IloConcertVersion::_ILO_NAME) + 100);
	n = snprintf(s = &version_[0], L,
		"AMPL/IBM ILOG CP Optimizer [%s %d.%d.%d]",
		IloConcertVersion::_ILO_NAME,
		IloConcertVersion::_ILO_MAJOR_VERSION,
		IloConcertVersion::_ILO_MINOR_VERSION,
		IloConcertVersion::_ILO_TECH_VERSION);
	DriverOptionInfo *doi = 0;
	oinfo_.reset(doi = new DriverOptionInfo());
	oinfo_->sname = CSTR("ilogcp");
	snprintf(oinfo_->bsname = s + n + 1, L-n, "ilogcp %d.%d.%d",
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

Driver::~Driver() {
   env_.end();
   ASL_free(reinterpret_cast<ASL**>(&asl));
}

char *Driver::use_cplex(Option_Info *oi, keyword *, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      d->options_[ILOGOPTTYPE] = CPLEX;
   return value;
}

char *Driver::use_cpoptimizer(Option_Info *oi, keyword *, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      d->options_[ILOGOPTTYPE] = CPOPTIMIZER;
   return value;
}

char *Driver::set_int_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      return skip_nonspace(value);
   keyword thiskw(*kw);
   thiskw.info = d->options_ + reinterpret_cast<size_t>(kw->info);
   return I_val(oi, &thiskw, value);
}

char *Driver::set_bool_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
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

void Driver::set_cp_option(keyword *kw, int value) {
  try {
     optimizer_->set_option(kw->info, value);
  } catch (const IloException &) {
     cerr << "Invalid value " << value << " for option " << kw->name << endl;
     ++n_badvals;
  }
}

char *Driver::set_cp_int_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      return skip_nonspace(value);
   if (d->get_option(ILOGOPTTYPE) != CPOPTIMIZER) {
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
           d->set_cp_option(kw, i);
           return end;
         }
       }
     }
     if (info->accepts_auto && strncmp(value, "auto", end - value) == 0) {
       d->set_cp_option(kw, IloCP::Auto);
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
   d->set_cp_option(kw, intval);
   return result;
}

char *Driver::set_cp_dbl_option(Option_Info *oi, keyword *kw, char *value) {
   Driver *d = static_cast<DriverOptionInfo*>(oi)->driver;
   if (!d->gotopttype)
      return skip_nonspace(value);
   if (d->get_option(ILOGOPTTYPE) != CPOPTIMIZER) {
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

bool Driver::parse_options(char **argv) {
   // Get optimizer type.
   gotopttype = false;
   oinfo_->option_echo &= ~ASL_OI_echo;
   if (getopts(argv, oinfo_.get()))
      return false;

   int &ilogopttype = options_[ILOGOPTTYPE];
   if (ilogopttype == DEFAULT_OPT)
      ilogopttype = nlo + nlc + n_lcon == 0 ? CPLEX : CPOPTIMIZER;
   if (ilogopttype == CPLEX)
      optimizer_.reset(new CPLEXOptimizer(env_, asl));
   else optimizer_.reset(new CPOptimizer(env_, asl));

   // Parse remaining options.
   gotopttype = true;
   n_badvals = 0;
   oinfo_->option_echo |= ASL_OI_echo;
   if (getopts(argv, oinfo_.get()) || n_badvals != 0)
      return false;
   return true;
}

int Driver::wantsol() const
{
   return oinfo_->wantsol;
}

/*----------------------------------------------------------------------

  Main Program

----------------------------------------------------------------------*/

int Driver::run(char **argv) {
   /*** Initialize timers ***/

   IloTimer timer(env_);
   timer.start();

   IloNum Times[5];
   Times[0] = timer.getTime();

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
   if (n_var_cont != 0 && get_option(ILOGOPTTYPE) == CPOPTIMIZER) {
      cerr << "CP Optimizer doesn't support continuous variables" << endl;
      return 1;
   }
   for (int j = 0; j < n_var_cont; j++)
      vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOFLOAT);
   for (int j = n_var_cont; j < n_var; j++)
      vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOINT);

   if (n_obj > 0) {
      IloExpr objExpr(env_, objconst0(asl));
      if (0 < nlo)
         objExpr += build_expr (obj_de[0].e);
      for (ograd *og = Ograd[0]; og; og = og->next)
         objExpr += (og -> coef) * vars_[og -> varno];
      IloObjective MinOrMax(env_, objExpr,
         objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
      optimizer_->set_obj(MinOrMax);
      IloAdd (mod_, MinOrMax);
   }

   IloRangeArray Con(optimizer_->cons());

   for (int i = 0; i < n_con; i++) {
      IloExpr conExpr(env_);
      for (cgrad *cg = Cgrad[i]; cg; cg = cg->next)
         conExpr += (cg -> coef) * vars_[cg -> varno];
      if (i < nlc)
         conExpr += build_expr (con_de[i].e);
      Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
   }

   IloConstraintArray LCon(env_,n_lcon);

   for (int i = 0; i < n_lcon; i++)
      LCon[i] = build_constr (lcon_de[i].e);

   if (n_con > 0) mod_.add (Con);
   if (n_lcon > 0) mod_.add (LCon);

   finish_building_numberof ();

   int timing = get_option(TIMING);
   Times[1] = timing ? timer.getTime() : 0;

   // Solve the problem.
   IloAlgorithm alg(optimizer_->algorithm());
   alg.extract (mod_);
   Times[2] = timing ? timer.getTime() : 0;
   IloBool successful = alg.solve();
   Times[3] = timing ? timer.getTime() : 0;

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
      Times[4] = timer.getTime();
      cerr << endl
           << "Define = " << Times[1] - Times[0] << endl
           << "Setup =  " << Times[2] - Times[1] << endl
           << "Solve =  " << Times[3] - Times[2] << endl
           << "Output = " << Times[4] - Times[3] << endl;
   }
   return 0;
}

#include <vector>
#include <climits>
#include <cfloat>
#include <algorithm>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "highsmpbackend.h"

extern "C" {
  #include "highsmp-ampls-c-api.h"    // Highs AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptHighs(void* prob) {
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateHighsBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()};
}


namespace mp {

/// Create Highs Model Manager
/// @param gc: the Highs common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return HighsModelMgr
std::unique_ptr<BasicModelManager>
CreateHighsModelMgr(HighsCommon&, Env&, pre::BasicValuePresolver*&);


HighsBackend::HighsBackend() {
  OpenSolver();

  pre::BasicValuePresolver* pPre;
  auto data = CreateHighsModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

HighsBackend::~HighsBackend() {
  CloseSolver();
}

const char* HighsBackend::GetBackendName()
  { return "HighsBackend"; }

std::string HighsBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", HIGHS_VERSION_MAJOR, 
    HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH);
}


bool HighsBackend::IsQCP() const {
  return false; 
}

ArrayRef<double> HighsBackend::PrimalSolution() {
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  int primal_solution_status;
  Highs_getIntInfoValue(lp(),
                        "primal_solution_status", &primal_solution_status);
  if (kHighsSolutionStatusFeasible == primal_solution_status)
    Highs_getSolution(lp(), x.data(), NULL, NULL, NULL);
  else
    x.clear();
  return x;
}

pre::ValueMapDbl HighsBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> HighsBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  int error = Highs_getSolution(lp(), NULL, NULL, NULL, pi.data());
  if (error)
    pi.clear();
  return pi;
}

double HighsBackend::ObjectiveValue() const {
  return Highs_getObjectiveValue(lp());
}

double HighsBackend::NodeCount() const {
  return getInt64Attr("mip_node_count");
}

double HighsBackend::SimplexIterations() const {
  return getIntAttr("simplex_iteration_count");
}

int HighsBackend::BarrierIterations() const {
  return getIntAttr("ipm_iteration_count");
}

void HighsBackend::DoWriteProblem(const std::string &file) {
  HIGHS_CCALL(Highs_writeModel(lp(), file.c_str()));
}

void HighsBackend::DoWriteSolution(const std::string &file) {
  HIGHS_CCALL(Highs_writeSolutionPretty(lp(), file.c_str()));
}


void HighsBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptHighs, lp());
}

void HighsBackend::Solve() {
  HIGHS_CCALL( Highs_run(lp()) );
  WindupHIGHSSolve();
}

void HighsBackend::WindupHIGHSSolve() { }

void HighsBackend::ReportResults() {
  ReportHIGHSResults();
  BaseBackend::ReportResults();
}

void HighsBackend::ReportHIGHSResults() {
  SetStatus( ConvertHIGHSStatus() );
  AddHIGHSMessages();
}

SolutionBasis HighsBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void HighsBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarConStatii(varstt, constt);
}


/// We also use this for MIP.
/// Attempting partial MIP start.
void HighsBackend::AddPrimalDualStart(Solution sol0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution(
    { sol0_unpres.primal, sol0_unpres.dual });
  auto ms = GetValuePresolver().PresolveGenericInt(
        { sol0_unpres.spars_primal } );
  auto x0 = mv.GetVarValues()();
  auto s0 = ms.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);
  /// If all variables provided or all missing are continuous,
  /// set warmstart, otherwise fix, solve, unfix, set warmstart
  if (s0.size() < (size_t)NumVars())
    s0.resize(NumVars());
  if (0==*std::min_element(s0.begin(), s0.end())) {
    bool fAllMissingAreRealVars = true;
    for (auto j=s0.size(); j--; ) {
      int integr;
      auto res = Highs_getColIntegrality(lp(), j, &integr);
      if (kHighsStatusOk != res)
        break;       // no information, it's an LP
      if (kHighsVarTypeContinuous != integr) {
        fAllMissingAreRealVars = false;
        break;
      }
    }
    if (!fAllMissingAreRealVars) {
      std::vector<double> costs(NumVars());
      std::vector<double> lb(NumVars());
      std::vector<double> ub(NumVars());
      int numnz, ncols;
      Highs_getColsByRange(lp(), 0, NumVars()-1, &ncols,
        costs.data(), lb.data(), ub.data(), &numnz, NULL, NULL, NULL);
      Highs_changeColsBoundsByMask(lp(), s0.data(), x0.data(), x0.data());
      Highs_run(lp());
      x0 = PrimalSolution();         // get new solution
      Highs_changeColsBoundsByMask(lp(), s0.data(), lb.data(), ub.data());
    }
  }
  HIGHS_CCALL(Highs_setSolution(lp(), x0.data(), NULL, NULL, pi0.data()));
}


ArrayRef<int> HighsBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  conStatiii_.resize(NumLinCons());
  HIGHS_CCALL(Highs_getBasis(lp(), vars.data(), conStatiii_.data()));
  for (auto& s : vars) {
    switch (s) {
    case kHighsBasisStatusBasic:
      s = (int)BasicStatus::bas;
      break;
    case kHighsBasisStatusLower:
      s = (int)BasicStatus::low;
      break;
    case kHighsBasisStatusUpper:
      s = (int)BasicStatus::upp;
      break;
    case kHighsBasisStatusNonbasic:
    case kHighsBasisStatusZero:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Highs VBasis value: {}", s));
    }
  }
  return vars;
}

ArrayRef<int> HighsBackend::ConStatii() {
  for (auto& s : conStatiii_) {
    switch (s) {
    case kHighsBasisStatusBasic:
      s = (int)BasicStatus::bas;
      break;
    case kHighsBasisStatusLower:
      s = (int)BasicStatus::low;
      break;
    case kHighsBasisStatusUpper:
      s = (int)BasicStatus::upp;
      break;
    case kHighsBasisStatusNonbasic:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Highs VBasis value: {}", s));
    }
  }
  return conStatiii_;
}

void HighsBackend::VarConStatii(ArrayRef<int> vst, ArrayRef<int> cst) {
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  std::vector<int> indicesOfMissing;
  for (size_t j = 0; j<stt.size(); j++) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = kHighsBasisStatusBasic;
      break;
    case BasicStatus::low:
    case BasicStatus::equ:
      s = kHighsBasisStatusLower;
      break;
    case BasicStatus::upp:
      s = kHighsBasisStatusUpper;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = kHighsBasisStatusNonbasic;
      break;
    case BasicStatus::none:
      indicesOfMissing.push_back(j);
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }

  if (indicesOfMissing.size() > 0)
  {
    /// 'none' is assigned to new variables. Compute low/upp/sup:
    /// Depending on where 0.0 is between bounds
    std::vector<double> lb(indicesOfMissing.size());
    std::vector<double> ub(indicesOfMissing.size());
    std::vector<int> di(indicesOfMissing.size(),   // dummy pointers
                        indicesOfMissing.size());
    std::vector<double> dd(indicesOfMissing.size());
    int numnz;
    Highs_getColsBySet(lp(), indicesOfMissing.size(), indicesOfMissing.data(),
      di.data(), dd.data(), lb.data(), ub.data(), &numnz, NULL, NULL, NULL);
    for (size_t i = 0; i < indicesOfMissing.size(); i++) {
        if (lb[i] >= -1e-6)
          stt[indicesOfMissing[i]] = kHighsBasisStatusLower;
        else if (ub[i] <= 1e-6)
          stt[indicesOfMissing[i]] = kHighsBasisStatusUpper;
        else
          stt[indicesOfMissing[i]] = kHighsBasisStatusNonbasic;
      }
  }
  std::vector<int> cstt(cst.data(), cst.data() + cst.size());
  for (auto& s : cstt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = kHighsBasisStatusBasic;
      break;
    case BasicStatus::none:
    case BasicStatus::upp:   
    case BasicStatus::sup:   
    case BasicStatus::low:    
    case BasicStatus::equ:    
    case BasicStatus::btw:    
      s = kHighsBasisStatusNonbasic;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  HIGHS_CCALL(Highs_setBasis(lp(), stt.data(), cstt.data()));
}

void HighsBackend::AddHIGHSMessages() {
  auto ni = SimplexIterations();
  if (true)
    AddToSolverMessage(
          fmt::format("{} simplex iterations\n", std::max(0.0, ni)));
  auto nbi = BarrierIterations();
  if (nbi > -1)
    AddToSolverMessage(
      fmt::format("{} barrier iterations\n", nbi));
  auto nnd = NodeCount();
  if (nnd > -1)
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> HighsBackend::ConvertHIGHSStatus() {
  namespace sol = mp::sol;
  int optstatus = Highs_getModelStatus(lp());
  auto obj = Highs_getObjectiveValue(lp());
  auto inf = Highs_getInfinity(lp());
  bool hasSol = (-inf < obj && obj < inf);
  switch (optstatus) {
  case kHighsModelStatusOptimal:
    return { sol::SOLVED, "optimal solution" };
  case kHighsModelStatusInfeasible:
    return { sol::INFEASIBLE, "infeasible problem" };
  case kHighsModelStatusUnbounded:
    if (hasSol)
      return { sol::UNBOUNDED_FEAS, "unbounded problem, feasible solution" };
    return { sol::UNBOUNDED_NO_FEAS, "unbounded problem, no solution" };
  case kHighsModelStatusUnboundedOrInfeasible:
    return { sol::LIMIT_INF_UNB, "unbounded or infeasible" };
  case kHighsModelStatusModelError:
  case kHighsModelStatusLoadError:
    return { sol::FAILURE, "solver error" };
  case kHighsModelStatusPresolveError:
  case kHighsModelStatusSolveError:
  case kHighsModelStatusPostsolveError:
    if (hasSol)
      return { sol::UNCERTAIN, "numeric issue, solution candidate returned" };
    return { sol::NUMERIC, "numeric issue" };
  case kHighsModelStatusTimeLimit:
    if (hasSol)
      return { sol::LIMIT_FEAS_TIME, "time limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_TIME, "time limit, no solution" };
  case kHighsModelStatusIterationLimit:
    if (hasSol)
      return { sol::LIMIT_FEAS_ITER, "iteration limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_ITER, "iteration limit, no solution" };
  case kHighsModelStatusSolutionLimit:
    assert (hasSol);
    return { sol::LIMIT_FEAS_NUMSOLS, "solution limit" };
  case kHighsModelStatusInterrupt:
    if (hasSol)
      return { sol::LIMIT_FEAS_INTERRUPT, "interrupt, feasible solution" };
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupt, no solution" };
  case kHighsModelStatusObjectiveBound:
    if (hasSol)
      return { sol::LIMIT_FEAS_BESTOBJ, "objective bound, feasible solution" };
    return { sol::LIMIT_NO_FEAS_BESTBND, "objective bound, no solution" };
  case kHighsModelStatusObjectiveTarget:
    if (hasSol)
      return { sol::LIMIT_FEAS_BESTOBJ, "objective target, feasible solution" };
    return { sol::LIMIT_NO_FEAS_CUTOFF, "objective target, no solution" };
  default:
    if (hasSol)
      return { sol::UNCERTAIN, "unknown, solution candidate returned" };
    return { sol::UNKNOWN, "unknown" };
  }
  return { sol::UNKNOWN, "not solved" };
}


void HighsBackend::FinishOptionParsing() {
  int v=-1;
  GetSolverOption("output_flag", v);
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo lp_values_method[] = {
  { "choose", "Automatic (default)", -1},
  { "simplex", "Simplex", 1},
  { "ipm", "Interior Point Method", 2},
};

static const mp::OptionValueInfo off_on_choose_values[] = {
  { "choose", "Automatic (default)", -1},
  { "off", "Off", 1},
  { "on", "On", 2},
};

static const mp::OptionValueInfo run_crossover_values[] = {
  { "choose", "Run if the results of IPM without crossover is imprecise", -1},
  { "off", "Off", 1},
  { "on", "On", 2},
};
static const mp::OptionValueInfo simplex_strategy_values_[] = {
  { "0", "Choose automatically (default)", 0},
  { "1", "Dual (serial)", 1},
  { "2", "Dual ('PAMI' - Parallelization Across Multiple Iterations)", 2},
  { "3", "Dual ('SIP' - Single Iteration Parallelism", 3},
  { "4", "Primal", 4}
};
static const mp::OptionValueInfo simplex_scale_strategy_values_[] = {
  { "0", "Off", 0},
  { "1", "Choose automatically (default)", 1},
  { "2", "Equilibration", 2},
  { "3", "Forced equilibration", 3},
  { "4", "Max value 0", 4},
  { "5", "Max value 1", 5}
};
static const mp::OptionValueInfo simplex_crash_strategy_values_[] = {
  { "0", "Off (default)", 0},
  { "1", "LTSSF", 1},
  { "2", "Bixby", 2}
};
static const mp::OptionValueInfo simplex_edge_weight_strategy_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Dantzig", 0},
  { "1", "Devex", 1},
  { "2", "Steepest", 2}
};

void HighsBackend::InitCustomOptions() {

  set_option_header(
      "HIGHS Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``highs_options``. For example::\n"
      "\n"
      "  ampl: option highs_options 'relgaptol=1e-6';\n");
 
  AddSolverOption("tech:outlev outlev",
    "0*/1: Whether to write HighS log lines (chatter) to stdout and to file.",
    "output_flag", 0, 1);

  AddSolverOption("tech:logfile logfile",
    "Log file name.", "log_file");

  const char* c = "choose";
  AddSolverOption("alg:method method lpmethod",
    "Which algorithm to use :\n"
    "\n.. value-table::\n", "solver", lp_values_method, c);

  AddSolverOption("alg:simplex simplex simplex_strategy",
    "Strategy for simplex solver :\n"
    "\n.. value-table::\n", "simplex_strategy", simplex_strategy_values_, 0);

  AddSolverOption("alg:simplexscale simplexscale simplex_scale_strategy",
    "Simplex scaling strategy :\n"
    "\n.. value-table::\n", "simplex_scale_strategy", 
    simplex_scale_strategy_values_, 1);

  AddSolverOption("alg:simplexcrash simplexcrash simplex_crash_strategy",
    "Simplex crash strategy :\n"
    "\n.. value-table::\n", "simplex_crash_strategy",
    simplex_crash_strategy_values_, 0);

  AddSolverOption("alg:simplexdualedge simplexdualedge simplex_dual_edge_weight_strategy",
    "Simplex dual edge weights strategy :\n"
    "\n.. value-table::\n", "simplex_dual_edge_weight_strategy",
    simplex_edge_weight_strategy_values_, 1);

  AddSolverOption("alg:simplexprimaledge simplexprimaledge simplex_primal_edge_weight_strategy",
    "Simplex primal edge weights strategy :\n"
    "\n.. value-table::\n", "simplex_primal_edge_weight_strategy",
    simplex_edge_weight_strategy_values_, 1);

  AddSolverOption("pre:solve presolve",
    "Whether to use presolve:\n"
    "\n.. value-table::\n",
    "presolve", off_on_choose_values, c);

  AddSolverOption("alg:parallel parallel",
    "Parallel option :\n"
    "\n.. value-table::\n", "parallel", off_on_choose_values, c);

  AddSolverOption("lim:time timelim timelimit time_limit",
    "Limit on solve time (in seconds; default: no limit).",
    "time_limit", 0.0, DBL_MAX);

  AddSolverOption("lim:simplexiterationlimit simplexiterationlimit simplex_iteration_limit",
    "Limit on simplex iterations (default: no limit).",
    "simplex_iteration_limit", 0, INT_MAX);

  AddSolverOption("lim:ipmiterationlimit ipmiterationlimit ipm_iteration_limit",
    "Limit on IPM iterations (default: no limit).",
    "ipm_iteration_limit", 0, INT_MAX);

  AddSolverOption("alg:infinitecost infinitecost infinite_cost",
    "Limit on cost coefficient : values larger than this will be treated as infinite (default: 1e20).",
    "infinite_cost", 1e15, Infinity());

  AddSolverOption("alg:infinitebound infinitebound infinite_bound",
    "Limit on |constraint bound|: values larger than this will be treated as infinite (default: 1e20).",
    "infinite_cost", 1e15, Infinity());

  AddSolverOption("alg:infinitecoeff infinitecoeff large_matrix_value",
    "Upper limit on |matrix entries|: values larger than this will be treated as infinite (default: 1e15).",
    "large_matrix_value", 1.0, Infinity());

  AddSolverOption("alg:zerocoeff zerocoeff small_matrix_value",
    "Lower limit on |matrix entries|: values smaller than this will be treated as zero (default: 1e-9).",
    "small_matrix_value", 1e-12, Infinity());

  AddSolverOption("alg:feastol feastol primal_feasibility_tolerance",
    "Primal feasibility tolerance (default 1e-7).",
    "primal_feasibility_tolerance", 1e-10, Infinity());

  AddSolverOption("alg:dualfeastol dualfeastol dual_feasibility_tolerance",
    "Dual feasibility tolerance (default 1e-7).",
    "dual_feasibility_tolerance", 1e-10, Infinity());

  AddSolverOption("alg:ipmopttol ipmopttol ipm_optimality_tolerance",
    "IPM optimality tolerance (default 1e-8).",
    "ipm_optimality_tolerance", 1e-12, Infinity());

  AddSolverOption("bar:crossover run_crossover",
    "Run crossover after IPM to get a basic solution",
    "run_crossover", run_crossover_values, c);

  AddSolverOption("tech:threads threads",
    "How many threads to use when using the barrier algorithm "
    "or solving MIP problems; default 0 ==> automatic choice.",
		"threads", 0, INT32_MAX);

  AddSolverOption("mip:detsimmetry detsimmetry mip_detect_symmetry",
    "Whether symmetry should be detected (default 1)",
    "mip_detect_symmetry",0, 1);

  AddSolverOption("lim:stallnodes stallnodelim stallnodelimit mip_max_stall_nodes",
    "Maximum MIP number of nodes where estimate is above cutoff bound (default: no limit).",
    "mip_max_stall_nodes", 0, INT_MAX);

  AddSolverOption("lim:leavenodes leaveslim mip_max_leaves",
    "Maximum MIP number of leaf nodes (default: no limit).",
    "mip_max_leaves", 0, INT_MAX);

  AddSolverOption("lim:nodes nodelim nodelimit mip_max_nodes",
    "Maximum MIP nodes to explore (default: no limit).",
    "mip_max_nodes", 0, INT_MAX);

  AddSolverOption("lim:improvingsols improvingsolslimit mip_max_improving_sols",
    "Maximum number of improving solutions found (default: no limit).",
    "mip_max_improving_sols", 1, INT_MAX);

  AddSolverOption("mip:lpagelimit lpagelimit mip_lp_age_limit",
    "Maximal age of dynamic LP rows before they are removed from the LP relaxation "
    "(default 10)",
    "mip_lp_age_limit", 0, INT_MAX);

  AddSolverOption("mip:poolsoftlimit poolsoftlimit mip_pool_soft_limit",
    "Soft limit on the number of rows in the cutpool for dynamic age adjustment"
    "(default 10000)",
    "mip_pool_soft_limit", 1, INT_MAX);

  AddSolverOption("mip:pscostreliability pscostreliability mip_pscost_minreliable",
    "Minimal number of observations before pseudo costs are considered reliable"
    "(default 8)",
    "mip_pscost_minreliable", 0, INT_MAX);

  AddSolverOption("mip:mincliquetable mincliquetable mip_min_cliquetable_entries_for_parallelism",
    "Minimal number of entries in the cliquetable before neighborhood queries of the conflict graph use parallel processing"
    "(default 100000)",
    "mip_min_cliquetable_entries_for_parallelism", 0, INT_MAX);

  AddSolverOption("tech:miploglev miploglev mip_report_level",
    "0/1*/2: MIP solver report level",
    "mip_report_level", 0, 2);

  AddSolverOption("mip:intfeastol intfeastol inttol mip_feasibility_tolerance",
    "Feasibility tolerance for integer variables (default 1e-06).",
    "mip_feasibility_tolerance", 1e-10, Infinity());

  AddSolverOption("mip:heureff heureff mip_heuristic_effort",
    "Fraction of time to spend in MIP heuristics (default 0.05).",
    "mip_heuristic_effort", 0.0, 1.0);

  AddSolverOption("mip:relgaptol relgaptol mip_rel_gap",
    "Tolerance on relative gap, | ub - lb|/|ub | , to determine whether optimality has been reached for a MIP instance "
    "(default 1e-04).",
    "mip_rel_gap", 0.0, Infinity());

  AddSolverOption("mip:absgaptol absgaptol mip_abs_gap",
    "Tolerance on absolute gap of MIP, |ub-lb|, to determine whether optimality has been reached for a MIP instance "
    "(default 1e-06).",
    "mip_abs_gap", 0.0, Infinity());
}

double HighsBackend::MIPGap() {
  // TODO Check if the following is always true
  if (BarrierIterations() == 0)
    return 0;
  return getDblAttr("mip_gap");
}
double HighsBackend::BestDualBound() {
  // TODO Check if the following is always true
  if (BarrierIterations() == 0)
    return 0;
  return getDblAttr("mip_dual_bound");
}

double HighsBackend::MIPGapAbs() {
  // TODO Check if the following is always true
  if (BarrierIterations() == 0)
    return 0;
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}
} // namespace mp


// AMPLs

AMPLS_MP_Solver* AMPLSOpenHighs(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()},
    cb);
}

void AMPLSCloseHighs(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetHighsmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::HighsBackend*>(AMPLSGetBackend(slv))->lp();
}

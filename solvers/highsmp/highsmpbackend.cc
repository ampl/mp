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
  LoadHighsLibrary(false);
  OpenSolver(); 
  pre::BasicValuePresolver* pPre;
  auto data = CreateHighsModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
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

  void HighsBackend::InputExtras() {
  BaseBackend::InputExtras();
  // Set the accumulated objectives
  // In case of native MO, accObjectives will contain all the problem's objectives - and their meta into
  // and will have the flag hadMultiObjective set to true
  accObjectives().setAllInHighs(lp());
  accObjectives().clear(); 
}
  void HighsBackend::FinishOptionParsing() {
    // If CUDA is specified, scrap the current library, load the cuda one
    // and replay all the options

    if(storedOptions_.lpmethod_== "pdlp-gpu")
    //if (storedOptions_.useGPU_)
    {
      #ifdef __APPLE__
            throw std::runtime_error("GPU support is not available on MacOS");
      #endif
      loader().Highs_destroy(lp());
      LoadHighsLibrary(true);
      OpenSolver();
      ReplaySolverOptions();
    }
    std::string method = storedOptions_.lpmethod_ == "pdlp-gpu" ? "pdlp" : storedOptions_.lpmethod_;
    SetSolverOption("solver", method);
    /// Copy env/lp to ModelAPI
    copy_common_info_to_other();
    int v = -1;
    GetSolverOption("output_flag", v);
    set_verbose_mode(v > 0);
  }



bool HighsBackend::IsQCP() const {
  return false; 
}

ArrayRef<double> HighsBackend::PrimalSolution() {
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  int primal_solution_status;
  loader().Highs_getIntInfoValue(lp(),
                        "primal_solution_status", &primal_solution_status);
  if (kHighsSolutionStatusFeasible == primal_solution_status)
    loader().Highs_getSolution(lp(), x.data(), NULL, NULL, NULL);
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
  int error = loader().Highs_getSolution(lp(), NULL, NULL, NULL, pi.data());
  if (error)
    pi.clear();
  return pi;
}

  ArrayRef<double> HighsBackend::GetObjectiveValues() {
    if (accObjectives().hadNativeMultiObj()) {
      AddToSolverMessage("Warning: HiGHS does not support returning the objective values for multiple\n"
        "objectives; AMPL will compute the actual values.");
      return std::vector<double>();
    }
    return std::vector<double> { ObjectiveValue() };
}

  void HighsBackend::ObjPriorities(ArrayRef<int> pri) {
    if (pri.size() > 0)
      SetSolverOption("blend_multi_objectives",
        0);
    accObjectives().setPriorities(pri);

}
  void HighsBackend::ObjWeights(ArrayRef<double> w) {
    accObjectives().setWeights(w);
}
  void HighsBackend::ObjAbsTol(ArrayRef<double> a) {
    accObjectives().setAbsTols(a);
  }
  void HighsBackend::ObjRelTol(ArrayRef<double> r) {
    accObjectives().setRelTols(r);
  }

double HighsBackend::ObjectiveValue() const {
  return loader().Highs_getObjectiveValue(lp());
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
int HighsBackend::PdlpIterations() const {
  return getIntAttr("pdlp_iteration_count");
}

void HighsBackend::DoWriteProblem(const std::string &file) {
  HIGHS_CCALL(loader().Highs_writeModel(lp(), file.c_str()));
}

void HighsBackend::DoWriteSolution(const std::string &file) {
  HIGHS_CCALL(loader().Highs_writeSolutionPretty(lp(), file.c_str()));
}


void HighsBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptHighs, lp());
}

void HighsBackend::Solve() {
  int status = loader().Highs_run(lp());
  if (status != kHighsStatusOk && status != kHighsStatusWarning)
    throw std::runtime_error(fmt::format("  Error {} while solving with HiGHS").c_str());
    // Mask warnings for pdlp
  if((status == kHighsStatusWarning) && (storedOptions_.lpmethod_!="pdlp-gpu"))
      fmt::print("  Warning code {} while solving with HiGHS\n", status);
  WindupHIGHSSolve();
}

void HighsBackend::WindupHIGHSSolve() { }

void HighsBackend::ReportResults() {
   ReportHIGHSResults();
  BaseBackend::ReportResults();
}

void HighsBackend::ReportHIGHSResults() {
  SetStatus( GetSolveResult() );
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
      auto res = loader().Highs_getColIntegrality(lp(), j, &integr);
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
      loader().Highs_getColsByRange(lp(), 0, NumVars()-1, &ncols,
        costs.data(), lb.data(), ub.data(), &numnz, NULL, NULL, NULL);
      loader().Highs_changeColsBoundsByMask(lp(), s0.data(), x0.data(), x0.data());
      loader().Highs_run(lp());
      x0 = PrimalSolution();         // get new solution
      loader().Highs_changeColsBoundsByMask(lp(), s0.data(), lb.data(), ub.data());
    }
  }
  HIGHS_CCALL(loader().Highs_setSolution(lp(), x0.data(), NULL, NULL, pi0.data()));
}


ArrayRef<int> HighsBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  conStatiii_.resize(NumLinCons());
  HIGHS_CCALL(loader().Highs_getBasis(lp(), vars.data(), conStatiii_.data()));
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
    loader().Highs_getColsBySet(lp(), indicesOfMissing.size(), indicesOfMissing.data(),
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
  HIGHS_CCALL(loader().Highs_setBasis(lp(), stt.data(), cstt.data()));
}

ArrayRef<double> HighsBackend::Ray() {
  HighsInt has_ray;
  std::vector<double> uray_pres(NumVars());
  auto res = loader().Highs_getPrimalRay(lp(), &has_ray, uray_pres.data());
  if (res)
    fmt::print("Error while getting primal ray");
  if (res || (!has_ray))
  {
    uray_pres.clear();
    return uray_pres;
  }
  auto mv = GetValuePresolver().PostsolveSolution({ uray_pres });
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> HighsBackend::DRay() {
  HighsInt has_ray;
  std::vector<double> dray_pres(NumLinCons());
  auto res = loader().Highs_getDualRay(lp(), &has_ray, dray_pres.data());
  if (res)
    fmt::print("Error while getting dual ray");
  if (res || (!has_ray))
  {
    dray_pres.clear();
    return dray_pres;
  }
  auto vm = GetValuePresolver().PostsolveSolution({
                                               {},
                                               {{{CG_Linear, std::move(dray_pres)}}}
    });
  return vm.GetConValues().MoveOut();        // need the vector itself
}



void HighsBackend::AddHIGHSMessages() {
  auto pdlp = PdlpIterations();
  if (pdlp > 0)
  {
    AddToSolverMessage(fmt::format("{} PDLP iterations\n", pdlp));
    return;
  }
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

std::pair<int, std::string> HighsBackend::GetSolveResult() {
  namespace sol = mp::sol;
  int optstatus = loader().Highs_getModelStatus(lp());
  auto obj = loader().Highs_getObjectiveValue(lp());
  auto inf = loader().Highs_getInfinity(lp());
  bool hasSol = (-inf < obj && obj < inf);
  int primal_solution_status;
  loader().Highs_getIntInfoValue(lp(),
                        "primal_solution_status", &primal_solution_status);
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
    return { sol::LIMIT_INF_UNB, "unbounded or infeasible. "
                                "Disable dual reductions "
                                "or run IIS finder for definitive answer." };
  case kHighsModelStatusModelError:
  case kHighsModelStatusLoadError:
    if (kHighsSolutionStatusInfeasible == primal_solution_status)   // HiGHS 7
      return { sol::INFEASIBLE, "infeasible problem" };
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



////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo lp_values_method[] = {
  { "choose", "Automatic (default)", -1},
  { "simplex", "Simplex", 1},
  { "ipm", "Interior Point Method", 2},
  { "pdlp", "cuPDLP-c solver", 3},
  { "pdlp-gpu", "cuPDLP-c solver on NVIDIA GPU. Requires CUDA v12, not available on MacOS", 3},
};

static const mp::OptionValueInfo off_on_choose_values[] = {
  { "choose", "Automatic (default)", -1},
  { "off", "Off", 1},
  { "on", "On", 2},
};


static const mp::OptionValueInfo pdlperestartmethod_values[] = {
  { "0", "None", 0},
  { "1", "GPU (default)", 1},
  { "2", "CPU", 2},
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
    "  ampl: option highs_options 'mip:gap=1e-6';\n");

  AddSolverOption("tech:outlev outlev",
    "0*/1: Whether to write HighS log lines (chatter) to stdout and to file.",
    "output_flag", 0, 1);

  AddSolverOption("tech:logfile logfile",
    "Log file name.", "log_file");

  std::string c;
  AddStoredOption("alg:method method lpmethod solver",
    "Which algorithm to use :\n"
    "\n.. value-table::\n", storedOptions_.lpmethod_, lp_values_method);

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

  AddSolverOption("pre:userboundscale user_bound_scale userboundscale",
    "Exponent of power-of-two bound scaling for model (default 0).",
    "user_bound_scale", 0, std::numeric_limits<HighsInt>::max());

  AddSolverOption("pre:usercostscale user_cost_scale usercostscale",
    "Exponent of power-of-two cost scaling for model (default 0).",
    "user_cost_scale", 0, std::numeric_limits<HighsInt>::max());

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

  AddSolverOption("lim:objectivebound objective_bound objectivebound",
    "Objective bound for termination of the dual simplex solver (default: no limit).",
    "objective_bound", 0.0, DBL_MAX);


  AddSolverOption("lim:objectivetarget objective_target objectivetarget",
    "Objective target for termination of the MIP solver (default: no limit).",
    "objective_target", 0.0, DBL_MAX);

  AddSolverOption("lim:pdlpnativetermination pdlp_native_termination pdlpnativetermination",
    "Use native termination for PDLP solver:\n"
    "\n.. value-table::\n", "pdlp_native_termination", values_01_noyes_0default_, 0);

  AddSolverOption("pre:pdlpscaling pdlp_scaling pdlpscaling",
    "Scaling option for PDLP solver:\n"
    "\n.. value-table::\n", "pdlp_scaling", values_01_noyes_1default_, 1);

  AddSolverOption("lim:pdlpiterationlimit pdlpiterationlimit pdlp_iteration_limit",
    "Iteration limit for PDLP solver (default: no limit).",
    "pdlp_iteration_limit", 0, INT_MAX);

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

  AddSolverOption("alg:dualrestol dualrestol dual_residual_tolerance",
    "Dual residual tolerance (default 1e-7).",
    "dual_residual_tolerance", 1e-10, Infinity());

  AddSolverOption("alg:ipmopttol ipmopttol ipm_optimality_tolerance",
    "IPM optimality tolerance (default 1e-8).",
    "ipm_optimality_tolerance", 1e-12, Infinity());

  AddSolverOption("alg:pdlperestartmethod pdlperestartmethod pdlp_e_restart_method",
    "Restart mode for PDLP solver (default 1).",
    "pdlp_e_restart_method", pdlperestartmethod_values, 1);

  AddSolverOption("alg:pdlpdgaptol pdlpdgaptol pdlp_d_gap_tol",
    "Duality gap tolerance for PDLP solver (default 1e-4).",
    "pdlp_d_gap_tol", 1e-12, Infinity());

  AddSolverOption("bar:crossover crossover run_crossover",
    "Run crossover after IPM to get a basic solution",
    "run_crossover", run_crossover_values, c);

  AddSolverOption("tech:threads threads",
    "How many threads to use when using the barrier algorithm "
    "or solving MIP problems; default 0 ==> automatic choice.",
    "threads", 0, INT32_MAX);
  AddSolverOption("mip:detsimmetry detsimmetry mip_detect_symmetry",
    "Whether symmetry should be detected (default 1)",
    "mip_detect_symmetry", 0, 1);

  AddSolverOption("mip:lifting lifting mip_lifting_for_probing",
    "Whether lifting for probing should be used (default -1)",
    "mip_lifting_for_probing", -1, INT_MAX);

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

  AddSolverOption("tech:seed seed random_seed",
    "Random number seed (default 0), affecting perturbations that "
    "may influence the solution path.",
    "random_seed", 0, INT_MAX);

  AddSolverOption("mip:intfeastol intfeastol inttol mip_feasibility_tolerance",
    "Feasibility tolerance for integer variables (default 1e-06).",
    "mip_feasibility_tolerance", 1e-10, Infinity());

  AddSolverOption("mip:heureff heureff mip_heuristic_effort",
    "Fraction of time to spend in MIP heuristics (default 0.05).",
    "mip_heuristic_effort", 0.0, 1.0);

  AddSolverOption("mip:gap mipgap mip:relgaptol relgaptol mip_rel_gap",
    "Tolerance on relative gap, | ub - lb|/|ub | , to determine whether optimality has been reached for a MIP instance "
    "(default 1e-04).",
    "mip_rel_gap", 0.0, Infinity());

  AddSolverOption("mip:gapabs mipgapabs mip:absgaptol absgaptol mip_abs_gap",
    "Tolerance on absolute gap of MIP, |ub-lb|, to determine whether optimality has been reached for a MIP instance "
    "(default 1e-06).",
    "mip_abs_gap", 0.0, Infinity());

  AddSolverOption("pre:centring run_centring centring",
    "Perform centring steps or not:\n"
    "\n.. value-table::\n", "run_centring", values_01_noyes_0default_, 0);

  AddSolverOption("pre:maxcentringsteps max_centring_steps maxcentringsteps",
    "Maximum number of steps to use when computing the analytic centre "
    "(default 5).",
    "max_centring_steps", 0, INT_MAX);

  AddSolverOption("pre:centringratiotolerance centring_ratio_tolerance centringratiotolerance",
    "Centring stops when the ratio max(x_j*s_j) / min(x_j*s_j) is below "
    "this tolerance (default 100).",
    "centring_ratio_tolerance", 0, INT_MAX);

  AddSolverOption("obj:blend blend_multi_objectives",
    "Whether to blend multiple objectives or apply lexicographical ordering",
    "blend_multi_objectives", values_01_noyes_1default_, 1
  );
}

double HighsBackend::MIPGap() {
  if (BarrierIterations() == 0)
    return 0;
  return getDblAttr("mip_gap");
}
double HighsBackend::BestDualBound() {
  if (BarrierIterations() == 0)
    return 0;
  return getDblAttr("mip_dual_bound");
}

double HighsBackend::MIPGapAbs() {
  if (BarrierIterations() == 0)
    return 0;
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}
} // namespace mp


// AMPLs

AMPLS_MP_Solver* Open_highs(CCallbacks cb = {}) {
  AMPLS_MP_Solver* slv =
    AMPLS__internal__Open(
      std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()},
      cb);
  return slv;
}

void AMPLSClose_highs(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_highs(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::HighsBackend*>(AMPLSGetBackend(slv))->lp();
}


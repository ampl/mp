#include <atomic>
#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mosekbackend.h"

extern "C" {
  #include "mosek-ampls-c-api.h"    // Mosek AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {

/// Set by the AMPL interrupt handler, polled by MOSEK's callback
std::atomic<bool> interrupt_requested{false};

bool InterruptMosek(void*) {
  interrupt_requested = true;
  return true;
}

/// MOSEK progress callback: a nonzero return stops the optimizer
int MSKAPI MosekProgressCallback(MSKtask_t, MSKuserhandle_t, MSKcallbackcodee,
                                 const MSKrealt*, const MSKint32t*, const MSKint64t*) {
  return interrupt_requested ? 1 : 0;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateMosekBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::MosekBackend()};
}


namespace mp {

/// Create Mosek Model Manager
/// @param gc: the Mosek common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return MosekModelMgr
std::unique_ptr<BasicModelManager>
CreateMosekModelMgr(MosekCommon&, Env&, pre::BasicValuePresolver*&);


MosekBackend::MosekBackend() {
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateMosekModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

MosekBackend::~MosekBackend() {
  CloseSolver();
}

static void MSKAPI printstr(void* handle, const char* str)
{
  fmt::print("{}", str);
  fflush(stdout);
}

void MosekBackend::OpenSolver() {
  int status = MSK_RES_OK;
  const auto& initialize = GetCallbacks().init;

  MSKenv_t env = NULL;
  MSKtask_t task;
  if (initialize) { // If an initialization callback is provided,
    // use it to create the environment
    MOSEK_CCALL(MSK_makeenv(&env, NULL));
    env = (MSKenv_t)initialize();
    set_env(env);
  }
  else {
    MSKrescodee res = MSK_makeenv(&env, NULL);
    set_env(env);
    // Check license here
    if (MSK_checkoutlicense(env, MSK_FEATURE_PTS))
    {
      const auto diag = GetCallbacks().diagnostics;
      if (diag) {
        // If a diagnostic function is provided, do not print more
        // information
        diag();
        exit(res);
      }
    }
  }
  status = MSK_maketask(env, 0, 0, &task);
  if (status)
    throw std::runtime_error(fmt::format(
      "Failed to create task, error code {}.", status));
  set_lp(task); // Assign it

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();

  /// Turn off verbosity by default
  MOSEK_CCALL(MSK_putintparam(task, MSK_IPAR_LOG, 0));
  
  // Register callback for console logging (controlled by the outlev param
  // in all AMPL solver drivers)
  MSK_linkfunctotaskstream(lp(), MSK_STREAM_LOG, NULL, printstr);
  // Progress callback, for interrupts (Ctrl-C)
  MOSEK_CCALL(MSK_putcallbackfunc(lp(), MosekProgressCallback, nullptr));

}

void MosekBackend::CloseSolver() {
  if (fixed_task_) MSK_deletetask(&fixed_task_);
  if (lp()) MSK_deletetask(&lp_ref());
  if (env()) MSK_deleteenv(&env_ref());
}

const char* MosekBackend::GetBackendName()
  { return "MosekBackend"; }

std::string MosekBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", MSK_VERSION_MAJOR,
    MSK_VERSION_MINOR, MSK_VERSION_REVISION);
}

void MosekBackend::InitOptionParsing() {
  OpenSolver();
}

bool MosekBackend::IsMIP() const {
  int num_int = 0;
  MOSEK_CCALL(MSK_getnumintvar(lp(), &num_int));
  return num_int != 0;
}

bool MosekBackend::IsQCP() const {
  MSKproblemtypee type;
  MOSEK_CCALL(MSK_getprobtype(lp(), &type));
  return MSK_PROBTYPE_QCQO == type;
}

ArrayRef<double> MosekBackend::PrimalSolution() {
  if (HasFixedLP()) return FixedLP().primal;
	std::vector<double> x;
  // MSKsolstae solst;
  // MSK_getsolsta(lp(), solToFetch_, &solst);
  if (true    // Mosek 11.0.8 seems to report UNKNOWN
      // even with feasible solutions when stalling
   //    || (MSK_SOL_STA_PRIM_INFEAS_CER!=solst
      // && MSK_SOL_STA_UNKNOWN!=solst)
      ) {
		int num_vars = NumVars();
		x.resize(num_vars);
		// TODO get appropriate solution?
		auto status = MSK_getxx(lp(), solToFetch_, x.data());
		if (MSK_RES_OK != status)
			x.clear();
	}
	return x;
}

pre::ValueMapDbl MosekBackend::DualSolution() {
	return {{
			{ CG_Algebraic, DualSolution_LP() }
		}};
}

ArrayRef<double> MosekBackend::DualSolution_LP() {
  if (IsMIP()) return HasFixedLP() ? ArrayRef<double>(FixedLP().dual) : ArrayRef<double>{};
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  // if (!IsMIP())
  {
    MSKrescodee error = MSK_gety(lp(), solToFetch_, pi.data());
    if (error != MSK_RES_OK)
      pi.clear();
  }
  return pi;
}

ArrayRef<double> MosekBackend::DualSolution_Cones() {
	// TODO can we return anything sensible?
	// Variable suffixes?
	return {};
}

double MosekBackend::ObjectiveValue() const {
  double v;
  // TODO get appropriate solution
  MOSEK_CCALL(MSK_getprimalobj(lp(), solToFetch_, &v));
  return v;
}

int MosekBackend::NodeCount() const {
  return getIntAttr(MSK_IINF_MIO_NUM_BRANCH);
}

int MosekBackend::SimplexIterations() const {
  return getIntAttr(MSK_IINF_SIM_PRIMAL_ITER) + getIntAttr(MSK_IINF_SIM_DUAL_ITER);
}

int MosekBackend::BarrierIterations() const {
  return getIntAttr(MSK_IINF_INTPNT_ITER);
}

std::map<std::string, std::variant<int, double, std::string>>
MosekBackend::SolutionStats() {
    std::map<std::string, std::variant<int, double, std::string>> stats;
    stats["simplex_iterations"] = SimplexIterations();
    stats["barrier_iterations"] = BarrierIterations();
    stats["node_count"] = NodeCount();
    return stats;
}

void MosekBackend::DoWriteProblem(const std::string &file) {
  MOSEK_CCALL(MSK_writedata(lp(), file.data()));
}

void MosekBackend::DoWriteSolution(const std::string& name) {
  MOSEK_CCALL(MSK_writesolution(lp(), solToFetch_, name.c_str()));
}


void MosekBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptMosek, lp());
}

/// The solution whose values and status keys start the optimizer
MSKsoltypee MosekBackend::GetSolutionTypeForStart() {
  if (IsMIP())
    return MSK_SOL_ITG;
  MSKproblemtypee type;
  MOSEK_CCALL(MSK_getprobtype(lp(), &type));
  return MSK_PROBTYPE_LO == type ? MSK_SOL_BAS : MSK_SOL_ITR;
}

MSKsoltypee MosekBackend::GetSolutionTypeToFetch() {
  int b = 0;

  // Handle MIP
  MOSEK_CCALL(MSK_solutiondef(lp(), MSK_SOL_ITG, &b));
  if (b) return MSK_SOL_ITG;

  // Handle LP
  MSKproblemtypee type;
  MOSEK_CCALL(MSK_getprobtype(lp(), &type));
  if (type == MSK_PROBTYPE_LO)
  {
    MOSEK_CCALL(MSK_solutiondef(lp(), MSK_SOL_BAS, &b));
    if (b) return MSK_SOL_BAS;
  }

  // All other case
  return MSK_SOL_ITR;
}

void MosekBackend::Solve() {
  auto res = MSK_optimizetrm(lp(), &termCode_);
  if (MSK_RES_ERR_CON_Q_NOT_PSD == res) {
    AddWarning("Nonconvex QC",
               "Mosek reported the problem as nonconvex QCP.\n"
               "If the constraints are in fact Second-Order Cones,\n"
               "make sure they have standard SOCP forms, and\n"
               "the objective is linear by moving the\n"
               "quadratic terms into auxiliary SOC constraints.\n"
               "See mp.ampl.com/modeling-expressions.html#conic-optimization.");

  }
  MOSEK_CCALL( res );
  solToFetch_ = GetSolutionTypeToFetch();
  MOSEK_CCALL(MSK_getprosta(lp(), solToFetch_, &proSta_));
  MOSEK_CCALL(MSK_getsolsta(lp(), solToFetch_, &solSta_));
  WindupMOSEKSolve();
}

MosekBackend::FixedLPResult MosekBackend::SolveFixedLP(ArrayRef<double> x) {
  if (solToFetch_ != MSK_SOL_ITG ||
      (solSta_ != MSK_SOL_STA_INTEGER_OPTIMAL && solSta_ != MSK_SOL_STA_PRIM_FEAS))
    MP_RAISE("no integer-feasible incumbent");
  if (fixed_task_) MSK_deletetask(&fixed_task_);
  MOSEK_CCALL(MSK_clonetask(lp(), &fixed_task_));
  int n = NumVars(), m = NumLinCons();
  double tol;
  GetSolverOption(MSK_DPAR_MIO_TOL_ABS_RELAX_INT, tol);
  for (int j = 0; j < n; ++j) {
    MSKvariabletypee type;
    MOSEK_CCALL(MSK_getvartype(lp(), j, &type));
    if (type == MSK_VAR_TYPE_INT) {
      MSKboundkeye key;
      double lb, ub;
      MOSEK_CCALL(MSK_getvarbound(lp(), j, &key, &lb, &ub));
      if (key == MSK_BK_UP || key == MSK_BK_FR) lb = -INFINITY;
      if (key == MSK_BK_LO || key == MSK_BK_FR) ub = INFINITY;
      double value = FixedIntegerValue(x[j], lb, ub, tol);
      MOSEK_CCALL(MSK_putvartype(fixed_task_, j, MSK_VAR_TYPE_CONT));
      MOSEK_CCALL(MSK_putvarbound(fixed_task_, j, MSK_BK_FX, value, value));
    }
  }
  MOSEK_CCALL(MSK_putintparam(fixed_task_, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_DUAL_SIMPLEX));
  MOSEK_CCALL(MSK_putdouparam(fixed_task_, MSK_DPAR_OPTIMIZER_MAX_TIME,
      fixed_model_time_limit() < 0 ? -1.0 : fixed_model_time_limit()));
  MOSEK_CCALL(MSK_putintparam(fixed_task_, MSK_IPAR_SIM_MAX_ITERATIONS,
      fixed_model_iter_limit() < 0 ? INT_MAX : fixed_model_iter_limit()));
  MSKrescodee termination;
  MOSEK_CCALL(MSK_optimizetrm(fixed_task_, &termination));
  MSKsolstae status;
  MOSEK_CCALL(MSK_getsolsta(fixed_task_, MSK_SOL_BAS, &status));
  if (status != MSK_SOL_STA_OPTIMAL)
    MP_RAISE(fmt::format("fixed LP not optimal (MOSEK status {}, termination {})", (int)status, (int)termination));
  FixedLPResult result;
  MOSEK_CCALL(MSK_getprimalobj(fixed_task_, MSK_SOL_BAS, &result.objective));
  result.primal.resize(n);
  result.dual.resize(m);
  std::vector<MSKstakeye> vst(n), cst(m);
  MOSEK_CCALL(MSK_getxx(fixed_task_, MSK_SOL_BAS, result.primal.data()));
  MOSEK_CCALL(MSK_gety(fixed_task_, MSK_SOL_BAS, result.dual.data()));
  MOSEK_CCALL(MSK_getskx(fixed_task_, MSK_SOL_BAS, vst.data()));
  MOSEK_CCALL(MSK_getskc(fixed_task_, MSK_SOL_BAS, cst.data()));
  auto convert = [](MSKstakeye status) {
    switch (status) {
    case MSK_SK_BAS: return (int)BasicStatus::bas;
    case MSK_SK_LOW: return (int)BasicStatus::low;
    case MSK_SK_UPR: return (int)BasicStatus::upp;
    case MSK_SK_FIX: return (int)BasicStatus::equ;
    case MSK_SK_SUPBAS: return (int)BasicStatus::sup;
    default: MP_RAISE(fmt::format("Unknown MOSEK fixed-LP basis status {}", (int)status));
    }
  };
  for (auto status : vst) result.varstt.push_back(convert(status));
  for (auto status : cst) result.constt.push_back(convert(status));
  int primal = 0, dual = 0;
  MOSEK_CCALL(MSK_getintinf(fixed_task_, MSK_IINF_SIM_PRIMAL_ITER, &primal));
  MOSEK_CCALL(MSK_getintinf(fixed_task_, MSK_IINF_SIM_DUAL_ITER, &dual));
  result.iterations = primal + dual;
  return result;
}

void MosekBackend::WindupMOSEKSolve() { }

void MosekBackend::ReportResults() {
  ReportMOSEKResults();
  ConsiderEmulatedFixedModel();
  BaseBackend::ReportResults();
}

void MosekBackend::ReportMOSEKResults() {
  SetStatus( GetSolveResult() );
  AddMOSEKMessages();
  if (need_multiple_solutions())
    ReportMOSEKPool();
}

std::vector<double> MosekBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
  // TODO get solutions in the pool 
  return vars;
}

double MosekBackend::getPoolObjective(int i)
{
	double obj=0.0;
  // TODO get objective value of solution i
  return obj;
}

void MosekBackend::ReportMOSEKPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  // TODO
  /*
  while (++iPoolSolution < getIntAttr(MOSEK_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void MosekBackend::AddMOSEKMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

/// Result for a termination code (limit, interrupt, numerical trouble)
std::pair<int, std::string> MosekBackend::GetLimitResult(bool has_solution) {
  namespace sol = mp::sol;
  auto result = [has_solution](int feas, int nofeas, const std::string& text) {
    return std::pair<int, std::string>{
      has_solution ? feas : nofeas,
      text + (has_solution ? ", feasible solution" : ", no solution") };
  };
  switch (termCode_) {
  case MSK_RES_TRM_MAX_ITERATIONS:
    return result(sol::LIMIT_FEAS_ITER, sol::LIMIT_NO_FEAS_ITER, "iteration limit");
  case MSK_RES_TRM_MAX_TIME:
  case MSK_RES_TRM_SERVER_MAX_TIME:
    return result(sol::LIMIT_FEAS_TIME, sol::LIMIT_NO_FEAS_TIME, "time limit");
  case MSK_RES_TRM_OBJECTIVE_RANGE:
    return result(sol::LIMIT_FEAS_BESTOBJ, sol::LIMIT_NO_FEAS_CUTOFF,
                  "objective range (lim:lowerobjcut, lim:upperobjcut)");
  case MSK_RES_TRM_MIO_NUM_BRANCHES:
    return result(sol::LIMIT_FEAS_NODES, sol::LIMIT_NO_FEAS_NODES, "branch limit");
  case MSK_RES_TRM_MIO_NUM_RELAXS:
    return result(sol::LIMIT_FEAS_NODES, sol::LIMIT_NO_FEAS_NODES, "relaxation limit");
  case MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS:
    return result(sol::LIMIT_FEAS_NUMSOLS, sol::LIMIT_NO_FEAS, "solution limit");
  case MSK_RES_TRM_USER_CALLBACK:
    return result(sol::LIMIT_FEAS_INTERRUPT, sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted");
  case MSK_RES_TRM_SERVER_MAX_MEMORY:
    return result(sol::LIMIT_FEAS_SOFTMEM, sol::LIMIT_NO_FEAS_SOFTMEM, "memory limit");
  case MSK_RES_TRM_STALL:
    return result(sol::UNCERTAIN, sol::NUMERIC, "stalling");
  case MSK_RES_TRM_NUMERICAL_PROBLEM:
  case MSK_RES_TRM_MAX_NUM_SETBACKS:
    return result(sol::UNCERTAIN, sol::NUMERIC, "numerical problems");
  default:
    return result(sol::LIMIT_FEAS, sol::LIMIT_NO_FEAS,
                  "termination code " + std::to_string(termCode_));
  }
}

std::pair<int, std::string> MosekBackend::GetSolveResult() {
  namespace sol = mp::sol;
  std::string term_info = ConvertMOSEKTermStatus();
  bool has_solution = MSK_SOL_STA_PRIM_FEAS == solSta_ ||
      MSK_SOL_STA_PRIM_AND_DUAL_FEAS == solSta_ ||
      MSK_SOL_STA_OPTIMAL == solSta_ || MSK_SOL_STA_INTEGER_OPTIMAL == solSta_;
  switch (solSta_) {
  case MSK_SOL_STA_OPTIMAL:
  case MSK_SOL_STA_INTEGER_OPTIMAL:
    return { sol::SOLVED, "optimal" + term_info };
  case MSK_SOL_STA_PRIM_FEAS:
  case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
    if (MSK_RES_OK != termCode_)
      return GetLimitResult(true);
    return { sol::UNCERTAIN, "feasible, not proven optimal" + term_info };
  case MSK_SOL_STA_DUAL_FEAS:
    if (MSK_RES_OK != termCode_)
      return GetLimitResult(false);
    return { sol::UNCERTAIN, "feasible dual" + term_info };
  case MSK_SOL_STA_PRIM_INFEAS_CER:
    return { sol::INFEASIBLE,
          "the solution is a primal infeasibility certificate" + term_info };
  case MSK_SOL_STA_DUAL_INFEAS_CER:
    return { sol::LIMIT_INF_UNB,
          "the solution is a dual infeasibility certificate" + term_info };
  case MSK_SOL_STA_PRIM_ILLPOSED_CER:
    return { sol::NUMERIC,
          "the solution is a primal ill-posedness certificate" + term_info };
  case MSK_SOL_STA_DUAL_ILLPOSED_CER:
    return { sol::NUMERIC,
          "the solution is a dual ill-posedness certificate" + term_info };
  case MSK_SOL_STA_UNKNOWN:
  default:
    switch (proSta_) {
    case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
      return
      { sol::UNCERTAIN, "primal and dual feasible" + term_info };
    case MSK_PRO_STA_PRIM_FEAS:
      return { sol::LIMIT_FEAS, "feasible primal" + term_info };
    case MSK_PRO_STA_DUAL_FEAS:
      return { sol::UNCERTAIN, "feasible dual" + term_info };
    case MSK_PRO_STA_PRIM_INFEAS:
      return { sol::INFEASIBLE, "primal infeasible" + term_info };
    case MSK_PRO_STA_DUAL_INFEAS:
      return { sol::LIMIT_INF_UNB, "dual infeasible" + term_info };
    case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
      return
      { sol::INFEASIBLE, "primal and dual infeasible" + term_info };
    case MSK_PRO_STA_ILL_POSED:
      return { sol::NUMERIC, "ill-posed problem" + term_info };
    case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      return { sol::INF_OR_UNB, "infeasible or unbounded" + term_info };
    case MSK_PRO_STA_UNKNOWN:
    default:
      if (MSK_RES_OK != termCode_)
        return GetLimitResult(has_solution);
      return { sol::UNKNOWN,
            "unknown status (solution status: " + std::to_string(solSta_)
            + ", problem status: " + std::to_string(proSta_)
            + ")" + term_info };
    }
  }
}

std::string MosekBackend::ConvertMOSEKTermStatus() {
  switch (termCode_) {
  case MSK_RES_OK:
      return "";
  case MSK_RES_TRM_MAX_ITERATIONS:
    return { ", max number of iterations reached" };
  case MSK_RES_TRM_MAX_TIME:
    return { ", maximum allowed time reached" };
  case MSK_RES_TRM_OBJECTIVE_RANGE:
    return { ", objective range" };
  case MSK_RES_TRM_STALL:
    return { ", stalling" };
  case MSK_RES_TRM_NUMERICAL_PROBLEM:
    return { ", numerical issues" };
  case MSK_RES_TRM_MAX_NUM_SETBACKS:
    return { ", numerical issues (too many setbacks)" };
  case MSK_RES_TRM_MIO_NUM_RELAXS:
    return { ", relaxation limit hit" };
  case MSK_RES_TRM_MIO_NUM_BRANCHES:
    return { ", branch limit hit" };
  case MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS:
    return { ", solution limit hit" };
  case MSK_RES_TRM_USER_CALLBACK:
    return { ", interrupted" };
  case MSK_RES_TRM_LOST_RACE:
    return { ", lost race" };
  case MSK_RES_TRM_INTERNAL:
  case MSK_RES_TRM_INTERNAL_STOP:
    return { ", internal stop" };
  case MSK_RES_TRM_SERVER_MAX_TIME:
    return { ", server time limit" };
  case MSK_RES_TRM_SERVER_MAX_MEMORY:
    return { ", server memory limit" };
  default:
    return { ", termination code "
          + std::to_string(termCode_) };
  }
  return { ", termination reason unknown" };
}

void MosekBackend::FinishOptionParsing() {
  int log_level = 0;
  GetSolverOption(MSK_IPAR_LOG, log_level);
  set_verbose_mode(log_level > 0);
  if (!logFile().empty()) {
    MOSEK_CCALL(MSK_linkfiletotaskstream(lp(), MSK_STREAM_LOG, logFile().c_str(), false));
    if (0 == log_level) {         // log file only: standard log, quiet console
      SetSolverOption(MSK_IPAR_LOG, 1);
      MOSEK_CCALL(MSK_unlinkfuncfromtaskstream(lp(), MSK_STREAM_LOG));
    }
  }
  // Native params
  if (paramfile_read().size())
    MOSEK_CCALL(
      MSK_readparamfile(lp(), paramfile_read().c_str()));
  /// Set advanced parameters
  for (const auto& prm : storedOptions_.inlineParams_) {
    auto sp1 = prm.find_first_of(' ');
    auto sp2 = prm.find_last_of(' ');
    if (sp1 != sp2 || std::string::npos == sp1
      || 0==sp1 || prm.back()==' ')
      MP_RAISE(
        fmt::format("Bad native option assignment:\n"
          "   '{}',\n"
          "   need 'name value'", prm));
    MOSEK_CCALL( MSK_putparam(lp(),
      prm.substr(0, sp1).c_str(), prm.substr(sp1+1, prm.size()-sp1-1).c_str()));
  }
  // Write native params
  if (paramfile_write().size())
    MOSEK_CCALL(
      MSK_writeparamfile(lp(), paramfile_write().c_str()));
  // Scaling works for both
  if (storedOptions_.PREScale_ != 0) {
    SetSolverOption(MSK_IPAR_SIM_SCALING, storedOptions_.PREScale_);
    SetSolverOption(MSK_IPAR_INTPNT_SCALING, storedOptions_.PREScale_);
  }
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo alg_values_method[] = {
  { "0", "Optimizer for conic constraints", MSK_OPTIMIZER_CONIC},
  { "1", "Dual simplex", MSK_OPTIMIZER_DUAL_SIMPLEX},
  { "2", "Automatic (default)", MSK_OPTIMIZER_FREE},
  { "3", "Free simplex", MSK_OPTIMIZER_FREE_SIMPLEX},
  { "4", "Interior-point method", MSK_OPTIMIZER_INTPNT},
  { "5", "Mixed-integer optimizer", MSK_OPTIMIZER_MIXED_INT},
  { "6", "New dual simplex (MOSEK 11)", MSK_OPTIMIZER_NEW_DUAL_SIMPLEX},
  { "7", "New primal simplex (MOSEK 11)", MSK_OPTIMIZER_NEW_PRIMAL_SIMPLEX},
  { "8", "Primal simplex", MSK_OPTIMIZER_PRIMAL_SIMPLEX}
};

static const mp::OptionValueInfo values_bar_intpnt_basis[] = {
    { "0", "Never", MSK_BI_NEVER},
    { "1", "Always (default)", MSK_BI_ALWAYS},
    { "2", "If no error detected", MSK_BI_NO_ERROR},
    { "3", "If primal and dual feasible", MSK_BI_IF_FEASIBLE}
};

static const mp::OptionValueInfo alg_values_presolve_use[] = {
  { "0", "Do not use presolve", 0},
  { "1", "Use presolve", 1},
  { "2", "Automatic (default)", 2}
};
static const mp::OptionValueInfo values_prefoldinguse_[] = {
  { "0", "Disabled", MSK_FOLDING_MODE_OFF},
  { "1", "The solver decides on the usage and amount of folding", MSK_FOLDING_MODE_FREE },
  { "2", "If only the interior-point solution is requested then the solver decides; if the basic solution is requested then folding is disabled (default)", MSK_FOLDING_MODE_FREE_UNLESS_BASIC},
  { "3", "Full folding is always performed regardless of workload", MSK_FOLDING_MODE_FORCE },
};
static const mp::OptionValueInfo values_premipfoldinguse_[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Disabled", 0},
  { "1", "Low amount", 1},
  { "2", "Medium amount", 2},
  { "3", "High amount", 3},
  { "4", "Extremely high amount", 4},
};

static const mp::OptionValueInfo values_mip_feaspump_[] = {
    { "-1", "Automatic (default)", -1},
    { "0", "The Feasibility Pump is disabled", 0},
    { "1", "The Feasibility Pump is enabled with an effort to improve solution quality", 1},
    { "2", "The Feasibility Pump is enabled with an effort to reach feasibility early", 2}
};

static const mp::OptionValueInfo values_mip_var_selection_types[] = {
  { "0", "Automatic (default)", 0},
  { "1", "Pseudocost variable selection", 1},
  { "2", "Strong branching selection", 2}
};

static const mp::OptionValueInfo values_mip_presolve_dual_ray[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Disabled", 0},
  { "1", "Low amount of analyis", 1},
  { "2", "Higher amount of analysis", 2}
};

static const mp::OptionValueInfo values_mio_optimizer[] = {
  { "0", "Optimizer for conic constraints", MSK_OPTIMIZER_CONIC},
  { "1", "Dual simplex", MSK_OPTIMIZER_DUAL_SIMPLEX},
  { "2", "Automatic (default)", MSK_OPTIMIZER_FREE},
  { "3", "Free simplex", MSK_OPTIMIZER_FREE_SIMPLEX},
  { "4", "Interior-point method", MSK_OPTIMIZER_INTPNT},
  { "6", "New dual simplex", MSK_OPTIMIZER_NEW_DUAL_SIMPLEX},
  { "7", "New primal simplex", MSK_OPTIMIZER_NEW_PRIMAL_SIMPLEX},
  { "8", "Primal simplex", MSK_OPTIMIZER_PRIMAL_SIMPLEX}
};

static const mp::OptionValueInfo values_mio_nodeselect[] = {
  { "0", "Automatic (default)", MSK_MIO_NODE_SELECTION_FREE},
  { "1", "Depth first (first)", MSK_MIO_NODE_SELECTION_FIRST},
  { "2", "Best bound", MSK_MIO_NODE_SELECTION_BEST},
  { "3", "Best estimate using pseudo costs", MSK_MIO_NODE_SELECTION_PSEUDO}
};

static const mp::OptionValueInfo values_mio_branchdir[] = {
  { "0", "Automatic (default)", MSK_BRANCH_DIR_FREE},
  { "1", "Up", MSK_BRANCH_DIR_UP},
  { "2", "Down", MSK_BRANCH_DIR_DOWN},
  { "3", "Nearest integer", MSK_BRANCH_DIR_NEAR},
  { "4", "Farthest integer", MSK_BRANCH_DIR_FAR},
  { "5", "Towards the root LP value", MSK_BRANCH_DIR_ROOT_LP},
  { "6", "Guided by the incumbent", MSK_BRANCH_DIR_GUIDED},
  { "7", "By pseudo costs", MSK_BRANCH_DIR_PSEUDOCOST}
};

static const mp::OptionValueInfo values_mio_cutselection[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Generated cuts are less likely to be added", 0},
  { "1", "Cuts are selected more aggressively", 1}
};

static const mp::OptionValueInfo values_mio_probing[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Disabled", 0},
  { "1", "Low amount of probing", 1},
  { "2", "Medium amount of probing", 2},
  { "3", "High amount of probing", 3}
};

static const mp::OptionValueInfo values_sim_selection[] = {
  { "0", "Automatic (default)", MSK_SIM_SELECTION_FREE},
  { "1", "Full pricing", MSK_SIM_SELECTION_FULL},
  { "2", "Approximate steepest edge", MSK_SIM_SELECTION_ASE},
  { "3", "Devex pricing", MSK_SIM_SELECTION_DEVEX},
  { "4", "Steepest edge", MSK_SIM_SELECTION_SE},
  { "5", "Partial pricing", MSK_SIM_SELECTION_PARTIAL}
};

static const mp::OptionValueInfo values_sim_hotstart[] = {
  { "0", "No hot start", MSK_SIM_HOTSTART_NONE},
  { "1", "Automatic (default)", MSK_SIM_HOTSTART_FREE},
  { "2", "Hot start from the status keys only", MSK_SIM_HOTSTART_STATUS_KEYS}
};

static const mp::OptionValueInfo values_sim_degen[] = {
  { "0", "No degeneration handling", MSK_SIM_DEGEN_NONE},
  { "1", "Automatic (default)", MSK_SIM_DEGEN_FREE},
  { "2", "Aggressive", MSK_SIM_DEGEN_AGGRESSIVE},
  { "3", "Moderate", MSK_SIM_DEGEN_MODERATE},
  { "4", "Minimum", MSK_SIM_DEGEN_MINIMUM}
};

static const mp::OptionValueInfo values_solve_form[] = {
  { "0", "Automatic (default)", MSK_SOLVE_FREE},
  { "1", "Solve the primal problem", MSK_SOLVE_PRIMAL},
  { "2", "Solve the dual problem", MSK_SOLVE_DUAL}
};

static const mp::OptionValueInfo values_bar_order[] = {
  { "0", "Automatic (default)", MSK_ORDER_METHOD_FREE},
  { "1", "Approximate minimum local fill-in", MSK_ORDER_METHOD_APPMINLOC},
  { "2", "Experimental", MSK_ORDER_METHOD_EXPERIMENTAL},
  { "3", "Try graph partitioning", MSK_ORDER_METHOD_TRY_GRAPHPAR},
  { "4", "Force graph partitioning", MSK_ORDER_METHOD_FORCE_GRAPHPAR},
  { "5", "No ordering", MSK_ORDER_METHOD_NONE}
};

/// MIP cut families, on/off
struct MosekCutFamily {
  const char* names;
  MSKiparame key;
  const char* descr;
};
static const MosekCutFamily mosek_cut_families[] = {
  { "cut:clique cliquecuts", MSK_IPAR_MIO_CUT_CLIQUE, "clique cuts (default 1)" },
  { "cut:mir mircuts cut:cmir", MSK_IPAR_MIO_CUT_CMIR, "complemented mixed-integer rounding cuts (default 1)" },
  { "cut:gomory gomorycuts cut:gmi", MSK_IPAR_MIO_CUT_GMI, "Gomory mixed-integer cuts (default 1)" },
  { "cut:implied impliedcuts cut:impliedbound", MSK_IPAR_MIO_CUT_IMPLIED_BOUND, "implied bound cuts (default 1)" },
  { "cut:cover covercuts cut:knapsackcover", MSK_IPAR_MIO_CUT_KNAPSACK_COVER, "knapsack cover cuts (default 1)" },
  { "cut:lipro liprocuts", MSK_IPAR_MIO_CUT_LIPRO, "lift-and-project cuts (default 0)" }
};

static const mp::OptionValueInfo values_bar_presolve_scaling[] = {
  {"0", "Mosek chooses scaling heuristic (default)", MSK_SCALING_FREE},
  {"1", "No scaling.", MSK_SCALING_NONE},
};


void MosekBackend::InitCustomOptions() {

  set_option_header(
      "MOSEK Optimizer Options for AMPL\n"
      "--------------------------------\n\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``mosek_options``. For example::\n"
      "\n"
      "  ampl: option mosek_options 'threads=3';\n");

  // ********** ALG OPTIONS **********
  AddMosekOption("alg:method method lpmethod simplex",
                  "Which algorithm to use for non-MIP problems or for the root node of MIP problems:\n"
                  "\n.. value-table::\n", MSK_IPAR_OPTIMIZER, alg_values_method, 2);

  AddMosekOption("alg:feastol basistolx",
                 "Primal feasibility tolerance of the simplex optimizers "
                 "(MSK_DPAR_BASIS_TOL_X, min 1e-9, default 1e-6).",
                 MSK_DPAR_BASIS_TOL_X, 1e-9, Infinity());

  AddMosekOption("alg:dualfeastol basistols",
                 "Dual feasibility tolerance of the simplex optimizers "
                 "(MSK_DPAR_BASIS_TOL_S, min 1e-9, default 1e-6).",
                 MSK_DPAR_BASIS_TOL_S, 1e-9, Infinity());

  // ********** SIMPLEX OPTIONS **********
  AddMosekOption("lp:primalselection primalselection",
    "Pricing strategy of the primal simplex optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_SIM_PRIMAL_SELECTION, values_sim_selection, 0);

  AddMosekOption("lp:dualselection dualselection",
    "Pricing strategy of the dual simplex optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_SIM_DUAL_SELECTION, values_sim_selection, 0);

  AddMosekOption("lp:hotstart simhotstart",
    "Hot start of the simplex optimizers:\n"
    "\n.. value-table::\n",
    MSK_IPAR_SIM_HOTSTART, values_sim_hotstart, 1);

  AddMosekOption("lp:degen simdegen",
    "Degeneration handling of the simplex optimizers:\n"
    "\n.. value-table::\n",
    MSK_IPAR_SIM_DEGEN, values_sim_degen, 1);

  AddMosekOption("lp:solveform simsolveform",
    "Whether the simplex optimizers solve the primal or the dual problem:\n"
    "\n.. value-table::\n",
    MSK_IPAR_SIM_SOLVE_FORM, values_solve_form, 0);

  // ********** BARRIER OPTIONS **********
  AddMosekOption("bar:iterlim bariterlim",
                 "Interior-point iteration limit (default 400).",
                 MSK_IPAR_INTPNT_MAX_ITERATIONS, 0, INT_MAX);

  AddMosekOption("bar:pfeastol intpnttolpfeas",
                 "Primal feasibility tolerance of the interior-point optimizer "
                 "for linear problems (0..1, default 1e-8).",
                 MSK_DPAR_INTPNT_TOL_PFEAS, 0.0, 1.0);

  AddMosekOption("bar:dfeastol intpnttoldfeas",
                 "Dual feasibility tolerance of the interior-point optimizer "
                 "for linear problems (0..1, default 1e-8).",
                 MSK_DPAR_INTPNT_TOL_DFEAS, 0.0, 1.0);

  AddMosekOption("bar:relgaptol intpnttolrelgap",
                 "Relative gap termination tolerance of the interior-point "
                 "optimizer for linear problems (min 1e-14, default 1e-8).",
                 MSK_DPAR_INTPNT_TOL_REL_GAP, 1e-14, Infinity());

  AddMosekOption("bar:copfeastol intpntcotolpfeas",
                 "Primal feasibility tolerance of the interior-point optimizer "
                 "for conic problems (0..1, default 1e-8).",
                 MSK_DPAR_INTPNT_CO_TOL_PFEAS, 0.0, 1.0);

  AddMosekOption("bar:codfeastol intpntcotoldfeas",
                 "Dual feasibility tolerance of the interior-point optimizer "
                 "for conic problems (0..1, default 1e-8).",
                 MSK_DPAR_INTPNT_CO_TOL_DFEAS, 0.0, 1.0);

  AddMosekOption("bar:corelgaptol intpntcotolrelgap",
                 "Relative gap termination tolerance of the interior-point "
                 "optimizer for conic problems (0..1, default 1e-8).",
                 MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 0.0, 1.0);

  AddMosekOption("bar:order intpntordermethod",
    "Ordering method of the interior-point optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_INTPNT_ORDER_METHOD, values_bar_order, 0);

  AddMosekOption("bar:solveform intpntsolveform",
    "Whether the interior-point optimizer solves the primal or the dual problem:\n"
    "\n.. value-table::\n",
    MSK_IPAR_INTPNT_SOLVE_FORM, values_solve_form, 0);

  // ********** PRESOLVER OPTIONS **********
  AddMosekOption("pre:aggregate aggregate", 
    "Whether to use aggregation in presolve:\n"
    "\n.. value-table::\n", 
    MSK_IPAR_MIO_PRESOLVE_AGGREGATOR_USE, values_01_noyes_1default_, 1);

  AddMosekOption("pre:solve presolve",
    "MIP presolve:\n"
    "\n.. value-table::\n",
    MSK_IPAR_PRESOLVE_USE, alg_values_presolve_use, 2);

  AddMosekOption("pre:passes prepasses",
    "Limit on the number of presolve passes; a negative value implies MOSEK decides:\n"
    "\n"
    "  | -1 - Automatic choice (default)\n"
    "  | n>=0 - At most n passes.",
    MSK_IPAR_PRESOLVE_MAX_NUM_PASS, -1, std::numeric_limits<int>::max());

  AddMosekOption("pre:dualray_analysis dualrayanalysis",
    "Controls the amount of dual ray analysis employed by the mixed-integer optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_DUAL_RAY_ANALYSIS_LEVEL, values_mip_presolve_dual_ray, -1);

  AddMosekOption("pre:lindep lindep",
    "0/1: whether presolve removes linearly dependent rows (default 1).",
    MSK_IPAR_PRESOLVE_LINDEP_USE, 0, 1);

  AddMosekOption("pre:folding folding foldinguse",
    "Whether to use folding in presolve (for MIP problems use mip:symmetry):\n"
    "\n.. value-table::\n",
    MSK_IPAR_FOLDING_USE, values_prefoldinguse_, 2);

  AddStoredOption("pre:scale scale",
    "Whether to use scaling in presolve. Applies to both simplex and interior point method:\n"
    "\n.. value-table::\n", storedOptions_.PREScale_, values_bar_presolve_scaling);

  AddMosekOption("mip:symmetry symmetry pre:mipfolding mipfolding miosimmetrylevel",
    "Amount of symmetry detection and handling by the mixed-integer optimizer "
    "(MSK_IPAR_MIO_SYMMETRY_LEVEL):\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_SYMMETRY_LEVEL, values_premipfoldinguse_, -1);

  AddMosekOption("bar:basis bar:crossover crossover",
    "Whether the interior-point optimizer also computes an optimal basis:\n"
    "\n.. value-table::\n", MSK_IPAR_INTPNT_BASIS, values_bar_intpnt_basis, 1);


  // ********** MIP OPTIONS **********
  AddStoredOption("mip:constructsol mipconstructsol",
      "Sets MSK_IPAR_MIO_CONSTRUCT_SOL. If set to MSK_ON and all integer variables "
      "have been given a value for which a feasible mixed integer solution exists, "
      "then MOSEK generates an initial solution to the mixed integer problem by "
      "fixing all integer values and solving the remaining problem."
      "Default = OFF",
      storedOptions_.MIPConstructSol_);

  AddMosekOption("mip:gap mipgap",
                  "Max. relative MIP optimality gap (default 1e-4).",
                  MSK_DPAR_MIO_TOL_REL_GAP, 0.0, DBL_MAX);

  AddMosekOption("mip:gapabs mipgapabs",
                  "Max. absolute MIP optimality gap (default 0.0).",
                  MSK_DPAR_MIO_TOL_ABS_GAP, 0.0, DBL_MAX);

  AddMosekOption("mip:feastol feastol",
                  "Feasibility tolerance for the mixed-integer optimizer "
                  "(1e-9..1e-3, default 1e-6). Note: in this driver, 'feastol' "
                  "is this MIP tolerance; the LP tolerances are alg:feastol "
                  "(simplex) and bar:pfeastol (interior-point).",
                  MSK_DPAR_MIO_TOL_FEAS, 1e-9, 1e-3);

  AddMosekOption("mip:inttol inttol intfeastol mip:intfeastol",
                  "MIP integrality tolerance (min 1e-9, default 1e-5).",
                  MSK_DPAR_MIO_TOL_ABS_RELAX_INT, 1e-9, Infinity());

  AddMosekOption("mip:heurlevel heurlevel",
                  "MIP heuristic level to find an initial good feasible solution.\n\n"
                  "-1: automatic (default);\n0: not used;\n"
                  "positive, larger - more effort (3-5 recommended).",
                  MSK_IPAR_MIO_HEURISTIC_LEVEL, -1, INT_MAX);

  AddMosekOption("mip:feaspump feaspump",
                  "MIP feasibility pump.\n"
                  "\n.. value-table::\n",
                  MSK_IPAR_MIO_FEASPUMP_LEVEL, values_mip_feaspump_, -1);


  AddMosekOption("mip:relgapconst miorelgapconst",
    "This value is used to compute the relative gap for the solution "
    "to an integer optimization problem (min 1e-15, default 1e-10).",
    MSK_DPAR_MIO_REL_GAP_CONST, 1e-15, DBL_MAX);

  AddMosekOption("mip:conic:outapprox conicoutapprox",
                  "0*/1: If this option is turned on outer approximation is used "
                  "when solving relaxations of conic problems; "
                  "otherwise interior point is used.",
                  MSK_IPAR_MIO_CONIC_OUTER_APPROXIMATION, 0, 1);


  AddMosekOption("mip:cuts cuts mip:cutselection",
    "Cut selection level of the mixed-integer optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_CUT_SELECTION_LEVEL, values_mio_cutselection, -1);

  for (const auto& family: mosek_cut_families)
    AddMosekOption(family.names,
      fmt::format("0/1: whether to generate {}", family.descr).c_str(),
      family.key, 0, 1);

  AddMosekOption("mip:nodeselect nodeselect mip:nodeselection",
    "Node selection strategy of the mixed-integer optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_NODE_SELECTION, values_mio_nodeselect, 0);

  AddMosekOption("mip:branchdir branchdir",
    "Branching direction of the mixed-integer optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_BRANCH_DIR, values_mio_branchdir, 0);

  AddMosekOption("mip:probing probing",
    "Amount of probing in the presolve of the mixed-integer optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_PROBING_LEVEL, values_mio_probing, -1);

  AddMosekOption("mip:rootmethod rootmethod mip:rootoptimizer",
    "Optimizer for the root relaxation of a mixed-integer problem:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_ROOT_OPTIMIZER, values_mio_optimizer, 2);

  AddMosekOption("mip:nodemethod nodemethod mip:nodeoptimizer",
    "Optimizer for the node relaxations of a mixed-integer problem:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_NODE_OPTIMIZER, values_mio_optimizer, 2);

  AddMosekOption("mip:rinsmaxnodes rinsmaxnodes",
    "Maximal number of nodes in each RINS heuristic run "
    "(-1: automatic, default; 0: RINS off).",
    MSK_IPAR_MIO_RINS_MAX_NODES, -1, INT_MAX);

  AddMosekOption("mip:restarts maxrestarts",
    "Maximal number of restarts of the mixed-integer optimizer (default 10).",
    MSK_IPAR_MIO_MAX_NUM_RESTARTS, 0, INT_MAX);

  AddMosekOption("mip:rootcutrounds rootcutrounds",
    "Maximal number of cut separation rounds at the root node (default 100).",
    MSK_IPAR_MIO_MAX_NUM_ROOT_CUT_ROUNDS, 0, INT_MAX);

  AddMosekOption("mip:djcbigm djcmaxbigm",
    "Maximal big-M value used when disjunctive constraints are "
    "reformulated by the mixed-integer optimizer (default 1e6).",
    MSK_DPAR_MIO_DJC_MAX_BIGM, 0.0, DBL_MAX);

  AddMosekOption("mip:varselection varselection",
    "Controls the variable selection strategy employed by the mixed-integer optimizer:\n"
    "\n.. value-table::\n",
    MSK_IPAR_MIO_VAR_SELECTION, values_mip_var_selection_types, 2);



  // *********** Tech options *********** 

  AddListOption("tech:optionnative optionnative optnative tech:param",
    "General way to specify values of both documented and "
    "undocumented Mosek parameters; value should be a quoted "
    "string (delimited by ' or \") containing a parameter name, a "
    "space, and the value to be assigned to the parameter.  Can "
    "appear more than once.  Cannot be used to query current "
    "parameter values.",
    storedOptions_.inlineParams_);

  AddStoredOption("tech:optionnativeread optionnativeread tech:param:read param:read",
    "Name of Mosek parameter file (surrounded by 'single' or "
    "\"double\" quotes if the name contains blanks) to be read. File format:\n"
    "\n"
    "  BEGIN MOSEK\n"
    "  MSK_DPAR_MIO_MAX_TIME 3\n"
    "  END MOSEK\n"
    "\n"
    "Parameter descriptions: docs.mosek.com/latest/cmdtools/param-groups.html.",
    storedOptions_.paramRead_);

  AddStoredOption("tech:optionnativewrite optionnativewrite tech:param:write param:write",
    "Name of Mosek parameter file (surrounded by 'single' or \"double\" quotes if the "
    "name contains blanks) to be written.",
    storedOptions_.paramWrite_);


  AddMosekOption("tech:seed seed tech:mipseed mipseed",
                  "Random number seed (default 42), "
                  "used for randomization in the mixed-integer optimizer, "
                  "may influence the solution path.",
                  MSK_IPAR_MIO_SEED, 0, INT_MAX);

  AddMosekOption("tech:lpseed lpseed",
                  "Sets the random seed (0..32749, default 23456) used for "
                  "randomization in the simplex optimizers.",
                  MSK_IPAR_SIM_SEED, 0, 32749);

  AddMosekOption("tech:threads threads",
    "Controls the number of threads employed by the optimizer. "
    "Default 0 ==> number of threads used will be equal to the number "
    "of cores detected on the machine.",
    MSK_IPAR_NUM_THREADS, 0, INT_MAX);

  AddMosekOption("tech:outlev outlev",
    "0*/1: Whether to write MOSEK log lines to stdout (MSK_IPAR_LOG; "
    "larger levels can be set with tech:optionnative).",
    MSK_IPAR_LOG, 0, INT_MAX);

  AddStoredOption("tech:logfile logfile",
    "Log file name. Without tech:outlev=1 the log goes to the file only.",
    storedOptions_.logFile_);
  // *********** Tech options end *********** 
  





  AddMosekOption("lim:time timelim timelimit time_limit",
    "Limit on solve time (in seconds; default -1: no limit).",
    MSK_DPAR_OPTIMIZER_MAX_TIME, -1.0, DBL_MAX);

  AddMosekOption("lim:iter iterlim simmaxiterations",
    "Simplex iteration limit (default 10000000).",
    MSK_IPAR_SIM_MAX_ITERATIONS, 0, INT_MAX);

  AddMosekOption("lim:miptime miptimelim",
    "Time limit for the mixed-integer optimizer in seconds "
    "(MSK_DPAR_MIO_MAX_TIME; default -1: no limit; lim:time applies too).",
    MSK_DPAR_MIO_MAX_TIME, -1.0, DBL_MAX);

  AddMosekOption("lim:lowerobjcut lowerobjcut",
    "Lower objective cut: MOSEK stops when it proves that the optimal "
    "objective value is below this (default -infinity). "
    "With lim:upperobjcut this is MOSEK's objective range, reported as "
    "'objective range' (solve_result 406).",
    MSK_DPAR_LOWER_OBJ_CUT, -DBL_MAX, DBL_MAX);

  AddMosekOption("lim:upperobjcut upperobjcut",
    "Upper objective cut: MOSEK stops when it proves that the optimal "
    "objective value is above this (default +infinity).",
    MSK_DPAR_UPPER_OBJ_CUT, -DBL_MAX, DBL_MAX);

  AddMosekOption("lim:nodes nodelim nodelimit",
    "Limit on the number of MIP branches (MSK_IPAR_MIO_MAX_NUM_BRANCHES; "
    "default -1: no limit).",
    MSK_IPAR_MIO_MAX_NUM_BRANCHES, -1, INT_MAX);

  AddMosekOption("lim:sol sollim sollimit solutionlimit",
    "Limit on the number of feasible MIP solutions found, causing early "
    "termination if exceeded; default -1 (no limit).",
    MSK_IPAR_MIO_MAX_NUM_SOLUTIONS, -1, INT_MAX);


}

double MosekBackend::MIPGap() {
  return getDblAttr(MSK_DINF_MIO_OBJ_REL_GAP);
}

double MosekBackend::BestDualBound() {
  return getDblAttr(MSK_DINF_MIO_OBJ_BOUND);
}

double MosekBackend::MIPGapAbs() {
  return getDblAttr(MSK_DINF_MIO_OBJ_ABS_GAP);
}

SensRangesPresolved MosekBackend::GetSensRangesPresolved()
{
  if (IsMIP() && !HasFixedLP()) return {};
  MSKtask_t task = HasFixedLP() ? fixed_task_ : lp();
  // Set sensitivity range parameters for constraints
  int lencon = NumLinCons();
  std::vector<MSKint32t> cindex(lencon);
  std::vector<MSKmarke> cmarklb(lencon);
  std::vector<MSKmarke> cmarkub(lencon);
  for (int i=0; i<lencon; i++)
  {
    cindex[i] = i;
    cmarklb[i] = MSK_MARK_LO;
    cmarkub[i] = MSK_MARK_UP;
  }
  std::vector<MSKrealt> cpricelblo(lencon);
  std::vector<MSKrealt> cpricelbhi(lencon);
  std::vector<MSKrealt> cpriceublo(lencon);
  std::vector<MSKrealt> cpriceubhi(lencon);
  std::vector<MSKrealt> crangelblo(lencon);
  std::vector<MSKrealt> crangelbhi(lencon);
  std::vector<MSKrealt> crangeublo(lencon);
  std::vector<MSKrealt> crangeubhi(lencon);

  // Set sensitivity range parameters for variables
  int lenvar = NumVars();
  std::vector<MSKint32t> vindex(lenvar);
  std::vector<MSKmarke> vmarklb(lenvar);
  std::vector<MSKmarke> vmarkub(lenvar);
  for (int i=0; i<lenvar; i++)
  {
    vindex[i] = i;
    vmarklb[i] = MSK_MARK_LO;
    vmarkub[i] = MSK_MARK_UP;
  }
  std::vector<MSKrealt> vpricelblo(lenvar);
  std::vector<MSKrealt> vpricelbhi(lenvar);
  std::vector<MSKrealt> vpriceublo(lenvar);
  std::vector<MSKrealt> vpriceubhi(lenvar);
  std::vector<MSKrealt> vrangelblo(lenvar);
  std::vector<MSKrealt> vrangelbhi(lenvar);
  std::vector<MSKrealt> vrangeublo(lenvar);
  std::vector<MSKrealt> vrangeubhi(lenvar);

  std::vector<MSKrealt> opricelo(lenvar);
  std::vector<MSKrealt> opricehi(lenvar);
  std::vector<MSKrealt> orangelo(lenvar);
  std::vector<MSKrealt> orangehi(lenvar);

  // Compute sensitivity ranges for lower bounds
  MOSEK_CCALL(
    MSK_primalsensitivity(
        task,
        lencon,
        cindex.data(),
        cmarklb.data(),
        lenvar,
        vindex.data(),
        vmarklb.data(),
        cpricelblo.data(),
        cpricelbhi.data(),
        crangelblo.data(),
        crangelbhi.data(),
        vpricelblo.data(),
        vpricelbhi.data(),
        vrangelblo.data(),
        vrangelbhi.data()
    )
  );

  // Compute sensitivity ranges for upper bounds
  MOSEK_CCALL(
    MSK_primalsensitivity(
        task,
        lencon,
        cindex.data(),
        cmarkub.data(),
        lenvar,
        vindex.data(),
        vmarkub.data(),
        cpriceublo.data(),
        cpriceubhi.data(),
        crangeublo.data(),
        crangeubhi.data(),
        vpriceublo.data(),
        vpriceubhi.data(),
        vrangeublo.data(),
        vrangeubhi.data()
    )
  );

  // Compute sensitivity ranges for objective coeffs
  MOSEK_CCALL(
    MSK_dualsensitivity(
        task,
        lenvar,
        vindex.data(),
        opricelo.data(),
        opricehi.data(),
        orangelo.data(),
        orangehi.data()
    )
  );

  // Get bounds and add them to sensitivity ranges
  std::vector<MSKboundkeye> vbk(lenvar);
  std::vector<MSKrealt> vbl(lenvar), vbu(lenvar);
  MOSEK_CCALL(MSK_getvarboundslice(task, (MSKint32t)0, (MSKint32t)lenvar, vbk.data(), vbl.data(), vbu.data()));
  std::transform(vrangelblo.begin(), vrangelblo.end(), vbl.begin(), vrangelblo.begin(), std::plus<MSKrealt>());
  std::transform(vrangelbhi.begin(), vrangelbhi.end(), vbl.begin(), vrangelbhi.begin(), std::plus<MSKrealt>());
  std::transform(vrangeublo.begin(), vrangeublo.end(), vbu.begin(), vrangeublo.begin(), std::plus<MSKrealt>());
  std::transform(vrangeubhi.begin(), vrangeubhi.end(), vbu.begin(), vrangeubhi.begin(), std::plus<MSKrealt>());
  for (int i=0; i<lenvar; i++)
  {
    if (vbk[i] == MSK_BK_UP || vbk[i] == MSK_BK_FR)
    {
      vrangelblo[i] = -MSK_INFINITY;
      vrangelbhi[i] = MSK_INFINITY;
    }
    if (vbk[i] == MSK_BK_LO || vbk[i] == MSK_BK_FR)
    {
      vrangeublo[i] = -MSK_INFINITY;
      vrangeubhi[i] = MSK_INFINITY;
    }
  }

  std::vector<MSKboundkeye> cbk(lencon);
  std::vector<MSKrealt> cbl(lencon), cbu(lencon);
  MOSEK_CCALL(MSK_getconboundslice(task, (MSKint32t)0, (MSKint32t)lencon, cbk.data(), cbl.data(), cbu.data()));
  std::transform(crangelblo.begin(), crangelblo.end(), cbl.begin(), crangelblo.begin(), std::plus<MSKrealt>());
  std::transform(crangelbhi.begin(), crangelbhi.end(), cbl.begin(), crangelbhi.begin(), std::plus<MSKrealt>());
  std::transform(crangeublo.begin(), crangeublo.end(), cbu.begin(), crangeublo.begin(), std::plus<MSKrealt>());
  std::transform(crangeubhi.begin(), crangeubhi.end(), cbu.begin(), crangeubhi.begin(), std::plus<MSKrealt>());
  std::vector<MSKrealt> rhslo(lencon, 0.0), rhshi(lencon, 0.0);
  for (int i=0; i<lencon; i++)
  {
    if (cbk[i] == MSK_BK_UP || cbk[i] == MSK_BK_FR)
    {
      crangelblo[i] = -MSK_INFINITY;
      crangelbhi[i] = MSK_INFINITY;
      rhslo[i] = crangeublo[i];
      rhshi[i] = crangeubhi[i];
    }
    if (cbk[i] == MSK_BK_LO || cbk[i] == MSK_BK_FR)
    {
      crangeublo[i] = -MSK_INFINITY;
      crangeubhi[i] = MSK_INFINITY;
      rhslo[i] = crangelblo[i];
      rhshi[i] = crangelbhi[i];
    }
  }

  std::vector<MSKrealt> c(lenvar);
  MOSEK_CCALL(MSK_getc(task, c.data()));
  std::transform(orangelo.begin(), orangelo.end(), c.begin(), orangelo.begin(), std::plus<MSKrealt>());
  std::transform(orangehi.begin(), orangehi.end(), c.begin(), orangehi.begin(), std::plus<MSKrealt>());

  /// Return sensitivity ranges
  SensRangesPresolved sensr;

  sensr.varlbhi = { {vrangelbhi} };
  sensr.varlblo = { {vrangelblo} };
  sensr.varubhi = { {vrangeubhi} };
  sensr.varublo = { {vrangeublo} };
  sensr.varobjhi = { {orangehi} };
  sensr.varobj = { {c} };
  sensr.varobjlo = { {orangelo} };
	sensr.conlbhi = { {}, {{{CG_Algebraic, crangelbhi}}} };
	sensr.conlblo = { {}, {{{CG_Algebraic, crangelblo}}} };
	sensr.conubhi = { {}, {{{CG_Algebraic, crangeubhi}}} };
	sensr.conublo = { {}, {{{CG_Algebraic, crangeublo}}} };
  sensr.conrhshi = { {}, {{{CG_Algebraic, rhshi}}} };
  sensr.conrhslo = { {}, {{{CG_Algebraic, rhslo}}} };

  return sensr;
}

ArrayRef<int> MosekBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  MOSEK_CCALL(MSK_getskx(lp(), solToFetch_, (MSKstakeye *)vars.data()));

  for (auto& s : vars)
  {
    switch (s)
    {
      case MSK_SK_UNK:
        s = (int)BasicStatus::none;
        break;
      case MSK_SK_BAS:
        s = (int)BasicStatus::bas;
        break;
      case MSK_SK_LOW:
        s = (int)BasicStatus::low;
        break;
      case MSK_SK_UPR:
        s = (int)BasicStatus::upp;
        break;
      case MSK_SK_SUPBAS:
        s = (int)BasicStatus::sup;
        break;
      case MSK_SK_FIX:
        s = (int)BasicStatus::equ;
        break;
      default:
        MP_RAISE(fmt::format("Unknown Mosek VBasis value: {}", s));
    }
  }

  return vars;
}

ArrayRef<int> MosekBackend::ConStatii() {
  std::vector<int> cons(NumLinCons());
  MOSEK_CCALL(MSK_getskc(lp(), solToFetch_, (MSKstakeye *)cons.data()));

  for (auto& s : cons)
  {
    switch (s)
    {
      case MSK_SK_UNK:
        s = (int)BasicStatus::none;
        break;
      case MSK_SK_BAS:
        s = (int)BasicStatus::bas;
        break;
      case MSK_SK_LOW:
        s = (int)BasicStatus::low;
        break;
      case MSK_SK_UPR:
        s = (int)BasicStatus::upp;
        break;
      case MSK_SK_SUPBAS:
        s = (int)BasicStatus::sup;
        break;
      case MSK_SK_FIX:
        s = (int)BasicStatus::equ;
        break;
      default:
        MP_RAISE(fmt::format("Unknown Mosek VBasis value: {}", s));
    }
  }
  return cons;
}

void MosekBackend::VarStatii(ArrayRef<int> vst) {
  // BasicStatus enum values: from 0 to 6: none, bas, sup, low, upp, equ, btw
  std::vector<MSKstakeye> skx(vst.size());

  // TODO: remove this if we are sure that it is the same.
  MSKint32t numvar;
  MOSEK_CCALL(MSK_getnumvar(lp(), &numvar));
  assert(numvar == vst.size());

  for (auto j = 0; j < skx.size(); j++)
  {
    switch ((BasicStatus)vst[j])
    {
      case BasicStatus::bas:
        skx[j] = MSK_SK_BAS;
        break;
      case BasicStatus::low:
        skx[j] = MSK_SK_LOW;
        break;
      case BasicStatus::equ:
        skx[j] = MSK_SK_FIX;
        break;
      case BasicStatus::upp:
        skx[j] = MSK_SK_UPR;
        break;
      case BasicStatus::sup:
      case BasicStatus::btw:
        skx[j] = MSK_SK_SUPBAS;
        break;
      case BasicStatus::none:
        /// 'none' is assigned to new variables (added after solving).
        /// Compute low/upp/fix/sup, depending on where initial value 0.0 is between the bounds
        MSKboundkeye bk;
        MSKrealt bl, bu, tolx;
        MOSEK_CCALL(MSK_getvarbound(lp(), (MSKint32t)j, &bk, &bl, &bu));
        MOSEK_CCALL(MSK_getdouparam(lp(), MSK_DPAR_DATA_TOL_X, &tolx));

        if (0.0 < bl + tolx) skx[j] = MSK_SK_LOW;
        else if (bu - tolx < 0.0) skx[j] = MSK_SK_UPR;
        else if (bu - bl < tolx) skx[j] = MSK_SK_FIX;
        else skx[j] = MSK_SK_SUPBAS;

        break;
      default:
        MP_RAISE(fmt::format("Unknown AMPL var status value: {}", skx[j]));
    }
  }

  // Input the variable status keys
  MOSEK_CCALL(MSK_putskx(lp(), solToFetch_, (MSKstakeye *)skx.data()));
}

void MosekBackend::ConStatii(ArrayRef<int> cst) {
  MSKint32t numcon;
  MOSEK_CCALL(MSK_getnumcon(lp(), &numcon));
	std::vector<int> skc(cst.data(), cst.data() + numcon);

  // The status key of a constraint is the status key of the logical (slack) variable assigned to it.
  // The slack is defined as: l <= a'x <= u rewritten as a'x - s = 0, l <= s <= u.
	for (size_t j = 0; j < skc.size(); j++)
  {
    switch ((BasicStatus)skc[j])
    {
      case BasicStatus::bas:
        skc[j] = MSK_SK_BAS;
        break;
      case BasicStatus::sup:
      case BasicStatus::btw:
        skc[j] = MSK_SK_SUPBAS;
        break;
      case BasicStatus::low:
        skc[j] = MSK_SK_LOW;
        break;
      case BasicStatus::upp:
        skc[j] = MSK_SK_UPR;
        break;
      case BasicStatus::equ:
        skc[j] = MSK_SK_FIX;
        break;
      case BasicStatus::none:
        /// 'none' is assigned to new constraints (added after solving).
        /// Compute low/upp/fix/sup, depending on where initial value 0.0 is between the bounds
        MSKboundkeye bk;
        MSKrealt bl, bu, tolx;
        MOSEK_CCALL(MSK_getconbound(lp(), (MSKint32t)j, &bk, &bl, &bu));
        // MSK_DPAR_DATA_TOL_X is both for variables and constraints
        MOSEK_CCALL(MSK_getdouparam(lp(), MSK_DPAR_DATA_TOL_X, &tolx));

        if (0.0 < bl + tolx) skc[j] = MSK_SK_LOW;
        else if (bu - tolx < 0.0) skc[j] = MSK_SK_UPR;
        else if (bu - bl < tolx) skc[j] = MSK_SK_FIX;
        else skc[j] = MSK_SK_SUPBAS;

        break;
      default:
        MP_RAISE(fmt::format("Unknown AMPL var status value: {}", skc[j]));
    }
  }

  // Input the variable status keys
  MOSEK_CCALL(MSK_putskc(lp(), solToFetch_, (MSKstakeye *)skc.data()));
}

SolutionBasis MosekBackend::GetBasis() {
  if (IsMIP()) return FixedLPBasis(CG_Algebraic);
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
				{{{ CG_Algebraic, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());       // not for constraints
  }
  return { std::move(varstt), std::move(constt) };
}

void MosekBackend::SetBasis(SolutionBasis basis) {
  solToFetch_ = GetSolutionTypeForStart();
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
	auto constt = mv.GetConValues()(CG_Algebraic);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

void MosekBackend::AddPrimalDualStart(Solution sol)
{
  solToFetch_ = GetSolutionTypeForStart();
  auto mv = GetValuePresolver().PresolveSolution(
        { sol.primal, sol.dual } );
  auto x0 = mv.GetVarValues()();
	auto pi0 = mv.GetConValues()(CG_Algebraic);
  MOSEK_CCALL(MSK_putxx(lp(), solToFetch_, (MSKrealt *)x0.data()));
  MOSEK_CCALL(MSK_puty(lp(), solToFetch_, (MSKrealt *)pi0.data()));
}

void MosekBackend::AddMIPStart(
    ArrayRef<double> x0_unpres, ArrayRef<int> sparsity)
{
  auto mv = GetValuePresolver().PresolveSolution( { x0_unpres } );
  auto x0 = mv.GetVarValues()();

  if (Mosek_mip_construct_sol())
  {
    // Use integer part of given solution to construct initial solution by solving for continuous part.
    MOSEK_CCALL(MSK_putintparam(lp(), MSK_IPAR_MIO_CONSTRUCT_SOL, MSK_ON));
  }
  MOSEK_CCALL(MSK_putxx(lp(), MSK_SOL_ITG, (MSKrealt *)x0.data()));
}

ArrayRef<double> MosekBackend::Ray()
{
  // Problem is checked to be dual infeasible at this point, so the ray is the primal solution variable.
  // (Primal can be unbounded or infeasible.)
  std::vector<double> xx(NumVars());
  MOSEK_CCALL(MSK_getxx(lp(), solToFetch_, (MSKrealt *)xx.data()));
  // Argument is a ModelValues<ValueMap>, which can be constructed from a ValueMap for only the variables.
  auto mv = GetValuePresolver().PostsolveSolution({xx});
  // We get the ValueMap for variables, and by calling that, the value itself.
  // (Similar to .MoveOut(0) for single key maps)
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> MosekBackend::DRay()
{
  // Problem is checked to be (primal) infeasible at this point, so the ray is the dual solution variable.
  // (Dual can be unbounded or infeasible.)
  std::vector<double> y(NumLinCons());
  if (!IsMIP())
    MOSEK_CCALL(MSK_gety(lp(), solToFetch_, (MSKrealt *)y.data()));
  // Argument is a ModelValues<ValueMap>, which is constructed now from two ValueMaps, an empty for variables, and
  // a second one for constraints. The constraints ValueMap is given as a std::map for only linear constraints
	auto mv = GetValuePresolver().
			PostsolveSolution({
													{},
													{ { {CG_Algebraic, std::move(y)} } }
												});
  // We get the ValueMap for constraints, and by calling MoveOut(), the value itself.
  // (Similar to .MoveOut(0) for single key maps)
  return mv.GetConValues().MoveOut();
}

} // namespace mp


// AMPLs

AMPLS_MP_Solver* Open_mosek(CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::MosekBackend()},
        cb);
}

void AMPLSClose_mosek(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_mosek(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::MosekBackend*>(AMPLSGetBackend(slv))->lp();
}


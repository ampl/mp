#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "scipbackend.h"

extern "C" {
  #include "scip-ampls-c-api.h"    // Scip AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptScip(void* scip) {
  return SCIPinterruptSolve((SCIP*)scip);
  //return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateScipBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::ScipBackend()};
}


namespace mp {

/// Create Scip Model Manager
/// @param gc: the Scip common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return ScipModelMgr
std::unique_ptr<BasicModelManager>
CreateScipModelMgr(ScipCommon&, Env&, pre::BasicValuePresolver*&);


ScipBackend::ScipBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateScipModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

ScipBackend::~ScipBackend() {
  CloseSolver();
}


const char* ScipBackend::GetBackendName()
  { return "SCIPBackend"; }

std::string ScipBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", SCIPmajorVersion(), SCIPminorVersion(), 
  SCIPtechVersion());
}


bool ScipBackend::IsMIP() const {
  // TODO
  //return getIntAttr(Solver::VARS_INT) > 0;
  //return getIntAttr(SCIP_INTATTR_ISMIP);
  return true;
}

bool ScipBackend::IsQCP() const {
  //return getIntAttr(Solver::CONS_QUAD) > 0;
// return getIntAttr(SCIP_INTATTR_QELEMS) > 0;
  return false;
}

ArrayRef<double> ScipBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  /*
  if (IsMIP()) 
    error = SCIP_GetSolution(lp(), x.data());
  else
    error = SCIP_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);
  if (error)
    x.clear();
    */
  return x;
}

pre::ValueMapDbl ScipBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> ScipBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
 // int error = SCIP_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 0;
  if (error)
    pi.clear();
  return pi;
}

double ScipBackend::ObjectiveValue() const {
  return SCIPgetPrimalbound(getSCIP());
}

double ScipBackend::NodeCount() const {
  return SCIPgetNNodes(getSCIP());
}

double ScipBackend::SimplexIterations() const {
  return SCIPgetNPrimalLPIterations(getSCIP()) + SCIPgetNDualLPIterations(getSCIP());
}

int ScipBackend::BarrierIterations() const {
  return SCIPgetNBarrierLPIterations(getSCIP());
}

void ScipBackend::ExportModel(const std::string &file) {
  SCIP_CCALL( SCIPwriteOrigProblem(getSCIP(), file.data(), NULL, FALSE) );
}


void ScipBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptScip, getSCIP());
  // TODO Check interrupter
  //SCIP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void ScipBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  SCIP_CCALL( SCIPsolve(getSCIP()) );
  WindupSCIPSolve();
}

void ScipBackend::WindupSCIPSolve() { }

void ScipBackend::ReportResults() {
  ReportSCIPResults();
  BaseBackend::ReportResults();
}

void ScipBackend::ReportSCIPResults() {
  SetStatus( ConvertSCIPStatus() );
  AddSCIPMessages();
  if (need_multiple_solutions())
    ReportSCIPPool();
}
std::vector<double> ScipBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // SCIP_CCALL(SCIP_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double ScipBackend::getPoolObjective(int i)
{
  double obj;
 // SCIP_CCALL(SCIP_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void ScipBackend::ReportSCIPPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(SCIP_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void ScipBackend::AddSCIPMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (!IsContinuous()) {
    auto nnd = NodeCount();
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
  }
}

std::pair<int, std::string> ScipBackend::ConvertSCIPStatus() {
  namespace sol = mp::sol;
  SCIP_STATUS status = SCIPgetStatus(getSCIP());
  switch (status) {
    case SCIP_STATUS_UNKNOWN:
      return { sol::UNKNOWN, "solving status not yet known" };
    case SCIP_STATUS_USERINTERRUPT:
      return { sol::INTERRUPTED, "unfinished" };
    case SCIP_STATUS_NODELIMIT:
      return { sol::LIMIT, "node limit reached" };
    case SCIP_STATUS_TOTALNODELIMIT:
      return { sol::LIMIT, "total node limit reached (incl. restarts)" };
    case SCIP_STATUS_STALLNODELIMIT:
      return { sol::LIMIT, "stalling node limit reached (no inprovement w.r.t. primal bound)" };
    case SCIP_STATUS_TIMELIMIT:
      return { sol::LIMIT, "time limit reached" };
    case SCIP_STATUS_MEMLIMIT:
      return { sol::LIMIT, "memory limit reached" };
    case SCIP_STATUS_GAPLIMIT:
      return { sol::LIMIT, "gap limit reached" };
    case SCIP_STATUS_SOLLIMIT:
      return { sol::LIMIT, "solution limit reached" };
    case SCIP_STATUS_BESTSOLLIMIT:
      return { sol::LIMIT, "solution improvement limit reached" };
    case SCIP_STATUS_RESTARTLIMIT:
      return { sol::LIMIT, "restart limit was reached" };
    case SCIP_STATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case SCIP_STATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case SCIP_STATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case SCIP_STATUS_INFORUNBD:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case SCIP_STATUS_TERMINATE:
      return { sol::FAILURE, "process received a SIGTERM signal" };
    }
  return { sol::UNKNOWN, "not solved" };
}


void ScipBackend::FinishOptionParsing() {
  int v=-1;
  GetSolverOption("display/verblevel", v);
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo lp_values_method[] = {
  { "-1", "Automatic (default)", -1},
  { "1", "Dual simplex", 1},
  { "2", "Barrier", 2},
  { "3", "Crossover", 3},
  { "4", "Concurrent (simplex and barrier simultaneously)", 4},
};


static const mp::OptionValueInfo alg_values_level[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Fast", 1},
  { "2", "Normal", 2},
  { "3", "Aggressive", 3}
}; 

static const mp::OptionValueInfo lp_dualprices_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Use Devex pricing algorithm", 0},
  { "1", "Using dual steepest-edge pricing algorithm", 1}
};

static const mp::OptionValueInfo lp_barorder_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Approximate Minimum Degree (AMD)", 0},
  { "1", "Nested Dissection (ND)", 1}
};

void ScipBackend::InitCustomOptions() {

  set_option_header(
      "SCIP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``scip_options``. For example::\n"
      "\n"
      "  ampl: option scip_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  AddSolverOption("lim:time timelim timelimit time_limit",
      "Limit on solve time (in seconds; default: 1e+20).",
      "limits/time", 0.0, SCIP_REAL_MAX);

  AddSolverOption("tech:threads threads",
     "How many threads to use when using the barrier algorithm "
     "or solving MIP problems; default 0 ==> automatic choice.",
     "lp/advanced/threads", 0, 128);

  AddStoredOption("tech:logfile logfile",
    "Log file name; note that the solver log will be written to the log "
    "regardless of the value of tech:outlev.",
    storedOptions_.logFile_);
}


double ScipBackend::MIPGap() {
  return SCIPgetGap(getSCIP())<Infinity() ? SCIPgetGap(getSCIP()) : AMPLInf();
}
double ScipBackend::BestDualBound() {
  return SCIPgetDualbound(getSCIP());
}

double ScipBackend::MIPGapAbs() {
  double gapabs = std::fabs(ObjectiveValue() - BestDualBound());
  return gapabs<Infinity() ? gapabs : AMPLInf();
}


ArrayRef<int> ScipBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  SCIP_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case SCIP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case SCIP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case SCIP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case SCIP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case SCIP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Scip VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> ScipBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  SCIP_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case SCIP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case SCIP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case SCIP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case SCIP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case SCIP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Scip VBasis value: {}", s));
    }
  }*/
  return cons;
}

void ScipBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = SCIP_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = SCIP_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = SCIP_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = SCIP_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = SCIP_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!SCIP_GetColInfo(lp(), SCIP_DBLINFO_LB, 1, index, &lb) && 
        !SCIP_GetColInfo(lp(), SCIP_DBLINFO_UB, 1, index, &ub))
      { 
        if (lb >= -1e-6)
          s = -1;
        else if (ub <= 1e-6)
          s = -2;
        else
          s = -3;  // or, leave at 0?
      }
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }
  SCIP_SetBasis(lp(), stt.data(), NULL);
  */
}

void ScipBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = SCIP_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = SCIP_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  SCIP_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis ScipBackend::GetBasis() {
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

void ScipBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void ScipBackend::AddMIPStart(ArrayRef<double> x0) {
  //SCIP_CCALL(SCIP_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));

}


} // namespace mp


// AMPLs
void* AMPLSOpenScip(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::ScipBackend()},
    slv_opt, cb);
}

void AMPLSCloseScip(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetScipmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::ScipBackend*>(AMPLSGetBackend(slv))->getSCIP();
}

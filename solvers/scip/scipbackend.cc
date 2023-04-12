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
  SCIPinterruptSolve(static_cast<SCIP*>(scip));
  return true;
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
  SCIP* scip = getSCIP();
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  for (int i = 0; i < num_vars; i++)
    x[i] = SCIPgetSolVal(scip, SCIPgetBestSol(scip), getPROBDATA()->vars[i]);
  return x;
}

pre::ValueMapDbl ScipBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> ScipBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
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
  SCIP* scip = getSCIP();
  int num_vars = NumVars();
  std::vector<double> vars(num_vars);
  for (int j = 0; j < num_vars; j++)
    vars[j] = SCIPgetSolVal(scip, SCIPgetSols(scip)[i], getPROBDATA()->vars[j]);
  return vars;
}
double ScipBackend::getPoolObjective(int i)
{
  assert(i < SCIPgetNSols(getSCIP()));
  double obj;
  obj = SCIPgetSolOrigObj(getSCIP(), SCIPgetSols(getSCIP())[i]);
  return obj;
}
void ScipBackend::ReportSCIPPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions = SCIPgetNSols(getSCIP());
  
  while (++iPoolSolution < nsolutions) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
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
  { "s", "automatic simplex (default)", -1},
  { "p", "primal simplex", 0},
  { "d", "dual simplex", 1},
  { "b", "barrier", 2},
  { "c", "barrier with crossover", 3},
};

static const mp::OptionValueInfo values_pricing[] = {
  { "l", "lpi default (default)", -1},
  { "a", "auto", 0},
  { "f", "full pricing", 1},
  { "p", "partial", 2},
  { "s", "steepest edge pricing", 3},
  { "q", "quickstart steepest edge pricing", 4},
  { "d", "devex pricing", 5}
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

  AddSolverOption("tech:outlev outlev",
    "0*/1/2/3/4/5: Whether to write SCIP log lines (chatter) to stdout and to file.",
    "display/verblevel", 0, 5);

  AddStoredOption("tech:exportfile writeprob writemodel",
    "Specifies the name of a file where to export the model before "
    "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
    "Default = \"\" (don't export the model).",
    storedOptions_.exportFile_);

  AddSolverOption("lim:time timelim timelimit time_limit",
    "Limit on solve time (in seconds; default: 1e+20).",
    "limits/time", 0.0, SCIP_REAL_MAX);

  AddStoredOption("tech:logfile logfile",
    "Log file name; note that the solver log will be written to the log "
    "regardless of the value of tech:outlev.",
    storedOptions_.logFile_);



  AddSolverOption("alg:method method lpmethod",
    "LP algorithm for solving initial LP relaxations:\n"
    "\n.. value-table::\n", "lp/initalgorithm", lp_values_method, "s");

  AddSolverOption("alg:remethod remethod relpmethod",
    "LP algorithm for resolving LP relaxations if a starting basis exists:\n"
    "\n.. value-table::\n", "lp/resolvealgorithm", lp_values_method, "s");

  ////////////////////////// LP //////////////////////////
  AddSolverOption("lp:pricing pricing",
    "Pricing strategy:\n"
    "\n.. value-table::",
    "lp/pricing", values_pricing, "l");
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

void ScipBackend::AddPrimalDualStart(Solution sol)
{
  auto mv = GetValuePresolver().PresolveSolution(
        { sol.primal, sol.dual } );
  auto x0 = mv.GetVarValues()();
	auto pi0 = mv.GetConValues()(CG_Linear);
  SCIP_SOL* solution;
  SCIP_Bool keep;
  SCIP_CCALL( SCIPcreateSol(getSCIP(), &solution, NULL) );

  SCIP_CCALL( SCIPsetSolVals(getSCIP(), solution, getPROBDATA()->nvars, getPROBDATA()->vars, x0.data()) );

  SCIP_CCALL( SCIPaddSolFree(getSCIP(), &solution, &keep) );
}

void ScipBackend::AddMIPStart(ArrayRef<double> x0) {
  SCIP_SOL* solution;
  SCIP_Bool keep;
  SCIP_CCALL( SCIPcreateSol(getSCIP(), &solution, NULL) );

  SCIP_CCALL( SCIPsetSolVals(getSCIP(), solution, getPROBDATA()->nvars, getPROBDATA()->vars, (double*)x0.data()) );

  SCIP_CCALL( SCIPaddSolFree(getSCIP(), &solution, &keep) );
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

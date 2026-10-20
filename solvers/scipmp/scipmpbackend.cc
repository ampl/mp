#include <algorithm>
#include <cstdio>
#include <memory>
#include <sstream>
#include <unordered_map>
#include "lpi/lpi.h"
#include "scip/cons_linear.h"
#include "scip/scip_iisfinder.h"
#include "scip/pub_iisfinder.h"
#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "scipmpbackend.h"

extern "C" {
  #include "scipmp-ampls-c-api.h"    // Scip AMPLS C API
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

std::string ScipBackend::set_external_libs() {
  std::string ext_libs = "";

  for (int i = 0; i < SCIPgetNExternalCodes(getSCIP()); i++)
    ext_libs = ext_libs + std::string(SCIPgetExternalCodeNames(getSCIP())[i]) + " " + std::string(SCIPgetExternalCodeDescriptions(getSCIP())[i]) + "\n";

  return ext_libs;
}

bool ScipBackend::IsMIP() const {
  return SCIPgetNOrigIntVars(getSCIP()) > 0 || SCIPgetNOrigBinVars(getSCIP()) > 0;
}

bool ScipBackend::IsQCP() const {
  for (int i = 0; i < SCIPgetNOrigConss(getSCIP()); i++) {
    SCIP_Bool isquadratic;
    SCIP_CCALL( SCIPcheckQuadraticNonlinear(getSCIP(), SCIPgetOrigConss(getSCIP())[i], &isquadratic) );
    if (isquadratic == true)
      return true;
  }
  return false;
}

ArrayRef<double> ScipBackend::PrimalSolution() {
  if (HasFixedLP()) return FixedLP().primal;
  SCIP* scip = getSCIP();
  int num_vars = NumVars();
  auto sol = SCIPgetBestSol(scip);

  if (sol == nullptr)
    return std::vector<double>();

  std::vector<double> x(num_vars);
  SCIPgetSolVals(scip, sol, num_vars, getPROBDATA()->vars, x.data());
  return x;
}

pre::ValueMapDbl ScipBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> ScipBackend::DualSolution_LP() {
  if (IsMIP()) return HasFixedLP() ? ArrayRef<double>(FixedLP().dual) : ArrayRef<double>{};
  // SCIP prints the reason (presolving performed, ...) if asked
  if (!SCIPisDualSolAvailable(getSCIP(), verbose_mode()))
    return {};
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  for (int i=0; i < num_cons; i++)
    SCIP_CCALL( SCIPgetDualSolVal(getSCIP(), getPROBDATA()->linconss[i], pi.data() + i, NULL) );
  return pi;
}

double ScipBackend::ObjectiveValue() const {
  return SCIPgetPrimalbound(getSCIP());
}

int ScipBackend::NodeCount() const {
  return SCIPgetNNodes(getSCIP());
}

int ScipBackend::SimplexIterations() const {
  return SCIPgetNPrimalLPIterations(getSCIP()) + SCIPgetNDualLPIterations(getSCIP());
}

int ScipBackend::BarrierIterations() const {
  return SCIPgetNBarrierLPIterations(getSCIP());
}
std::map<std::string, std::variant<int, double, std::string>>
ScipBackend::SolutionStats() {
    std::map<std::string, std::variant<int, double, std::string>> stats;
    stats["simplex_iterations"] = SimplexIterations();
    stats["barrier_iterations"] = BarrierIterations();
    stats["node_count"] = NodeCount();
    return stats;
}

void ScipBackend::ExportModel(const std::string &file) {
  SCIP_CCALL( SCIPwriteOrigProblem(getSCIP(), file.data(), NULL, FALSE) );
}


void ScipBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptScip, getSCIP());
}

void ScipBackend::ApplyParamSettings() {
  auto* scip = getSCIP();
  if (storedOptions_.emphasis_ != 0)
    SCIP_CCALL( SCIPsetEmphasis(scip, (SCIP_PARAMEMPHASIS)storedOptions_.emphasis_, TRUE) );
  if (storedOptions_.heuristics_ != 0)
    SCIP_CCALL( SCIPsetHeuristics(scip, (SCIP_PARAMSETTING)storedOptions_.heuristics_, TRUE) );
  if (storedOptions_.cuts_ != 0)
    SCIP_CCALL( SCIPsetSeparating(scip, (SCIP_PARAMSETTING)storedOptions_.cuts_, TRUE) );
  if (storedOptions_.presolvings_ != 0)
    SCIP_CCALL( SCIPsetPresolving(scip, (SCIP_PARAMSETTING)storedOptions_.presolvings_, TRUE) );
  if (0 == storedOptions_.presolve_)
    SCIP_CCALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
}

bool ScipBackend::UseConcurrentSolve() const {
  if (storedOptions_.concurrent_ < 0)
    return storedOptions_.threads_ > 1;
  return storedOptions_.concurrent_ > 0;
}

void ScipBackend::Solve() {
  ApplyParamSettings();
  if (!storedOptions_.paramWrite_.empty())
    SCIP_CCALL( SCIPwriteParams(getSCIP(), storedOptions_.paramWrite_.c_str(), FALSE, FALSE) );
  if (UseConcurrentSolve())
    SCIP_CCALL( SCIPsolveConcurrent(getSCIP()) );
  else
    SCIP_CCALL( SCIPsolve(getSCIP()) );
  WindupSCIPSolve();
}

ScipBackend::FixedLPResult ScipBackend::SolveFixedLP(ArrayRef<double> x) {
  if (!SCIPgetBestSol(getSCIP())) MP_RAISE("no integer-feasible incumbent");
  int n = NumVars(), m = NumLinCons();
  auto* vars = getPROBDATA()->vars;
  std::unordered_map<SCIP_VAR*, int> columns;
  std::vector<double> obj(n), lb(n), ub(n), lhs(m), rhs(m);
  std::vector<std::vector<std::pair<int, double>>> entries(n);
  for (int j = 0; j < n; ++j) {
    columns[vars[j]] = j;
    obj[j] = SCIPvarGetObj(vars[j]);
    lb[j] = SCIPvarGetLbOriginal(vars[j]); ub[j] = SCIPvarGetUbOriginal(vars[j]);
    if (SCIPvarGetType(vars[j]) != SCIP_VARTYPE_CONTINUOUS)
      lb[j] = ub[j] = FixedIntegerValue(x[j], lb[j], ub[j], SCIPfeastol(getSCIP()));
  }
  for (int i = 0; i < m; ++i) {
    auto* con = getPROBDATA()->linconss[i];
    lhs[i] = SCIPgetLhsLinear(getSCIP(), con); rhs[i] = SCIPgetRhsLinear(getSCIP(), con);
    auto* rowvars = SCIPgetVarsLinear(getSCIP(), con);
    auto* values = SCIPgetValsLinear(getSCIP(), con);
    for (int k = 0; k < SCIPgetNVarsLinear(getSCIP(), con); ++k) {
      auto it = columns.find(rowvars[k]);
      if (it == columns.end()) MP_RAISE("fixed LP row contains an unmapped variable");
      entries[it->second].push_back({i, values[k]});
    }
  }
  std::vector<int> start(n), index;
  std::vector<double> value;
  for (int j = 0; j < n; ++j) {
    start[j] = index.size();
    for (auto entry : entries[j]) { index.push_back(entry.first); value.push_back(entry.second); }
  }
  SCIP_LPI* model = nullptr;
  SCIP_CCALL(SCIPlpiCreate(&model, nullptr, "mip:basis", static_cast<SCIP_OBJSEN>(SCIPgetObjsense(getSCIP()))));
  auto destroy = [](SCIP_LPI* model) { SCIPlpiFree(&model); };
  std::unique_ptr<SCIP_LPI, decltype(destroy)> fixed(model, destroy);
  SCIP_CCALL(SCIPlpiLoadColLP(model, static_cast<SCIP_OBJSEN>(SCIPgetObjsense(getSCIP())), n, obj.data(),
      lb.data(), ub.data(), nullptr, m, lhs.data(), rhs.data(), nullptr,
      index.size(), start.data(), index.data(), value.data()));
  if (fixed_model_time_limit() == 0) MP_RAISE("fixed LP time limit 0");
  if (fixed_model_time_limit() > 0)
    SCIP_CCALL(SCIPlpiSetRealpar(model, SCIP_LPPAR_LPTILIM, fixed_model_time_limit()));
  if (fixed_model_iter_limit() >= 0)
    SCIP_CCALL(SCIPlpiSetIntpar(model, SCIP_LPPAR_LPITLIM, fixed_model_iter_limit()));
  SCIP_CCALL(SCIPlpiSolveDual(model));
  if (!SCIPlpiIsOptimal(model)) MP_RAISE("fixed LP not optimal (limit or numerical failure)");
  FixedLPResult result;
  result.primal.resize(n);
  result.dual.resize(m);
  result.varstt.resize(n);
  result.constt.resize(m);
  SCIP_CCALL(SCIPlpiGetSol(model, &result.objective, result.primal.data(), result.dual.data(), nullptr, nullptr));
  result.objective = SCIPgetOrigObjscale(getSCIP()) * result.objective + SCIPgetOrigObjoffset(getSCIP());
  SCIP_CCALL(SCIPlpiGetBase(model, result.varstt.data(), result.constt.data()));
  auto convert = [](int status) {
    switch (status) {
    case SCIP_BASESTAT_BASIC: return (int)BasicStatus::bas;
    case SCIP_BASESTAT_LOWER: return (int)BasicStatus::low;
    case SCIP_BASESTAT_UPPER: return (int)BasicStatus::upp;
    case SCIP_BASESTAT_ZERO: return (int)BasicStatus::sup;
    default: MP_RAISE(fmt::format("Unknown SCIP fixed-LP basis status {}", (int)status));
    }
  };
  for (auto& status : result.varstt) status = convert(status);
  for (auto& status : result.constt) status = convert(status);
  SCIP_CCALL(SCIPlpiGetIterations(model, &result.iterations));
  return result;
}

void ScipBackend::WindupSCIPSolve() { }

void ScipBackend::ReportResults() {
  ReportSCIPResults();
  ConsiderEmulatedFixedModel();
  BaseBackend::ReportResults();
}

void ScipBackend::ReportSCIPResults() {
  SetStatus( GetSolveResult() );
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
bool ScipBackend::AcceptPoolSolution(double obj, double best) const {
  double gap = std::fabs(obj - best);
  if (storedOptions_.poolGapAbs_ >= 0.0 && gap > storedOptions_.poolGapAbs_)
    return false;
  if (storedOptions_.poolGap_ >= 0.0 && gap > storedOptions_.poolGap_ * std::fabs(best))
    return false;
  return true;
}

void ScipBackend::ReportSCIPPool() {
  if (!IsMIP())
    return;
  int nsolutions = SCIPgetNSols(getSCIP());        // sorted best-first
  if (storedOptions_.poolLimit_ >= 0)
    nsolutions = std::min(nsolutions, storedOptions_.poolLimit_);
  double best = nsolutions ? getPoolObjective(0) : 0.0;
  for (int i = 0; i < nsolutions; ++i) {
    double obj = getPoolObjective(i);
    if (!AcceptPoolSolution(obj, best))
      break;                                       // sorted: the rest is worse
    ReportIntermediateSolution({ getPoolSolution(i), {}, { obj } });
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

std::pair<int, std::string> ScipBackend::GetSolveResult() {
  namespace sol = mp::sol;
  SCIP_STATUS status = SCIPgetStatus(getSCIP());
  auto solu = SCIPgetBestSol(getSCIP());
  bool hasSol = (solu != nullptr);
  switch (status) {
  case SCIP_STATUS_UNKNOWN:
    return { sol::UNKNOWN, "solving status not yet known" };
  case SCIP_STATUS_USERINTERRUPT:
    if (hasSol)
      return { sol::LIMIT_FEAS_INTERRUPT, "interrupted, feasible solution" };
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted, no solution" };
  case SCIP_STATUS_NODELIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS_NODES, "node limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_NODES, "node limit, no solution" };
  case SCIP_STATUS_TOTALNODELIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS_NODES,
            "total node limit (incl. restarts), feasible solution" };
    return { sol::LIMIT_NO_FEAS_NODES,
          "total node limit (incl. restarts), no solution" };
  case SCIP_STATUS_STALLNODELIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS_NODES,
            "stalling node limit (no improvement w.r.t. primal bound), feasible solution" };
    return { sol::LIMIT_NO_FEAS_NODES,
          "stalling node limit (no improvement w.r.t. primal bound), no solution" };
  case SCIP_STATUS_TIMELIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS_TIME, "time limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_TIME, "time limit, no solution" };
  case SCIP_STATUS_MEMLIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS_SOFTMEM, "memory limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_SOFTMEM, "memory limit, no solution" };
  case SCIP_STATUS_GAPLIMIT:
    return { sol::LIMIT_FEAS_GAP, "gap limit" };
  case SCIP_STATUS_PRIMALLIMIT:
    return { sol::LIMIT_FEAS_BESTOBJ, "primal bound limit (lim:bestobjstop)" };
  case SCIP_STATUS_DUALLIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS_BESTBND, "dual bound limit (lim:bestbndstop), feasible solution" };
    return { sol::LIMIT_NO_FEAS_BESTBND, "dual bound limit (lim:bestbndstop), no solution" };
  case SCIP_STATUS_SOLLIMIT:
    return { sol::LIMIT_FEAS_NUMSOLS, "solution limit" };
  case SCIP_STATUS_BESTSOLLIMIT:
    return { sol::LIMIT_FEAS_NUMSOLS, "solution improvement limit" };
  case SCIP_STATUS_RESTARTLIMIT:
    if (hasSol)
      return { sol::LIMIT_FEAS, "restart limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS, "restart limit, no solution" };
  case SCIP_STATUS_OPTIMAL:
    return { sol::SOLVED, "optimal solution" };
  case SCIP_STATUS_INFEASIBLE:
    return { sol::INFEASIBLE, "infeasible problem" };
  case SCIP_STATUS_UNBOUNDED:
    if (hasSol)
      return { sol::UNBOUNDED_FEAS, "unbounded problem, feasible solution" };
    return { sol::UNBOUNDED_NO_FEAS, "unbounded problem, no solution" };
  case SCIP_STATUS_INFORUNBD:
    return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
  case SCIP_STATUS_TERMINATE:
    if (hasSol)
      return { sol::LIMIT_FEAS_INTERRUPT,
            "process received a SIGTERM signal, feasible solution" };
    return { sol::LIMIT_NO_FEAS_INTERRUPT,
          "process received a SIGTERM signal, no solution" };
  default:
    if (hasSol)
      return { sol::UNCERTAIN, "unknown, solution candidate returned" };
    return { sol::UNKNOWN, "unknown" };
  }
  return { sol::UNKNOWN, "not solved" };
}


void ScipBackend::FinishOptionParsing() {
  if (!storedOptions_.paramRead_.empty())
    SCIP_CCALL( SCIPreadParams(getSCIP(), storedOptions_.paramRead_.c_str()) );
  ApplyInlineParams();
  ApplyThreads();
  SetUpLogging();
}

void ScipBackend::ApplyInlineParams() {
  for (const auto& text: storedOptions_.inlineParams_) {
    std::string pairs = text;
    std::replace(pairs.begin(), pairs.end(), '=', ' ');
    std::istringstream in(pairs);
    std::string name, value;
    while (in >> name) {
      if (!(in >> value))
        MP_RAISE(fmt::format(
            "tech:optionnative: no value for SCIP parameter '{}'", name));
      SetNativeParam(name, value);
    }
  }
}

void ScipBackend::ApplyThreads() {
  if (storedOptions_.threads_ <= 0)
    return;
  SetSolverOption("parallel/maxnthreads", storedOptions_.threads_);
  SetSolverOption("lp/threads", storedOptions_.threads_);
}

void ScipBackend::SetUpLogging() {
  if (storedOptions_.outlev_ == 1)
    SetSolverOption("display/verblevel", 4);
  int verblevel = 0;
  GetSolverOption("display/verblevel", verblevel);
  bool console_log = verblevel > 0;
  if (!storedOptions_.logFile_.empty()) {
    SCIPsetMessagehdlrLogfile(getSCIP(), storedOptions_.logFile_.c_str());
    if (!console_log)             // log file only: full log, quiet console
      SetSolverOption("display/verblevel", 4);
  }
  // Quiet also hides messages ignoring verblevel (concurrent solve);
  // errors go to stderr regardless.
  SCIPsetMessagehdlrQuiet(getSCIP(), !console_log);
  set_verbose_mode(console_log);
}


////////////////////////////// OPTIONS /////////////////////////////////
static const mp::OptionValueInfo lp_values_method[] = {
  { "s", "automatic simplex (default)", -1},
  { "p", "primal simplex", 0},
  { "d", "dual simplex", 1},
  { "b", "barrier (SoPlex has none: runs the dual simplex)", 2},
  { "c", "barrier with crossover (SoPlex: dual simplex)", 3},
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

static const mp::OptionValueInfo estimation_method[] = {
  { "c", "completion", 0},
  { "e", "ensemble", 1},
  { "g", "time series forecasts on either gap", 2},
  { "l", "leaf frequency", 3},
  { "o", "open nodes", 4},
  { "w", "tree weight (default)", 5},
  { "s", "ssg", 6},
  { "t", "tree profile", 7},
  { "b", "weighted backtrack estimate (wbe)", 8}
};

static const mp::OptionValueInfo estimation_completion[] = {
  { "a", "auto (default)", 0},
  { "g", "gap", 1},
  { "w", "tree weight", 2},
  { "m", "monotone regression", 3},
  { "r", "regression forest", 4},
  { "s", "ssg", 5}
};

static const mp::OptionValueInfo values_lpscaling[] = {
  { "0", "none", 0},
  { "1", "normal (default)", 1},
  { "2", "aggressive", 2}
};

static const mp::OptionValueInfo values_symmetry[] = {
  { "0", "off", 0},
  { "1", "constraint-based (orbitopes, symresacks; lexicographic and orbitopal reduction if dynamic)", 1},
  { "2", "orbital reduction", 2},
  { "3", "orbitopes and symresacks, and lexicographic/orbital reduction", 3},
  { "4", "Schreier Sims cuts", 4},
  { "5", "Schreier Sims cuts, orbitopes, symresacks, and/or lexicographic reduction", 5},
  { "6", "Schreier Sims cuts, orbital reduction", 6},
  { "7", "Schreier Sims cuts, orbitopes, symresacks, and/or lexicographic/orbital reduction (default)", 7}
};

static const mp::OptionValueInfo values_clocktype[] = {
  { "1", "CPU user seconds", 1},
  { "2", "wall clock time (default)", 2}
};

static const mp::OptionValueInfo values_emphasis[] = {
  { "0", "default settings (default)", 0},
  { "1", "cpsolver: CP-like search (e.g., no LP relaxation)", 1},
  { "2", "easycip: solve easy problems fast", 2},
  { "3", "feasibility: detect feasibility fast", 3},
  { "4", "hardlp: be capable to handle hard LPs", 4},
  { "5", "optimality: prove optimality fast", 5},
  { "6", "counter: feasible and fast counting process", 6},
  { "7", "phasefeas: feasibility phase of the 3-phase approach", 7},
  { "8", "phaseimprove: improvement phase of the 3-phase approach", 8},
  { "9", "phaseproof: proof phase of the 3-phase approach", 9},
  { "10", "numerics: increased numerical safety", 10},
  { "11", "benchmark: do not try to avoid running into the memory limit", 11}
};

/// Cut families exposed as cut:<family> = separation frequency
struct CutFamily {
  const char* names;
  const char* key;
  const char* descr;
};
static const CutFamily cut_families[] = {
  { "cut:gomory gomorycuts", "separating/gomory/freq",
    "Frequency of Gomory mixed-integer cuts (default 10)" },
  { "cut:gomorymi gomorymicuts", "separating/gomorymi/freq",
    "Frequency of Gomory mixed-integer cuts from the LP tableau (default 10)" },
  { "cut:strongcg strongcgcuts", "separating/strongcg/freq",
    "Frequency of strong Chvatal-Gomory cuts (default 10)" },
  { "cut:mir mircuts cut:cmir", "separating/cmir/freq",
    "Frequency of complemented mixed-integer rounding (c-MIR) cuts (default 10)" },
  { "cut:flowcover flowcovercuts", "separating/flowcover/freq",
    "Frequency of flow cover cuts (default 10)" },
  { "cut:cover covercuts cut:knapsackcover", "separating/knapsackcover/freq",
    "Frequency of knapsack cover cuts (default 10)" },
  { "cut:clique cliquecuts", "separating/clique/freq",
    "Frequency of clique cuts (default 0: root node only)" },
  { "cut:zerohalf zerohalfcuts", "separating/zerohalf/freq",
    "Frequency of {0,1/2}-cuts (default 10)" },
  { "cut:implied impliedcuts cut:impliedbounds", "separating/impliedbounds/freq",
    "Frequency of implied bound cuts (default 10)" },
  { "cut:aggregation aggregationcuts", "separating/aggregation/freq",
    "Frequency of aggregated (c-MIR and flow cover) cuts (default 10)" },
  { "cut:oddcycle oddcyclecuts", "separating/oddcycle/freq",
    "Frequency of odd cycle cuts (default -1: off)" },
  { "cut:mcf mcfcuts", "separating/mcf/freq",
    "Frequency of multi-commodity-flow network cuts (default 0: root node only)" },
  { "cut:disjunctive disjunctivecuts", "separating/disjunctive/freq",
    "Frequency of disjunctive cuts (default 0: root node only)" },
  { "cut:rlt rltcuts", "separating/rlt/freq",
    "Frequency of reformulation-linearization-technique cuts (default 0: root node only)" }
};

static const mp::OptionValueInfo childsel[] = {
  { "d", "down", 0},
  { "u", "up", 1},
  { "p", "pseudo costs", 2},
  { "i", "inference", 3},
  { "l", "lp value", 4},
  { "r", "root LP value difference", 5},
  { "h", "hybrid inference/root LP value difference (default)", 6}
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

  AddStoredOption("tech:outlev outlev",
    "0*/1: Whether to write SCIP log lines (chatter) to stdout and to file.",
    storedOptions_.outlev_);

  AddScipOption("tech:outlevnative outlevnative tech:outlev-native outlev-native",
    "0*/1/2/3/4/5: Native SCIP verbosity level (display/verblevel; the SCIP "
    "default is 4, the driver's is 0). tech:outlev=1 selects 4.",
    "display/verblevel", 0, 5);

  AddStoredOption("tech:logfile logfile",
    "Log file name.",
    storedOptions_.logFile_);

  AddStoredOption("tech:optionnativeread optionnativeread tech:param:read param:read",
    "Filename of SCIP parameter file (as path)."
    "The suffix on a parameter file should be .set.\n",
    storedOptions_.paramRead_);

  AddStoredOption("tech:optionnativewrite optionnativewrite tech:param:write param:write",
    "Name of a SCIP parameter file (.set) to write before solving, "
    "with all parameters.",
    storedOptions_.paramWrite_);

  AddListOption("tech:optionnative optionnative optnative tech:param param",
    "General way to set any SCIP parameter; the value should be a quoted "
    "string (delimited by ' or \") containing one or more "
    "'name=value' or 'name value' pairs separated by spaces, for example "
    "'separating/maxrounds=5 heuristics/rins/freq -1'. Can appear more "
    "than once. Cannot be used to query current parameter values.",
    storedOptions_.inlineParams_);

  AddStoredOption("tech:threads threads",
    "Number of threads to use (default 0: SCIP defaults). A value > 1 "
    "runs SCIP's concurrent solve with at most this many threads "
    "(unless alg:concurrent=0) and sets lp:threads for a multi-threaded "
    "LP solver; the bundled SoPlex is single-threaded. "
    "SCIP's branch-and-bound itself is sequential.",
    storedOptions_.threads_, 0, 64);

  AddScipOption("tech:clocktype clocktype",
    "Clock type for time limits and reporting:\n"
    "\n.. value-table::\n",
    "timing/clocktype", values_clocktype, 2);

  AddScipOption("alg:method method lpmethod",
    "LP algorithm for solving initial LP relaxations:\n"
    "\n.. value-table::\n", "lp/initalgorithm", lp_values_method, "s");

  AddScipOption("alg:remethod remethod relpmethod",
    "LP algorithm for resolving LP relaxations if a starting basis exists:\n"
    "\n.. value-table::\n", "lp/resolvealgorithm", lp_values_method, "s");

  AddStoredOption("alg:concurrent concurrent",
    "-1/0/1: whether to solve the problem using concurrent solvers "
    "(several SCIP instances with different settings, "
    "options par:*, tech:threads)"
    "\n"
    "  | -1 - Concurrent solvers are used if tech:threads > 1 (default)\n"
    "  |  0 - No concurrent solvers\n"
    "  |  1 - Concurrent solvers are used, par:maxnthreads at most.",
    storedOptions_.concurrent_, -1, 1);

  AddStoredOption("alg:emphasis emphasis",
    "Predefined parameter settings (SCIPsetEmphasis), applied before "
    "solving, after individual parameter settings:\n"
    "\n.. value-table::\n",
    storedOptions_.emphasis_, values_emphasis);

  /////////////////////// BRANCHING ///////////////////////
  AddScipOption("branch:preferbinary preferbinary",
    "0/1: whether branching on binary variables should be preferred"
    "\n"
    "  | 0 - Binary variables should not be preferred on branching (default)\n"
    "  | 1 - Binary variables should be preferred on branching.",
    "branching/preferbinary", 0, 1);

  //////////////////////// CUTS //////////////////////////
  AddScipOption("cut:dircutoffdistweight dircutoffdistweight",
    "Weight of directed cutoff distance in cut score calculation (default: 0.0)",
    "cutselection/hybrid/dircutoffdistweight", 0.0, SCIP_REAL_MAX);

  AddScipOption("cut:efficacyweight efficacyweight",
    "Weight of efficacy in cut score calculation (default: 1.0)",
    "cutselection/hybrid/efficacyweight", 0.0, SCIP_REAL_MAX);

  AddScipOption("cut:intsupportweight intsupportweight",
    "Weight of integral support in cut score calculation (default: 0.1)",
    "cutselection/hybrid/intsupportweight", 0.0, SCIP_REAL_MAX);

  AddScipOption("cut:minortho minortho",
    "Minimal orthogonality for a cut to enter the LP (0..1; default: 0.9)",
    "cutselection/hybrid/minortho", 0.0, 1.0);

  AddScipOption("cut:minorthoroot minorthoroot",
    "Minimal orthogonality for a cut to enter the LP in the root node (0..1; default: 0.9)",
    "cutselection/hybrid/minorthoroot", 0.0, 1.0);

  AddScipOption("cut:objparalweight objparalweight",
    "Weight of objective parallelism in cut score calculation (default: 0.1)",
    "cutselection/hybrid/objparalweight", 0.0, SCIP_REAL_MAX);


  AddScipOption("cut:maxcuts maxcuts",
    "Maximal number of cuts separated per separation round (0: disable local separation; default: 100)",
    "separating/maxcuts", 0, INT_MAX);

  AddScipOption("cut:maxcutsroot maxcutsroot",
    "Maximal number of separated cuts per separation round at the root node (0: disable root node separation; default: 2000)",
    "separating/maxcutsroot", 0, INT_MAX);

  AddScipOption("cut:maxrounds",
    "Maximal number of separation rounds per node (default: -1: unlimited)",
    "separating/maxrounds", -1, INT_MAX);

  AddScipOption("cut:maxroundsroot",
    "Maximal number of separation rounds in the root node (default: -1: unlimited)",
    "separating/maxroundsroot", -1, INT_MAX);

  AddScipOption("cut:maxstallrounds",
    "Maximal number of consecutive separation rounds without objective or integrality improvement in local nodes (-1: no additional restriction; default: 1)",
    "separating/maxstallrounds", -1, INT_MAX);

  AddScipOption("cut:maxstallroundsroot",
    "Maximal number of consecutive separation rounds without objective or integrality improvement in the root node (-1: no additional restriction; default: 10)",
    "separating/maxstallroundsroot", -1, INT_MAX);

  AddScipOption("cut:minactivityquot",
    "Minimum cut activity quotient to convert cuts into constraints during a restart (0.0: all cuts are converted; default: 0.8)",
    "separating/minactivityquot", 0.0, 1.0);

  AddScipOption("cut:minefficacy",
    "Minimal efficacy for a cut to enter the LP (default: 0.0001)",
    "separating/minefficacy", 0.0, SCIP_REAL_MAX);

  AddScipOption("cut:minefficacyroot",
    "Minimal efficacy for a cut to enter the LP in the root node (default: 0.0001)",
    "separating/minefficacyroot", 0.0, SCIP_REAL_MAX);

  AddScipOption("cut:poolfreq",
    "Separation frequency for the global cut pool (-1: disable global cut pool; 0: only separate pool at the root; default: 10)",
    "separating/poolfreq", -1, SCIP_MAXTREEDEPTH);

  for (const auto& family: cut_families)
    AddScipOption(family.names,
      fmt::format("{}: separation frequency (-1: off, 0: root node only, "
                  "n: every n-th depth level)", family.descr).c_str(),
      family.key, -1, SCIP_MAXTREEDEPTH);

  AddStoredOption("cut:settings",
    "0/1/2/3: sets cuts settings"
    "\n"
    "  | 0 - Sets cuts default (default)\n"
    "  | 1 - Sets cuts aggressive\n"
    "  | 2 - Sets cuts fast\n"
    "  | 3 - Sets cuts off.",
    storedOptions_.cuts_, 0, 3);

  ////////////////////// HEURISTICS //////////////////////
  AddStoredOption("heu:settings",
    "0/1/2/3: sets heuristics settings"
    "\n"
    "  | 0 - Sets heuristics default (default)\n"
    "  | 1 - Sets heuristics aggressive\n"
    "  | 2 - Sets heuristics fast\n"
    "  | 3 - Sets heuristics off.",
    storedOptions_.heuristics_, 0, 3);

  //////////////////////// LIMITS ////////////////////////
  AddScipOption("mip:gapabs mipgapabs lim:absgap absgap",
    "Solving stops, if the absolute gap = |primalbound - dualbound| is below the given value (default: 0.0)",
    "limits/absgap", 0.0, SCIP_REAL_MAX);

  AddScipOption("lim:autorestartnodes",
    "If solve exceeds this number of nodes for the first time, an automatic restart is triggered (default: -1: no automatic restart)",
    "limits/autorestartnodes", -1, INT_MAX);

  AddScipOption("lim:bestbndstop bestbndstop lim:dual",
    "Solving stops, if the dual bound is at least as good as this value (default: no limit)",
    "limits/dual", -SCIP_REAL_MAX, SCIP_REAL_MAX);

  AddScipOption("lim:bestobjstop bestobjstop lim:primal",
    "Solving stops, if the primal bound is at least as good as this value (default: no limit)",
    "limits/primal", -SCIP_REAL_MAX, SCIP_REAL_MAX);

  AddScipOption("lim:iisnodes iisnodes",
    "Maximal number of branch-and-bound nodes for the IIS finders (default: -1: no limit)",
    "iis/nodes", -1, INT_MAX);

  AddScipOption("lim:iistime iistime",
    "Maximal time in seconds for the IIS finders (default: 1e+20)",
    "iis/time", 0.0, 1e+20);

  AddScipOption("lim:iter iterlim lp:iterlim",
    "Iteration limit for each single LP solve (default: -1: no limit). "
    "Unlike in other solvers, this is not a total limit.",
    "lp/iterlim", -1, INT_MAX);

  AddScipOption("lim:rootiter rootiterlim lp:rootiterlim",
    "Iteration limit for the initial root LP solve (default: -1: no limit)",
    "lp/rootiterlim", -1, INT_MAX);

  AddScipOption("lim:bestsol",
    "Solving stops, if the given number of solution improvements were found (default: -1: no limit)",
    "limits/bestsol", -1, INT_MAX);

  AddScipOption("mip:gap mipgap lim:gap gap",
    "Solving stops, if the relative gap = |primal - dual|/MIN(|dual|,|primal|) is below the given value, the gap is 'Infinity', if primal and dual bound have opposite signs (default: 0.0)",
    "limits/gap", 0.0, SCIP_REAL_MAX);

  AddScipOption("lim:maxorigsol",
    "Maximal number of solutions candidates to store in the solution storage of the original problem (default: 10)",
    "limits/maxorigsol", 0, INT_MAX);

  AddScipOption("lim:maxsol",
    "Maximal number of solutions to store in the solution storage (default: 100)",
    "limits/maxsol", 1, INT_MAX);

  AddScipOption("lim:memory memory memlim",
    "#maximal memory usage in MB; reported memory usage is lower than real memory usage! (default: 8796093022207.0)",
    "limits/memory", 0.0, (SCIP_Real)SCIP_MEM_NOLIMIT);

  AddScipOption("lim:nodes nodelim nodelimit",
    "Maximal number of nodes to process (default: -1: no limit)",
    "limits/nodes", -1, INT_MAX);

  AddScipOption("lim:restarts",
    "Solving stops, if the given number of restarts was triggered (default: -1: no limit)",
    "limits/restarts", -1, INT_MAX);

  AddScipOption("lim:softtime softtime",
    "Soft time limit which should be applied after first solution was found (default: -1.0: disabled)",
    "limits/softtime", -1.0, SCIP_REAL_MAX);

  AddScipOption("lim:sol sollim solutionlimit lim:solutions",
    "Solving stops, if the given number of solutions were found; this limit is first checked in presolving (default: -1: no limit)",
    "limits/solutions", -1, INT_MAX);

  AddScipOption("lim:stallnodes",
    "Solving stops, if the given number of nodes was processed since the last improvement of the primal solution value (default: -1: no limit)",
    "limits/stallnodes", -1, INT_MAX);

  AddScipOption("lim:time timelim timelimit time_limit",
    "Limit on solve time (in seconds; default: 1e+20).",
    "limits/time", 0.0, 1e+20);

  AddScipOption("lim:totalnodes",
    "Maximal number of total nodes (incl. restarts) to process (default: -1: no limit)",
    "limits/totalnodes", -1, INT_MAX);

  ///////////////////////// IIS //////////////////////////
  AddScipOption("iis:bounds iisbounds",
    "0/1: whether variable bounds take part in the IIS analysis "
    "(SCIP parameter iis/removebounds; driver default 1: a variable in "
    "the IIS gets .iis = low/upp/fix, with 0 all stay non)",
    "iis/removebounds", 0, 1);

  AddScipOption("iis:irreducible iisirreducible",
    "0/1: whether the infeasible subsystem should be irreducible "
    "(an IIS, default 1) rather than any infeasible subsystem",
    "iis/irreducible", 0, 1);

  ////////////////////////// LP //////////////////////////
  AddScipOption("lp:pricing pricing",
    "Pricing strategy:\n"
    "\n.. value-table::",
    "lp/pricing", values_pricing, "l");

  AddScipOption("lp:presolving",
    "0/1: whether presolving of LP solver should be used"
    "\n"
    "  | 0 - Presolving of LP solver should not be used\n"
    "  | 1 - Presolving of LP solver should be used (default).",
    "lp/presolving", 0, 1);

  AddScipOption("lp:scaling pre:scale scale",
    "LP scaling:\n"
    "\n.. value-table::\n",
    "lp/scaling", values_lpscaling, 1);

  AddScipOption("lp:threads",
    "Number of threads used for solving the LP (default: 0: automatic)",
    "lp/threads", 0, 64);

  AddScipOption("lp:alwaysgetduals alwaysgetfarkasduals alwaysgetduals",
    "0/1: whether the Farkas duals should always be collected when an LP is found to be infeasible"
    "\n"
    "  | 0 - The Farkas duals should not always be collected when an LP is found to be infeasible (default)\n"
    "  | 1 - The Farkas duals should always be collected when an LP is found to be infeasible.",
    "lp/alwaysgetduals", 0, 1);

  AddScipOption("lp:solvedepth",
    "Maximal depth for solving LP at the nodes (default: -1: no depth limit)",
    "lp/solvedepth", -1, SCIP_MAXTREEDEPTH);

  AddScipOption("lp:solvefreq",
    "Frequency for solving LP at the nodes (-1: never; 0: only root LP; default: 1)",
    "lp/solvefreq", -1, SCIP_MAXTREEDEPTH);

  ///////////////////////// MISC /////////////////////////
  AddScipOption("misc:allowstrongdualreds allowstrongdualreds",
    "0/1: whether strong dual reductions should be allowed in propagation and presolving"
    "\n"
    "  | 0 - Strong dual reductions should not be allowed in propagation and presolving\n"
    "  | 1 - Strong dual reductions should be allowed in propagation and presolving (default).",
    "misc/allowstrongdualreds", 0, 1);

  AddScipOption("misc:allowweakdualreds allowweakdualreds",
    "0/1: whether weak dual reductions should be allowed in propagation and presolving"
    "\n"
    "  | 0 - Weak dual reductions should not be allowed in propagation and presolving\n"
    "  | 1 - Weak dual reductions should be allowed in propagation and presolving (default).",
    "misc/allowweakdualreds", 0, 1);

  AddScipOption("misc:scaleobj scaleobj",
    "0/1: whether the objective function should be scaled so that it is always integer"
    "\n"
    "  | 0 - The objective function should not be scaled so that it is always integer\n"
    "  | 1 - The objective function should be scaled so that it is always integer (default).",
    "misc/scaleobj", 0, 1);

  AddScipOption("mip:symmetry symmetry misc:usesymmetry",
    "Bitset describing the symmetry handling techniques:\n"
    "\n.. value-table::\n",
    "misc/usesymmetry", values_symmetry, 7);

  ////////////////////////// NLP /////////////////////////
  AddScipOption("nlp:disable",
    "0/1: whether the NLP relaxation should be always disabled (also for NLPs/MINLPs)"
    "\n"
    "  | 0 - NLP relaxation should not be always disabled (default)\n"
    "  | 1 - NLP relaxation should be always disabled.",
    "nlp/disable", 0, 1);

  ////////////////////// NUMERICS ////////////////////////
  AddScipOption("num:checkfeastolfac checkfeastolfac",
    "Feasibility tolerance factor; for checking the feasibility of the best solution (default: 1.0)",
    "numerics/checkfeastolfac", 0.0, SCIP_REAL_MAX);

  AddScipOption("alg:dualfeastol dualfeastol num:dualfeastol",
    "Feasibility tolerance for reduced costs in LP solution (default: 1e-07)",
    "numerics/dualfeastol", SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON);

  AddScipOption("num:epsilon epsilon",
    "Absolute values smaller than this are considered zero (default: 1e-09)",
    "numerics/epsilon", SCIP_MINEPSILON, SCIP_MAXEPSILON);

  AddScipOption("alg:feastol feastol num:feastol",
    "Feasibility tolerance for constraints (default: 1e-06)",
    "numerics/feastol", SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON);

  AddScipOption("num:infinity infinity",
    "Values larger than this are considered infinity (default: 1e+20)",
    "numerics/infinity", 1e+10, SCIP_INVALID/10.0);

  AddScipOption("num:lpfeastolfactor lpfeastolfactor",
    "Factor w.r.t. primal feasibility tolerance that determines default (and maximal) primal feasibility tolerance of LP solver (default: 1.0)",
    "numerics/lpfeastolfactor", 1e-6, 1.0);

  AddScipOption("num:sumepsilon sumepsilon",
    "Absolute values of sums smaller than this are considered zero (default: 1e-06)",
    "numerics/sumepsilon", SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON);

  //////////////////// NODESELECTION /////////////////////
  AddScipOption("nod:childsel",
    "Child selection rule:\n"
    "\n.. value-table::\n", "nodeselection/childsel", childsel, "h");

  ////////////////////// PARALLEL ////////////////////////
  AddScipOption("par:maxnthreads maxnthreads",
    "Maximum number of threads used during parallel solve (default: 8)",
    "parallel/maxnthreads", 0, 64);

  AddScipOption("par:minnthreads minnthreads",
    "Minimum number of threads used during parallel solve (default: 1)",
    "parallel/minnthreads", 0, 64);

  AddScipOption("par:mode mode",
    "0/1: Parallel optimisation mode"
    "\n"
    "  | 0 - Opportunistic\n"
    "  | 1 - Deterministic (default)",
    "parallel/mode", 0, 1);

  ////////////////////// PRESOLVE ////////////////////////
  AddScipOption("pre:abortfac abortfac",
    "Abort presolve, if at most this fraction of the problem was changed in last presolve round (default: 0.0008)",
    "presolving/abortfac", 0.0, 1.0);

  AddScipOption("pre:clqtablefac clqtablefac",
    "Limit on number of entries in clique table relative to number of problem nonzeros (default: 2.0)",
    "presolving/clqtablefac", 0.0, SCIP_REAL_MAX);

  AddScipOption("pre:donotaggr donotaggr",
    "0/1: whether aggregation of variables should be forbidden"
    "\n"
    "  | 0 - Aggregation of variables should not be forbidden (default)\n"
    "  | 1 - Aggregation of variables should be forbidden.",
    "presolving/donotaggr", 0, 1);

  AddScipOption("pre:donotmultaggr donotmultaggr",
    "0/1: whether multi-aggregation of variables should be forbidden"
    "\n"
    "  | 0 - Multi-aggregation of variables should not be forbidden (default)\n"
    "  | 1 - Multi-aggregation of variables should be forbidden.",
    "presolving/donotmultaggr", 0, 1);

  AddScipOption("pre:immrestartfac immrestartfac",
    "Fraction of integer variables that were fixed in the root node triggering an immediate restart with preprocessing (default: 0.05)",
    "presolving/immrestartfac", 0.0, 1.0);

  AddScipOption("pre:restartfac restartfac",
    "Fraction of integer variables that were fixed in the root node triggering a restart with preprocessing after root node evaluation (default: 0.025)",
    "presolving/restartfac", 0.0, 1.0);

  AddScipOption("pre:restartminred restartminred",
    "Minimal fraction of integer variables removed after restart to allow for an additional restart (default: 0.05)",
    "presolving/restartminred", 0.0, 1.0);

  AddScipOption("pre:subrestartfac subrestartfac",
    "Fraction of integer variables that were globally fixed during the solving process triggering a restart with preprocessing (default: 1.0)",
    "presolving/subrestartfac", 0.0, 1.0);

  AddScipOption("pre:maxrestarts",
    "Maximal number of restarts (default: -1: unlimited)",
    "presolving/maxrestarts", -1, INT_MAX);

  AddScipOption("pre:maxrounds",
    "Maximal number of presolving rounds (default: -1: unlimited; 0: off)",
    "presolving/maxrounds", -1, INT_MAX);

  AddStoredOption("pre:settings",
    "0/1/2/3: sets presolvings settings"
    "\n"
    "  | 0 - Sets presolvings default (default)\n"
    "  | 1 - Sets presolvings aggressive\n"
    "  | 2 - Sets presolvings fast\n"
    "  | 3 - Sets presolvings off.",
    storedOptions_.presolvings_, 0, 3);

  AddStoredOption("pre:solve presolve",
    "0/1: whether to use SCIP presolving (default 1); "
    "0 is the same as pre:settings=3.",
    storedOptions_.presolve_, 0, 1);

  ////////////////////// PROPAGATE ///////////////////////
  AddScipOption("pro:abortoncutoff",
    "0/1: whether propagation should be aborted immediately (setting this to 0 could help conflict analysis to produce more conflict constraints)"
    "\n"
    "  | 0 - Propagation should not be aborted immediately\n"
    "  | 1 - Propagation should be aborted immediately (default).",
    "propagating/abortoncutoff", 0, 1);

  AddScipOption("pro:maxrounds",
    "Maximal number of propagation rounds per node (-1: unlimited; 0: off; default: 100)",
    "propagating/maxrounds", -1, INT_MAX);

  AddScipOption("pro:maxroundsroot",
    "Maximal number of propagation rounds in the root node (-1: unlimited; 0: off; default: 1000)",
    "propagating/maxroundsroot", -1, INT_MAX);
  
  //////////////////// RANDOMIZATION /////////////////////
  AddScipOption("ran:permuteconss permuteconss",
    "0/1: whether the order of constraints should be permuted (depends on permutationseed)? "
    "\n"
    "  | 0 - Order of constraints should not be permuted\n"
    "  | 1 - Order of constraints should be permuted (default).",
    "randomization/permuteconss", 0, 1);

  AddScipOption("ran:permutevars permutevars",
    "0/1: whether the order of variables should be permuted (depends on permutationseed)? "
    "\n"
    "  | 0 - Order of variables should not be permuted (default)\n"
    "  | 1 - Order of variables should be permuted.",
    "randomization/permutevars", 0, 1);

  AddScipOption("tech:lpseed lpseed ran:lpseed",
    "Random seed for LP solver, e.g. for perturbations in the simplex (default: 0: LP default)",
    "randomization/lpseed", 0, INT_MAX);

  AddScipOption("tech:permseed permseed ran:permutationseed permutationseed",
    "Seed value for permuting the problem after reading/transformation (default: 0: no permutation) ",
    "randomization/permutationseed", 0, INT_MAX);

  AddScipOption("tech:seed seed ran:randomseedshift randomseedshift",
    "Global shift of all random seeds in the plugins and the LP random seed (default: 0) ",
    "randomization/randomseedshift", 0, INT_MAX);

	///////////////////////// SOLUTION POOL ///////////////////////////////
  AddStoredOption("sol:poollimit poollimit ams_limit solnlimit",
    "Limit on the number of alternative solutions reported "
    "(options sol:stub, sol:count); default -1: all solutions stored by "
    "SCIP, at most lim:maxsol.",
    storedOptions_.poolLimit_, -1, INT_MAX);

  AddStoredOption("sol:poolgap poolgap ams_eps",
    "Relative tolerance for reporting alternative solutions "
    "(default: no limit).",
    storedOptions_.poolGap_, -1.0, SCIP_REAL_MAX);

  AddStoredOption("sol:poolgapabs poolgapabs ams_epsabs",
    "Absolute tolerance for reporting alternative solutions "
    "(default: no limit).",
    storedOptions_.poolGapAbs_, -1.0, SCIP_REAL_MAX);

  ///////////////////////// TREE /////////////////////////
  AddScipOption("est:method",
    "Tree size estimation method:\n"
    "\n.. value-table::\n", "estimation/method", estimation_method, "w");

  AddScipOption("est:completiontype",
    "Approximation of search tree completion:\n"
    "\n.. value-table::\n",
                  "estimation/completiontype", estimation_completion, "a");
}


double ScipBackend::MIPGap() {
  return
      SCIPgetGap(getSCIP())<Infinity() ?
             SCIPgetGap(getSCIP()) : AMPLInf();
}
double ScipBackend::BestDualBound() {
  return SCIPgetDualbound(getSCIP());
}

double ScipBackend::MIPGapAbs() {
  double gapabs =
      std::fabs(ObjectiveValue() - BestDualBound());
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
  if (IsMIP()) return FixedLPBasis(CG_Linear);
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

void ScipBackend::AddMIPStart(ArrayRef<double> x0, ArrayRef<int> sparsity) {
  auto mv = GetValuePresolver().PresolveSolution({x0});
  auto ms = GetValuePresolver().PresolveGenericInt({sparsity});
  auto x = mv.GetVarValues()();
  auto present = ms.GetVarValues()();
  SCIP_SOL* solution = nullptr;
  bool complete = present.size() == x.size() &&
      std::all_of(present.begin(), present.end(), [](int value) { return value != 0; });
  if (complete)
    SCIP_CCALL(SCIPcreateSol(getSCIP(), &solution, nullptr));
  else
    SCIP_CCALL(SCIPcreatePartialSol(getSCIP(), &solution, nullptr));
  for (int j = 0; j < NumVars(); ++j)
    if (j < (int)present.size() && present[j])
      SCIP_CCALL(SCIPsetSolVal(getSCIP(), solution, getPROBDATA()->vars[j], x[j]));
  SCIP_Bool stored;
  SCIP_CCALL(SCIPaddSolFree(getSCIP(), &solution, &stored));
}

/// The IIS finder leaves a copy of the infeasible subsystem,
/// with the surviving constraints and bounds under their original names.
void ScipBackend::ComputeIIS() {
  iisVars_.clear();
  iisCons_.clear();
  if (SCIPparamIsDefault(Param("iis/removebounds")))   // driver default: analyze bounds
    SetSolverOption("iis/removebounds", 1);
  SCIP_CCALL( SCIPgenerateIIS(getSCIP()) );
  SCIP_IIS* iis = SCIPgetIIS(getSCIP());
  SCIP* sub = iis ? SCIPiisGetSubscip(iis) : nullptr;
  if (!sub || !SCIPiisIsSubscipInfeasible(iis)) {
    AddToSolverMessage("SCIP found no infeasible subsystem (see lim:iistime, lim:iisnodes).\n");
    return;
  }
  std::unordered_map<std::string, int> conIndex, varIndex;
  for (int i = 0; i < NumLinCons(); ++i)
    conIndex[SCIPconsGetName(getPROBDATA()->linconss[i])] = i;
  for (int j = 0; j < NumVars(); ++j)
    varIndex[SCIPvarGetName(getPROBDATA()->vars[j])] = j;
  iisVars_.assign(NumVars(), (int)IISStatus::non);
  iisCons_.assign(NumLinCons(), (int)IISStatus::non);
  int ncons = 0, nbounds = 0;
  for (int i = 0; i < SCIPgetNOrigConss(sub); ++i) {
    auto it = conIndex.find(SCIPconsGetName(SCIPgetOrigConss(sub)[i]));
    if (it != conIndex.end()) {
      iisCons_[it->second] = (int)IISStatus::mem;
      ++ncons;
    }
  }
  for (int k = 0; k < SCIPgetNOrigVars(sub); ++k) {
    auto* var = SCIPgetOrigVars(sub)[k];
    auto it = varIndex.find(SCIPvarGetName(var));
    if (it == varIndex.end())
      continue;
    bool lb = !SCIPisInfinity(sub, -SCIPvarGetLbOriginal(var));
    bool ub = !SCIPisInfinity(sub, SCIPvarGetUbOriginal(var));
    if (lb || ub)
      ++nbounds;
    iisVars_[it->second] = (int)(lb && ub ? IISStatus::fix
                                 : lb ? IISStatus::low
                                 : ub ? IISStatus::upp : IISStatus::non);
  }
  AddToSolverMessage(fmt::format(
      "IIS with {} constraints and {} variable bounds{}.\n", ncons, nbounds,
      SCIPiisIsSubscipIrreducible(iis) ? "" : " (not proven irreducible)"));
}

IIS ScipBackend::GetIIS() {
  if (iisVars_.empty() && iisCons_.empty())
    return {};
  auto mv = GetValuePresolver().PostsolveIIS(
        { iisVars_, {{{ CG_Linear, iisCons_ }}} });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

void ScipBackend::DoWriteProblem(const std::string& name) {
  ExportModel(name);
}

void ScipBackend::DoWriteSolution(const std::string& name) {
  std::unique_ptr<FILE, int(*)(FILE*)> file(std::fopen(name.c_str(), "w"), std::fclose);
  if (!file)
    MP_RAISE(fmt::format("Cannot open '{}' for writing the solution", name));
  SCIP_CCALL( SCIPprintBestSol(getSCIP(), file.get(), FALSE) );
}

} // namespace mp


// AMPLs
AMPLS_MP_Solver* Open_scip(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::ScipBackend()},
    cb);
}

void AMPLSClose_scip(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_scip(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::ScipBackend*>(AMPLSGetBackend(slv))->getSCIP();
}

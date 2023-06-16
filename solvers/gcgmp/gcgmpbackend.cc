#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "gcgmpbackend.h"

#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.hpp"

extern "C" {
  #include "gcgmp-ampls-c-api.h"    // Gcg AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptGcg(void* prob) {
  SCIPinterruptSolve(static_cast<SCIP*>(prob));
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateGcgBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::GcgBackend()};
}


namespace mp {

/// Create Gcg Model Manager
/// @param gc: the Gcg common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return GcgModelMgr
std::unique_ptr<BasicModelManager>
CreateGcgModelMgr(GcgCommon&, Env&, pre::BasicValuePresolver*&);


GcgBackend::GcgBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateGcgModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

GcgBackend::~GcgBackend() {
  CloseSolver();
}


const char* GcgBackend::GetBackendName()
  { return "GGCGBackend"; }

std::string GcgBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", GCGmajorVersion(), 
    GCGminorVersion(), GCGtechVersion());
}


bool GcgBackend::IsMIP() const {
  return SCIPgetNOrigIntVars(getSCIP()) > 0 || SCIPgetNOrigBinVars(getSCIP()) > 0;
}

bool GcgBackend::IsQCP() const {
  return false;
}

ArrayRef<double> GcgBackend::PrimalSolution() {
  SCIP* scip = getSCIP();
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  for (int i = 0; i < num_vars; i++)
    x[i] = SCIPgetSolVal(scip, SCIPgetBestSol(scip), getPROBDATA()->vars[i]);
  return x;
}

pre::ValueMapDbl GcgBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> GcgBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
 // int error = GCG_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 0;
  if (error)
    pi.clear();
  return pi;
}

double GcgBackend::ObjectiveValue() const {
  return SCIPgetPrimalbound(getSCIP());
}

double GcgBackend::NodeCount() const {
  return SCIPgetNNodes(getSCIP());
}

double GcgBackend::SimplexIterations() const {
  return SCIPgetNPrimalLPIterations(getSCIP()) + SCIPgetNDualLPIterations(getSCIP());
}

int GcgBackend::BarrierIterations() const {
  return SCIPgetNBarrierLPIterations(getSCIP());
}

void GcgBackend::ExportModel(const std::string &file) {
  GCG_CCALL( SCIPwriteOrigProblem(getSCIP(), file.data(), NULL, FALSE) );
}


void GcgBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptGcg, getSCIP());
  // TODO Check interrupter
  //GCG_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void GcgBackend::Solve() {
  if (!storedOptions_.exportFile_.empty())
    ExportModel(storedOptions_.exportFile_);
  if (!storedOptions_.paramRead_.empty())
    GCG_CCALL( SCIPreadParams(getSCIP(), storedOptions_.paramRead_.c_str()) );
  if (storedOptions_.heuristics_ != 0)
    GCG_CCALL( SCIPsetHeuristics(getSCIP(), (SCIP_PARAMSETTING)storedOptions_.heuristics_, TRUE) );
  if (storedOptions_.cuts_ != 0)
    GCG_CCALL( SCIPsetSeparating(getSCIP(), (SCIP_PARAMSETTING)storedOptions_.cuts_, TRUE) );
  if (storedOptions_.presolvings_ != 0)
    GCG_CCALL( SCIPsetSeparating(getSCIP(), (SCIP_PARAMSETTING)storedOptions_.presolvings_, TRUE) );

  InputDecomposition();
  
  GCG_CCALL(GCGsolve(getSCIP()));

  WindupGCGSolve();
}

void GcgBackend::WindupGCGSolve() { }

void GcgBackend::ReportResults() {
  ReportGCGResults();
  BaseBackend::ReportResults();
}

void GcgBackend::ReportGCGResults() {
  SetStatus( ConvertGCGStatus() );
  AddGCGMessages();
  if (need_multiple_solutions())
    ReportGCGPool();
}
std::vector<double> GcgBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // GCG_CCALL(GCG_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double GcgBackend::getPoolObjective(int i)
{
  double obj;
 // GCG_CCALL(GCG_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void GcgBackend::ReportGCGPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(GCG_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void GcgBackend::AddGCGMessages() {
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

std::pair<int, std::string> GcgBackend::ConvertGCGStatus() {
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


void GcgBackend::FinishOptionParsing() {
  if (storedOptions_.outlev_ == 1)
     SetSolverOption("display/verblevel", 4);
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

static const mp::OptionValueInfo childsel[] = {
  { "d", "down", 0},
  { "u", "up", 1},
  { "p", "pseudo costs", 2},
  { "i", "inference", 3},
  { "l", "lp value", 4},
  { "r", "root LP value difference", 5},
  { "h", "hybrid inference/root LP value difference (default)", 6}
};

static const mp::OptionValueInfo scoretype[] = {
  { "maxwhi", "max white", 0},
  { "border", "border area", 1},
  { "classi", "classic", 2},
  { "forswh", "max foreseeing white", 3},
  { "spfwh", "ppc-max-white (default)", 4},
  { "fawh", "max foreseeing white with aggregation info", 5},
  { "spfawh", "ppc-max-white with aggregation info", 6},
  { "bender", "experimental benders score", 7},
  { "strong", "strong decomposition score", 8}
};

static const mp::OptionValueInfo mode[] = {
  { "0", "Dantzig-Wolfe (default)", 0},
  { "1", "Benders' decomposition", 1},
  { "2", "No decomposition", 2}
};

void GcgBackend::InitCustomOptions() {

  set_option_header(
      "GCG Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``gcg_options``. For example::\n"
      "\n"
      "  ampl: option gcg_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:outlev outlev",
    "0*/1: Whether to write GCG log lines (chatter) to stdout and to file.",
    storedOptions_.outlev_);

  AddSolverOption("tech:outlev-native outlev-native",
    "0*/1/2/3/4/5: Whether to write GCG log lines (chatter) to stdout and to file (native output level of SCIP).",
    "display/verblevel", 0, 5);

  AddStoredOption("tech:exportfile writeprob writemodel",
    "Specifies the name of a file where to export the model before "
    "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
    "Default = \"\" (don't export the model).",
    storedOptions_.exportFile_);

  AddStoredOption("tech:logfile logfile",
    "Log file name; note that the solver log will be written to the log "
    "regardless of the value of tech:outlev.",
    storedOptions_.logFile_);

  AddStoredOption("tech:param:read param:read paramfile",
    "Filename of GCG parameter file (as path)."
    "The suffix on a parameter file should be .set.\n",
    storedOptions_.paramRead_);

  AddSolverOption("alg:method method lpmethod",
    "LP algorithm for solving initial LP relaxations:\n"
    "\n.. value-table::\n", "lp/initalgorithm", lp_values_method, "s");

  AddSolverOption("alg:remethod remethod relpmethod",
    "LP algorithm for resolving LP relaxations if a starting basis exists:\n"
    "\n.. value-table::\n", "lp/resolvealgorithm", lp_values_method, "s");

  /////////////////////// BRANCHING ///////////////////////
  AddSolverOption("branch:preferbinary preferbinary",
    "0/1: whether branching on binary variables should be preferred"
    "\n"
    "  | 0 - Binary variables should not be preferred on branching (default)\n"
    "  | 1 - Binary variables should be preferred on branching.",
    "branching/preferbinary", 0, 1);

  //////////////////////// CUTS //////////////////////////
  AddSolverOption("cut:dircutoffdistweight dircutoffdistweight",
    "Weight of directed cutoff distance in cut score calculation (default: 0.0)",
    "cutselection/hybrid/dircutoffdistweight", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:efficacyweight efficacyweight",
    "Weight of efficacy in cut score calculation (default: 1.0)",
    "cutselection/hybrid/efficacyweight", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:intsupportweight intsupportweight",
    "Weight of integral support in cut score calculation (default: 0.1)",
    "cutselection/hybrid/intsupportweight", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:minortho minortho",
    "Minimal orthogonality for a cut to enter the LP (default: 0.9)",
    "cutselection/hybrid/minortho", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:minorthoroot minorthoroot",
    "Minimal orthogonality for a cut to enter the LP in the root node (default: 0.1)",
    "cutselection/hybrid/minorthoroot", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:objparalweight objparalweight",
    "Weight of objective parallelism in cut score calculation (default: 0.1)",
    "cutselection/hybrid/objparalweight", 0.0, SCIP_REAL_MAX);


  AddSolverOption("cut:maxcuts maxcuts",
    "Maximal number of cuts separated per separation round (0: disable local separation; default: 100)",
    "separating/maxcuts", 0, INT_MAX);

  AddSolverOption("cut:maxcutsroot maxcutsroot",
    "Maximal number of separated cuts at the root node (0: disable root node separation; default: 2000)",
    "separating/maxcutsroot", 0, INT_MAX);

  AddSolverOption("cut:maxrounds",
    "Maximal number of separation rounds per node (default: -1: unlimited)",
    "separating/maxrounds", -1, INT_MAX);

  AddSolverOption("cut:maxroundsroot",
    "Maximal number of separation rounds in the root node (default: -1: unlimited)",
    "separating/maxroundsroot", -1, INT_MAX);

  AddSolverOption("cut:maxstallrounds",
    "Maximal number of consecutive separation rounds without objective or integrality improvement in local nodes (-1: no additional restriction; default: 1)",
    "separating/maxstallrounds", -1, INT_MAX);

  AddSolverOption("cut:maxstallroundsroot",
    "Maximal number of consecutive separation rounds without objective or integrality improvement in the root node (-1: no additional restriction; default: 10)",
    "separating/maxstallroundsroot", -1, INT_MAX);

  AddSolverOption("cut:minactivityquot",
    "Minimum cut activity quotient to convert cuts into constraints during a restart (0.0: all cuts are converted; default: 0.8)",
    "separating/minactivityquot", 0.0, 1.0);

  AddSolverOption("cut:minefficacy",
    "Minimal efficacy for a cut to enter the LP (default: 0.0001)",
    "separating/minefficacy", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:minefficacyroot",
    "Minimal efficacy for a cut to enter the LP in the root node (default: 0.0001)",
    "separating/minefficacyroot", 0.0, SCIP_REAL_MAX);

  AddSolverOption("cut:poolfreq",
    "Separation frequency for the global cut pool (-1: disable global cut pool; 0: only separate pool at the root; default: 10)",
    "separating/poolfreq", -1, SCIP_MAXTREEDEPTH);

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
  AddSolverOption("lim:absgap absgap",
    "Solving stops, if the absolute gap = |primalbound - dualbound| is below the given value (default: 0.0)",
    "limits/absgap", 0.0, SCIP_REAL_MAX);

  AddSolverOption("lim:autorestartnodes",
    "If solve exceeds this number of nodes for the first time, an automatic restart is triggered (default: -1: no automatic restart)",
    "limits/autorestartnodes", -1, INT_MAX);

  AddSolverOption("lim:bestsol",
    "Solving stops, if the given number of solution improvements were found (default: -1: no limit)",
    "limits/bestsol", -1, INT_MAX);

  AddSolverOption("lim:gap gap",
    "Solving stops, if the relative gap = |primal - dual|/MIN(|dual|,|primal|) is below the given value, the gap is 'Infinity', if primal and dual bound have opposite signs (default: 0.0)",
    "limits/gap", 0.0, SCIP_REAL_MAX);

  AddSolverOption("lim:maxorigsol",
    "Maximal number of solutions candidates to store in the solution storage of the original problem (default: 10)",
    "limits/maxorigsol", 0, INT_MAX);

  AddSolverOption("lim:maxsol",
    "Maximal number of solutions to store in the solution storage (default: 100)",
    "limits/maxsol", 1, INT_MAX);

  AddSolverOption("lim:memory memory",
    "#maximal memory usage in MB; reported memory usage is lower than real memory usage! (default: 8796093022207.0)",
    "limits/memory", 0.0, (SCIP_Real)SCIP_MEM_NOLIMIT);

  AddSolverOption("lim:nodes",
    "Maximal number of nodes to process (default: -1: no limit)",
    "limits/nodes", -1, (int)SCIP_LONGINT_MAX);

  AddSolverOption("lim:restarts",
    "Solving stops, if the given number of restarts was triggered (default: -1: no limit)",
    "limits/restarts", -1, INT_MAX);

  AddSolverOption("lim:solutions",
    "Solving stops, if the given number of solutions were found (default: -1: no limit)",
    "limits/solutions", -1, INT_MAX);

  AddSolverOption("lim:stallnodes",
    "Solving stops, if the given number of nodes was processed since the last improvement of the primal solution value (default: -1: no limit)",
    "limits/stallnodes", -1, (int)SCIP_LONGINT_MAX);

  AddSolverOption("lim:time timelim timelimit time_limit",
    "Limit on solve time (in seconds; default: 1e+20).",
    "limits/time", 0.0, 1e+20);

  AddSolverOption("lim:totalnodes",
    "Maximal number of total nodes (incl. restarts) to process (default: -1: no limit)",
    "limits/totalnodes", -1, (int)SCIP_LONGINT_MAX);

  ////////////////////////// LP //////////////////////////
  AddSolverOption("lp:pricing pricing",
    "Pricing strategy:\n"
    "\n.. value-table::",
    "lp/pricing", values_pricing, "l");

  AddSolverOption("lp:presolving",
    "0/1: whether presolving of LP solver should be used"
    "\n"
    "  | 0 - Presolving of LP solver should not be used\n"
    "  | 1 - Presolving of LP solver should be used (default).",
    "lp/advanced/presolving", 0, 1);

  AddSolverOption("lp:threads",
    "Number of threads used for solving the LP (default: 0: automatic)",
    "lp/advanced/threads", 0, 64);

  AddSolverOption("lp:alwaysgetduals alwaysgetfarkasduals alwaysgetduals",
    "0/1: whether the Farkas duals should always be collected when an LP is found to be infeasible"
    "\n"
    "  | 0 - The Farkas duals should not always be collected when an LP is found to be infeasible (default)\n"
    "  | 1 - The Farkas duals should always be collected when an LP is found to be infeasible.",
    "lp/alwaysgetduals", 0, 1);

  AddSolverOption("lp:solvedepth",
    "Maximal depth for solving LP at the nodes (default: -1: no depth limit)",
    "lp/solvedepth", -1, SCIP_MAXTREEDEPTH);

  AddSolverOption("lp:solvefreq",
    "Frequency for solving LP at the nodes (-1: never; 0: only root LP; default: 0)",
    "lp/solvefreq", -1, SCIP_MAXTREEDEPTH);

  ///////////////////////// MISC /////////////////////////
  AddSolverOption("misc:allowstrongdualreds allowstrongdualreds",
    "0/1: whether strong dual reductions should be allowed in propagation and presolving"
    "\n"
    "  | 0 - Strong dual reductions should not be allowed in propagation and presolving\n"
    "  | 1 - Strong dual reductions should be allowed in propagation and presolving (default).",
    "misc/allowstrongdualreds", 0, 1);

  AddSolverOption("misc:allowweakdualreds allowweakdualreds",
    "0/1: whether weak dual reductions should be allowed in propagation and presolving"
    "\n"
    "  | 0 - Weak dual reductions should not be allowed in propagation and presolving\n"
    "  | 1 - Weak dual reductions should be allowed in propagation and presolving (default).",
    "misc/allowweakdualreds", 0, 1);

  AddSolverOption("misc:scaleobj scaleobj",
    "0/1: whether the objective function should be scaled so that it is always integer"
    "\n"
    "  | 0 - The objective function should not be scaled so that it is always integer\n"
    "  | 1 - The objective function should be scaled so that it is always integer (default).",
    "misc/scaleobj", 0, 1);

  ////////////////////// NUMERICS ////////////////////////
  AddSolverOption("num:checkfeastolfac checkfeastolfac",
    "Feasibility tolerance factor; for checking the feasibility of the best solution (default: 1.0)",
    "numerics/checkfeastolfac", 0.0, SCIP_REAL_MAX);

  AddSolverOption("num:dualfeastol dualfeastol",
    "Feasibility tolerance for reduced costs in LP solution (default: 1e-07)",
    "numerics/dualfeastol", SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON);

  AddSolverOption("num:epsilon epsilon",
    "Absolute values smaller than this are considered zero (default: 1e-09)",
    "numerics/epsilon", SCIP_MINEPSILON, SCIP_MAXEPSILON);

  AddSolverOption("num:feastol feastol",
    "Feasibility tolerance for constraints (default: 1e-06)",
    "numerics/feastol", SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON);

  AddSolverOption("num:infinity infinity",
    "Values larger than this are considered infinity (default: 1e+20)",
    "numerics/infinity", 1e+10, SCIP_INVALID/10.0);

  AddSolverOption("num:lpfeastolfactor lpfeastolfactor",
    "Factor w.r.t. primal feasibility tolerance that determines default (and maximal) primal feasibility tolerance of LP solver (default: 1.0)",
    "numerics/lpfeastolfactor", 1e-6, 1.0);

  AddSolverOption("num:sumepsilon sumepsilon",
    "Absolute values of sums smaller than this are considered (default: 1e-06)",
    "numerics/sumepsilon", SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON);

  //////////////////// NODESELECTION /////////////////////
  AddSolverOption("nod:childsel",
    "Child selection rule:\n"
    "\n.. value-table::\n", "nodeselection/childsel", childsel, "h");

  ////////////////////// PARALLEL ////////////////////////
  AddSolverOption("par:maxnthreads maxnthreads",
    "Maximum number of threads used during parallel solve (default: 8)",
    "parallel/maxnthreads", 0, 64);

  AddSolverOption("par:minnthreads minnthreads",
    "Minimum number of threads used during parallel solve (default: 1)",
    "parallel/minnthreads", 0, 64);

  AddSolverOption("par:mode mode",
    "0/1: Parallel optimisation mode"
    "\n"
    "  | 0 - Opportunistic\n"
    "  | 1 - Deterministic (default)",
    "parallel/mode", 0, 1);

  ////////////////////// PRESOLVE ////////////////////////
  AddSolverOption("pre:abortfac abortfac",
    "Abort presolve, if at most this fraction of the problem was changed in last presolve round (default: 0.0008)",
    "presolving/advanced/abortfac", 0.0, 1.0);

  AddSolverOption("pre:clqtablefac clqtablefac",
    "Limit on number of entries in clique table relative to number of problem nonzeros (default: 2.0)",
    "presolving/advanced/clqtablefac", 0.0, SCIP_REAL_MAX);

  AddSolverOption("pre:donotaggr donotaggr",
    "0/1: whether aggregation of variables should be forbidden"
    "\n"
    "  | 0 - Aggregation of variables should not be forbidden (default)\n"
    "  | 1 - Aggregation of variables should be forbidden.",
    "presolving/advanced/donotaggr", 0, 1);

  AddSolverOption("pre:donotmultaggr donotmultaggr",
    "0/1: whether multi-aggregation of variables should be forbidden"
    "\n"
    "  | 0 - Multi-aggregation of variables should not be forbidden (default)\n"
    "  | 1 - Multi-aggregation of variables should be forbidden.",
    "presolving/advanced/donotmultaggr", 0, 1);

  AddSolverOption("pre:immrestartfac immrestartfac",
    "Fraction of integer variables that were fixed in the root node triggering an immediate restart with preprocessing (default: 0.1)",
    "presolving/advanced/immrestartfac", 0.0, 1.0);

  AddSolverOption("pre:restartfac restartfac",
    "Fraction of integer variables that were fixed in the root node triggering a restart with preprocessing after root node evaluation (default: 0.025)",
    "presolving/advanced/restartfac", 0.0, 1.0);

  AddSolverOption("pre:restartminred restartminred",
    "Minimal fraction of integer variables removed after restart to allow for an additional restart (default: 0.1)",
    "presolving/advanced/restartminred", 0.0, 1.0);

  AddSolverOption("pre:subrestartfac subrestartfac",
    "Fraction of integer variables that were globally fixed during the solving process triggering a restart with preprocessing (default: 1.0)",
    "presolving/advanced/subrestartfac", 0.0, 1.0);

  AddSolverOption("pre:maxrestarts",
    "Maximal number of restarts (default: -1: unlimited)",
    "presolving/maxrestarts", -1, INT_MAX);

  AddSolverOption("pre:maxrounds",
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

  ////////////////////// PROPAGATE ///////////////////////
  AddSolverOption("pro:abortoncutoff",
    "0/1: whether propagation should be aborted immediately (setting this to 0 could help conflict analysis to produce more conflict constraints)"
    "\n"
    "  | 0 - Propagation should not be aborted immediately\n"
    "  | 1 - Propagation should be aborted immediately (default).",
    "propagating/abortoncutoff", 0, 1);

  AddSolverOption("pro:maxrounds",
    "Maximal number of propagation rounds per node (-1: unlimited; 0: off; default: 100)",
    "propagating/maxrounds", -1, INT_MAX);

  AddSolverOption("pro:maxroundsroot",
    "Maximal number of propagation rounds in the root node (-1: unlimited; 0: off; default: 1000)",
    "propagating/maxroundsroot", -1, INT_MAX);
  
  //////////////////// RANDOMIZATION /////////////////////
  AddSolverOption("ran:permuteconss permuteconss",
    "0/1: whether the order of constraints should be permuted (depends on permutationseed)? "
    "\n"
    "  | 0 - Order of constraints should not be permuted\n"
    "  | 1 - Order of constraints should be permuted (default).",
    "randomization/advanced/permuteconss", 0, 1);

  AddSolverOption("ran:permutevars permutevars",
    "0/1: whether the order of variables should be permuted (depends on permutationseed)? "
    "\n"
    "  | 0 - Order of variables should not be permuted (default)\n"
    "  | 1 - Order of variables should be permuted.",
    "randomization/advanced/permutevars", 0, 1);

  AddSolverOption("ran:lpseed lpseed",
    "Random seed for LP solver, e.g. for perturbations in the simplex (default: 0: LP default)",
    "randomization/lpseed", 0, INT_MAX);

  AddSolverOption("ran:permutationseed permutationseed",
    "Seed value for permuting the problem after reading/transformation (default: 0: no permutation) ",
    "randomization/permutationseed", 0, INT_MAX);

  AddSolverOption("ran:randomseedshift randomseedshift",
    "Global shift of all random seeds in the plugins and the LP random seed (default: 0) ",
    "randomization/randomseedshift", 0, INT_MAX);


  ////////////////////// DETECTION ///////////////////////
  AddSolverOption("det:enabled",
    "0/1: whether detection should be enabled? "
    "\n"
    "  | 0 - Detection not enabled\n"
    "  | 1 - Detection enabled (default).",
    "detection/enabled", 0, 1);

  AddSolverOption("det:postprocess postprocess",
    "0/1: whether postprocessing of complete decompositions should be enabled? "
    "\n"
    "  | 0 - Postprocessing of complete decompositions not enabled\n"
    "  | 1 - Postprocessing of complete decompositions enabled (default).",
    "detection/enabled", 0, 1);

  AddSolverOption("det:maxrounds",
    "Maximum number of detection loop rounds (default: 1) ",
    "detection/maxrounds", 0, INT_MAX);

  AddSolverOption("det:maxtime",
    "Maximum detection time in seconds (default: 600) ",
    "detection/maxrounds", 0, INT_MAX);

  AddSolverOption("det:scoretype scoretype",
    "Score calculation for comparing (partial) decompositions:\n"
    "\n.. value-table::\n", "detection/scores/selected", scoretype, 4);

  AddSolverOption("det:origprob-classificationenabled origprob-classificationenabled",
    "0/1: whether classification for the original problem should be enabled? "
    "\n"
    "  | 0 - Classification for the original problem not enabled\n"
    "  | 1 - Classification for the original problem enabled (default).",
    "detection/origprob/classificationenabled", 0, 1);

  AddSolverOption("det:origprob-enabled origprob-enabled",
    "0/1: whether detection for the original problem should be enabled? "
    "\n"
    "  | 0 - Detection for the original problem not enabled\n"
    "  | 1 - Detection for the original problem enabled (default).",
    "detection/origprob/enabled", 0, 1);

  AddSolverOption("det:classification-enabled classification-enabled",
    "0/1: whether classification should be enabled? "
    "\n"
    "  | 0 - Classification not enabled\n"
    "  | 1 - Classification enabled (default).",
    "detection/classification/enabled", 0, 1);

  AddSolverOption("det:maxnclassesperpartition maxnclassesperpartition",
    "Maximum number of classes per partition (default: 9)",
    "detection/classification/maxnclassesperpartition", 0, INT_MAX);

  AddSolverOption("det:maxnclassesperpartitionforlargeprobs maxnclassesperpartitionforlargeprobs",
    "Maximum number of classes per partition for large problems (nconss + nvars >= 50000) (default: 5)",
    "detection/classification/maxnclassesperpartitionforlargeprobs", 0, INT_MAX);
  
  AddSolverOption("det:benders-enabled benders-enabled",
    "0/1: whether benders detection should be enabled? "
    "\n"
    "  | 0 - Benders detection not enabled (default)\n"
    "  | 1 - Benders detection enabled.",
    "detection/benders/enabled", 0, 1);

  AddSolverOption("det:benders-onlybinmaster benders-onlybinmaster",
    "0/1: whether only decomposition with only binary variables in the master are searched? "
    "\n"
    "  | 0 - Not only decomposition with only binary variables in the master are searched (default)\n"
    "  | 1 - Only decomposition with only binary variables in the master are searched.",
    "detection/benders/onlybinmaster", 0, 1);

  AddSolverOption("det:benders-onlycontsubpr benders-onlycontsubpr",
    "0/1: whether only decomposition with only continiuous variables in the subproblems are searched? "
    "\n"
    "  | 0 - Not only decomposition with only continiuous variables in the subproblems are searched (default)\n"
    "  | 1 - Only decomposition with only continiuous variables in the subproblems are searched.",
    "detection/benders/onlycontsubpr", 0, 1);

  //////////////////// RELAXING-GCG //////////////////////
  AddSolverOption("gcg:bliss-enabled bliss-enabled",
    "0/1: whether bliss should be used to check for identical blocks? "
    "\n"
    "  | 0 - Bliss should not be used to check for identical blocks\n"
    "  | 1 - Bliss should be used to check for identical blocks (default).",
    "relaxing/gcg/bliss/enabled", 0, 1);

  AddSolverOption("gcg:aggregation",
    "0/1: whether identical blocks should be aggregated (only for discretization approach)? "
    "\n"
    "  | 0 - Identical blocks should not be aggregated (only for discretization approach)\n"
    "  | 1 - Identical blocks should be aggregated (only for discretization approach) (default).",
    "relaxing/gcg/aggregation", 0, 1);

  AddSolverOption("gcg:discretization discretization",
    "0/1: whether discretization (TRUE) or convexification (FALSE) approach should be used? "
    "\n"
    "  | 0 - convexification approach should be used\n"
    "  | 1 - discretization approach should be used (default).",
    "relaxing/gcg/discretization", 0, 1);

  AddSolverOption("gcg:mipdiscretization mipdiscretization",
    "0/1: whether discretization (TRUE) or convexification (FALSE) approach should be used in mixed-integer programs?"
    "\n"
    "  | 0 - convexification approach should be used in mixed-integer programs\n"
    "  | 1 - discretization approach should be used in mixed-integer programs (default).",
    "relaxing/gcg/discretization", 0, 1);

  AddSolverOption("gcg:mode mode",
    "The decomposition mode that GCG will use:\n"
    "\n.. value-table::\n", "relaxing/gcg/mode", mode, 0);
}

void GcgBackend::InputDecomposition() {
  bool is_presolved = SCIPgetStage(getSCIP()) >= SCIP_STAGE_PRESOLVED;
  GCG_CCALL( SCIPallocClearMemory(scip, &getPROBDATA()->decomp) );
  gcg::PARTIALDECOMP* decomp = new gcg::PARTIALDECOMP(getSCIP(), !is_presolved);
  getPROBDATA()->decomp = decomp;

  if (auto block0 = ReadModelSuffixInt({"block", suf::Kind::CON_BIT | suf::Kind::VAR_BIT})) {
    auto block = GetValuePresolver().PresolveGenericInt( block0 );

    ArrayRef<int> blockconss = block.GetConValues()(CG_Linear);
    for (size_t i = 0; i < blockconss.size(); i++) {
      if (blockconss[i] != 0)
        getPROBDATA()->decomp->fixConsToBlock(getPROBDATA()->linconss[i], blockconss[i]-1);
    }

    ArrayRef<int> blockvars = block.GetVarValues()();
    for (size_t i = 0; i < blockvars.size(); i++) {
      if (blockvars[i] != 0) {
        int varindex = getPROBDATA()->decomp->getDetprobdata()->getIndexForVar(getPROBDATA()->vars[i]);
        getPROBDATA()->decomp->fixVarToBlock(varindex, blockvars[i]);
      }
    }
  }

  if (auto master0 = ReadModelSuffixInt({"master", suf::Kind::CON_BIT | suf::Kind::VAR_BIT})) {
    auto master = GetValuePresolver().PresolveGenericInt( master0 );

    ArrayRef<int> masterconss = master.GetConValues()(CG_Linear);
    for (size_t i = 0; i < masterconss.size(); i++) {
      if (masterconss[i] == 1)
        getPROBDATA()->decomp->fixConsToMaster(getPROBDATA()->linconss[i]);
    }

    ArrayRef<int> mastervars = master.GetVarValues()();
    for (size_t i = 0; i < mastervars.size(); i++) {
      if (mastervars[i] == 1) {
        int varindex = getPROBDATA()->decomp->getDetprobdata()->getIndexForVar(getPROBDATA()->vars[i]);
        getPROBDATA()->decomp->fixVarToMaster(varindex);
      }
    }
  }

  if (auto linking0 = ReadModelSuffixInt({"linking", suf::Kind::VAR_BIT})) {
    auto linking = GetValuePresolver().PresolveGenericInt( linking0 );

    ArrayRef<int> linkingvars = linking.GetConValues()(CG_Linear);
    for (size_t i = 0; i < linkingvars.size(); i++) {
      if (linkingvars[i] == 1) {
        int varindex = getPROBDATA()->decomp->getDetprobdata()->getIndexForVar(getPROBDATA()->vars[i]);
        getPROBDATA()->decomp->fixVarToLinking(varindex);
      }
    }
  }

  if (getPROBDATA()->decomp->getNOpenconss() != getPROBDATA()->decomp->getNConss() || getPROBDATA()->decomp->getNOpenvars() != getPROBDATA()->decomp->getNVars()) {
    getPROBDATA()->decomp->prepare();
    if (getPROBDATA()->decomp->isComplete())
      getPROBDATA()->decomp->setUsergiven(gcg::USERGIVEN::COMPLETE);
    else
      getPROBDATA()->decomp->setUsergiven(gcg::USERGIVEN::PARTIAL);
    GCGconshdlrDecompAddPreexisitingPartialDec(getSCIP(), getPROBDATA()->decomp);
    GCGpresolve(getSCIP());
  }
}

double GcgBackend::MIPGap() {
  return SCIPgetGap(getSCIP())<Infinity() ? SCIPgetGap(getSCIP()) : AMPLInf();
}
double GcgBackend::BestDualBound() {
  return SCIPgetDualbound(getSCIP());
}

double GcgBackend::MIPGapAbs() {
  double gapabs = std::fabs(ObjectiveValue() - BestDualBound());
  return gapabs<Infinity() ? gapabs : AMPLInf();
}


ArrayRef<int> GcgBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  GCG_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case GCG_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case GCG_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case GCG_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case GCG_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case GCG_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Gcg VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> GcgBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  GCG_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case GCG_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case GCG_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case GCG_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case GCG_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case GCG_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Gcg VBasis value: {}", s));
    }
  }*/
  return cons;
}

void GcgBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = GCG_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = GCG_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = GCG_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = GCG_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = GCG_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!GCG_GetColInfo(lp(), GCG_DBLINFO_LB, 1, index, &lb) && 
        !GCG_GetColInfo(lp(), GCG_DBLINFO_UB, 1, index, &ub))
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
  GCG_SetBasis(lp(), stt.data(), NULL);
  */
}

void GcgBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = GCG_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = GCG_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  GCG_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis GcgBackend::GetBasis() {
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

void GcgBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void GcgBackend::ComputeIIS() {
  //GCG_CCALL(GCG_ComputeIIS(lp()));
  SetStatus(ConvertGCGStatus());   // could be new information
}

IIS GcgBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> GcgBackend::VarsIIS() {
  return ArrayRef<int>();
//  return getIIS(lp(), NumVars(), GCG_GetColLowerIIS, GCG_GetColUpperIIS);
}
pre::ValueMapInt GcgBackend::ConsIIS() {
  /*auto iis_lincon = getIIS(lp(), NumLinCons(), GCG_GetRowLowerIIS, GCG_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  GCG_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  GCG_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };
      */
  return { {{ 0, std::vector<int>()}} };
}

void GcgBackend::AddMIPStart(ArrayRef<double> x0, ArrayRef<int> sparsity) {
  //GCG_CCALL(GCG_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));

}


} // namespace mp


// AMPLs
AMPLS_MP_Solver* Open_gcg(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::GcgBackend()},
    cb);
}

void AMPLSClose_gcg(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_gcg(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::GcgBackend*>(AMPLSGetBackend(slv))->getSCIP();
}

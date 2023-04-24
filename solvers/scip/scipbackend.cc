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
}

void ScipBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }

  if (storedOptions_.concurrent_)
    SCIP_CCALL( SCIPsolveConcurrent(getSCIP()) );
  else
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

  AddStoredOption("alg:concurrent concurrent",
    "0/1: whether to solve the problem using concurrent solvers"
    "\n"
    "  | 0 - No concurrent solvers are used to solve the problem (default)\n"
    "  | 1 - Concurrent solvers are used to solve the problem.",
    storedOptions_.concurrent_, 0, 1);

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
    "Maximal number of separated cuts at the root node (0: disable root node separation; default: 5000)",
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

  AddSolverOption("lim:softtime softtime",
    "Soft time limit which should be applied after first solution was found (default: -1.0: disabled)",
    "limits/softtime", -1.0, SCIP_REAL_MAX);

  AddSolverOption("lim:solutions",
    "Solving stops, if the given number of solutions were found (default: -1: no limit)",
    "limits/solutions", -1, INT_MAX);

  AddSolverOption("lim:stallnodes",
    "Solving stops, if the given number of nodes was processed since the last improvement of the primal solution value (default: -1: no limit)",
    "limits/stallnodes", -1, (int)SCIP_LONGINT_MAX);

  AddSolverOption("lim:time time",
    "Maximal time in seconds to run (default: 1e+20)",
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
    "Frequency for solving LP at the nodes (-1: never; 0: only root LP; default: 1)",
    "lp/solvefreq", -1, SCIP_MAXTREEDEPTH);

  /////////////////////// PRESOLVE ///////////////////////
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

#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "coptbackend.h"

extern "C" {
  #include "copt-ampls-c-api.h"    // Copt AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptCopt(void* prob) {
  return COPT_Interrupt((copt_prob*)prob);
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCoptBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CoptBackend()};
}


namespace mp {

/// Create Copt Model Manager
/// @param gc: the Copt common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return CoptModelMgr
std::unique_ptr<BasicModelManager>
CreateCoptModelMgr(CoptCommon&, Env&, pre::BasicValuePresolver*&);


void CoptBackend::OpenSolver() {
  int status = 0;
  const auto create_fn = GetCallbacks().init;
  if (create_fn)
    set_env((copt_env*)create_fn());
  else
    COPT_CCALL(COPT_CreateEnv(&env_ref()));
  if (env() == NULL) {
    throw std::runtime_error(
      fmt::format("Could not open COPT environment.\n{}", status));
  }
  
  /* Create an empty model */
  status = COPT_CreateProb(env(), &lp_ref());
  if (status)
    throw std::runtime_error(fmt::format(
      "Failed to create problem, error code {}.", status));
  COPT_CCALL(COPT_SetIntParam(lp(), "Logging", 0));
  /* Copy handlers to ModelAPI */
  copy_common_info_to_other();
  SetSolverOption(COPT_INTPARAM_LOGGING, 0);
}

void CoptBackend::CloseSolver() {
  if (lp() != NULL) {
    COPT_CCALL(COPT_DeleteProb(&lp_ref()));
  }
  if (env() != NULL) {
    COPT_CCALL(COPT_DeleteEnv(&env_ref()));
  }
}

CoptBackend::CoptBackend() {
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateCoptModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

CoptBackend::~CoptBackend() {
  CloseSolver();
}

const char* CoptBackend::GetBackendName()
  { return "CoptBackend"; }

std::string CoptBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", COPT_VERSION_MAJOR, 
    COPT_VERSION_MINOR, COPT_VERSION_TECHNICAL);
}


bool CoptBackend::IsMIP() const {
  return getIntAttr(COPT_INTATTR_ISMIP);
}

bool CoptBackend::IsQCP() const {
  return getIntAttr(COPT_INTATTR_QCONSTRS) > 0;
}
void CoptBackend::InitOptionParsing() {
  OpenSolver();
}


void CoptBackend::InputExtras() {
  MIPBackend::InputExtras();
  InputCOPTExtras();
}

void CoptBackend::InputCOPTExtras() {
  if (feasrelax())
    DoCOPTFeasRelax();
  if (need_ray_primal() || need_ray_dual())
    SetSolverOption(COPT_INTPARAM_REQFARKASRAY, 1);
}


void CoptBackend::DoCOPTFeasRelax() {
  int copttype;
  int map[] = { 0, 4, 2, 1, 5, 3 };
  copttype = map[feasrelax()-1];
  /// only relax linear constraints
  auto mv = GetValuePresolver().PresolveSolution({
                                              {},
                                              feasrelax().rhspen()
    });
  const auto& rhspen = mv.GetConValues()(CG_Linear);
  /// Account for new variables (TODO: proper presolve)
  std::vector<double> lbpen = feasrelax().lbpen();
  if (lbpen.size() && lbpen.size() < (size_t)NumVars())
    lbpen.resize(NumVars());
  std::vector<double> ubpen = feasrelax().ubpen();
  if (ubpen.size() && ubpen.size() < (size_t)NumVars())
    ubpen.resize(NumVars());

  COPT_CCALL(COPT_FeasRelax(lp(),
    (double*)data_or_null(lbpen),
    (double*)data_or_null(ubpen),
    (double*)data_or_null(rhspen),
    NULL));
}

ArrayRef<double> CoptBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  if (IsMIP()) 
    error = COPT_GetSolution(lp(), x.data());
  else
    error = COPT_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);


  if (error)
    x.clear();
  return x;
}

pre::ValueMapDbl CoptBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> CoptBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  if (!IsMIP()) {
    int error = COPT_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
    if (error)
      pi.clear();
  }
  return pi;
}

double CoptBackend::ObjectiveValue() const {
  if (IsMIP())
    return getDblAttr(COPT_DBLATTR_BESTOBJ);
  else
    return getDblAttr(COPT_DBLATTR_LPOBJVAL);
}

double CoptBackend::NodeCount() const {
  return getIntAttr(COPT_INTATTR_NODECNT);
}

double CoptBackend::SimplexIterations() const {
  return getIntAttr(COPT_INTATTR_SIMPLEXITER);
}

int CoptBackend::BarrierIterations() const {
  return getIntAttr(COPT_INTATTR_BARRIERITER);
}

inline bool ends_with(std::string const& value, std::string const& ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}
void CoptBackend::DoWriteProblem(const std::string& name) {
  if (ends_with(name, ".lp"))
    COPT_CCALL(COPT_WriteLp(lp(), name.c_str()));
  else if (ends_with(name, ".mps"))
    COPT_CCALL(COPT_WriteMps(lp(), name.c_str()));
  else
    throw std::runtime_error("Can only export '.lp' or '.mps' files.");
}


void CoptBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCopt, lp());
}

void CoptBackend::Solve() {
  COPT_CCALL(COPT_Solve(lp()));
  WindupCOPTSolve();
}

void CoptBackend::WindupCOPTSolve() { }

void CoptBackend::ReportResults() {
  ReportCOPTResults();
  BaseBackend::ReportResults();
}

void CoptBackend::ReportCOPTResults() {
  SetStatus( ConvertCOPTStatus() );
  AddCOPTMessages();
  if (need_multiple_solutions())
    ReportCOPTPool();
}
std::vector<double> CoptBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
  COPT_CCALL(COPT_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double CoptBackend::getPoolObjective(int i)
{
  double obj;
  COPT_CCALL(COPT_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void CoptBackend::ReportCOPTPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  
  while (++iPoolSolution < getIntAttr(COPT_INTATTR_POOLSOLS)) {
    auto mv = GetValuePresolver().PostsolveSolution(  // only single-objective with pool
          { { getPoolSolution(iPoolSolution) },
            {},                                       // no duals
            std::vector<double>{ getPoolObjective(iPoolSolution) } } );
    ReportIntermediateSolution(
          { mv.GetVarValues()(), mv.GetConValues()(), mv.GetObjValues()() });
  }
}


void CoptBackend::AddCOPTMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> CoptBackend::ConvertCOPTStatus() {
  namespace sol = mp::sol;
  if (IsMIP())
  {
    int optstatus = getIntAttr(COPT_INTATTR_MIPSTATUS);
    int solstatus
        = getIntAttr(COPT_INTATTR_HASMIPSOL)
        | getIntAttr(COPT_INTATTR_HASFEASRELAXSOL);
    switch (optstatus) {
    case COPT_MIPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case COPT_MIPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case COPT_MIPSTATUS_INF_OR_UNB:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case COPT_MIPSTATUS_UNBOUNDED:
      if (solstatus)
        return { sol::UNBOUNDED_FEAS, "unbounded problem, feasible solution" };
      return { sol::UNBOUNDED_NO_FEAS, "unbounded problem, no solution" };
    case COPT_MIPSTATUS_TIMEOUT:
    case COPT_MIPSTATUS_NODELIMIT:
    case COPT_MIPSTATUS_INTERRUPTED:
    case COPT_MIPSTATUS_UNSTARTED:
      if (solstatus)
        return { sol::LIMIT_FEAS, "interrupted, feasible solution" };
      return { sol::LIMIT_NO_FEAS, "interrupted, no solution" };
    case COPT_MIPSTATUS_UNFINISHED:
      return { sol::NUMERIC, "failure, numeric issues" };
    default:
      return { sol::UNKNOWN, "unknown" };
    }
  }
  else {
    int optstatus = getIntAttr(COPT_INTATTR_LPSTATUS);
    int solstatus
        = getIntAttr(COPT_INTATTR_HASLPSOL)
        | getIntAttr(COPT_INTATTR_HASFEASRELAXSOL);
    switch (optstatus) {
    case COPT_LPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case COPT_LPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case COPT_LPSTATUS_UNBOUNDED:
      if (solstatus)
        return { sol::UNBOUNDED_FEAS, "unbounded problem, feasible solution" };
      return { sol::UNBOUNDED_NO_FEAS, "unbounded problem, no solution" };
    case COPT_LPSTATUS_TIMEOUT:
    case COPT_LPSTATUS_INTERRUPTED:
    case COPT_LPSTATUS_UNSTARTED:
      if (solstatus)
        return { sol::LIMIT_FEAS, "interrupted, feasible solution" };
      return { sol::LIMIT_NO_FEAS, "interrupted, no solution" };
    case COPT_LPSTATUS_UNFINISHED:
      return { sol::NUMERIC, "failure, numeric issues" };
    case COPT_LPSTATUS_IMPRECISE:
      return { sol::UNCERTAIN, "solution is imprecise" };
    case COPT_LPSTATUS_NUMERICAL:
      if (solstatus)
        return { sol::UNCERTAIN, "solution returned but error likely" };
      return { sol::NUMERIC, "failure, numeric issues" };
    default:
      return { sol::UNKNOWN, "unknown" };
    }
  }
}


void CoptBackend::FinishOptionParsing() {
  int v = 1;
  if (!storedOptions_.logFile_.empty())
  {
    SetSolverOption(COPT_INTPARAM_LOGGING, 1);
    COPT_CCALL(COPT_SetLogFile(lp(), storedOptions_.logFile_.data()));
  }
  else
    GetSolverOption(COPT_INTPARAM_LOGGING, v);
  set_verbose_mode(v>0);
  
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo lp_values_method[] = {
  { "-1", "Automatic (default)", -1},
  { "1", "Dual simplex", 1},
  { "2", "Barrier", 2},
  { "3", "Crossover", 3},
  { "4", "Concurrent (simplex and barrier simultaneously)", 4},
  { "5", "Choose between simplex and barrier automatically", 5}
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

static const mp::OptionValueInfo values_iismethod[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Find smaller IIS", 0},
  { "1", "Find IIS quickly", 1},
};

static const mp::OptionValueInfo values_crossover[] = {
  {"-1", "Choose automatically - only run crossover when LP solution is not primal-dual feasible", -1},
  {"0", "No", 0},
  {"1", "Yes", 1}
};

void CoptBackend::InitCustomOptions() {

  set_option_header(
      "COPT Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``copt_options``. For example::\n"
      "\n"
      "  ampl: option copt_options 'mipgap=1e-6';\n");

  AddSolverOption("tech:outlev outlev",
      "0-1: output logging verbosity. "
      "Default = 0 (no logging).",
    COPT_INTPARAM_LOGGING, 0, 1);

  AddStoredOption("tech:logfile logfile",
    "Log file name.", storedOptions_.logFile_);

  AddSolverOption("lp:dualprice dualprice",
    "Specifies the dual simplex pricing algorithm:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DUALPRICE,
    lp_dualprices_values_, -1);

  AddSolverOption("lp:dualperturb dualperturb",
    "Whether to allow the objective function perturbation when using "
    "the dual simplex method:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DUALPERTURB,
    values_autonoyes_, -1);

  AddSolverOption("lp:barhomogeneous barhomogeneous",
    "Whether to use homogeneous self-dual form in barrier:\n"
    "\n.. value-table::\n", COPT_INTPARAM_BARHOMOGENEOUS,
    values_autonoyes_, -1);

  AddSolverOption("lp:barorder barorder",
    "Barrier ordering algorithm:\n"
    "\n.. value-table::\n", COPT_INTPARAM_BARORDER,
    lp_barorder_values_, -1);

  
  AddSolverOption("bar:iterlim BarIterLimit",
    "Limit on the number of barrier iterations (default 500).",
    COPT_INTPARAM_BARITERLIMIT, 0, INT_MAX);

  AddSolverOption("mip:cutlevel cutlevel",
    "Level of cutting-planes generation:\n"
    "\n.. value-table::\n", COPT_INTPARAM_CUTLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:intfeastol intfeastol inttol",
    "Feasibility tolerance for integer variables (default 1e-06).",
    COPT_DBLPARAM_INTTOL, 1e-9, 0.1);

  AddSolverOption("mip:rootcutlevel rootcutlevel",
    "Level of cutting-planes generation of root node:\n"
    "\n.. value-table::\n", COPT_INTPARAM_ROOTCUTLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:treecutlevel treecutlevel",
    "Level of cutting-planes generation of search tree:\n"
    "\n.. value-table::\n", COPT_INTPARAM_TREECUTLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:rootcutrounds rootcutrounds",
    "Rounds of cutting-planes generation of root node;\n" 
    "default -1 ==> automatic.",
    COPT_INTPARAM_ROOTCUTROUNDS, -1, INT_MAX);

  AddSolverOption("mip:nodecutrounds nodecutrounds",
    "Rounds of cutting-planes generation of search tree node;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_NODECUTROUNDS, -1, INT_MAX);

  AddSolverOption("mip:heurlevel heurlevel",
    "Level of heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_HEURLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:roundingheurlevel roundingheurlevel",
    "Level of rounding heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_ROUNDINGHEURLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:divingheurlevel divingheurlevel",
    "Level of diving heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DIVINGHEURLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:submipheurlevel submipheurlevel",
    "Level of Sub-MIP heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_SUBMIPHEURLEVEL,
    alg_values_level, -1);


  AddSolverOption("mip:strongbranching strongbranching",
    "Level of strong branching:\n"
    "\n.. value-table::\n", COPT_INTPARAM_STRONGBRANCHING,
    alg_values_level, -1); 

    AddSolverOption("mip:conflictanalysis conflictanalysis",
      "Whether to perform conflict analysis:\n"
      "\n.. value-table::\n", COPT_INTPARAM_CONFLICTANALYSIS,
      values_autonoyes_, -1);

  AddSolverOption("mip:gap mipgap",
      "Relative optimality gap, default 1e-4.\n",
        COPT_DBLPARAM_RELGAP, 0.0, DBL_MAX);

  AddSolverOption("pre:dualize dualize",
    "Whether to dualize the problem before solving it:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DUALIZE,
    values_autonoyes_, -1);

  AddSolverOption("pre:solve presolve",
    "Level of presolving to perform before solving the problem:\n"
    "\n.. value-table::\n", COPT_INTPARAM_PRESOLVE, 
    alg_values_level, -1);

  AddSolverOption("pre:scale scale",
    "Whether to scale the problem:\n"
    "\n.. value-table::\n"
    "Scaling typically reduces solution times, but it may lead "
    "to larger constraint violations in the original, unscaled "
    "model. Choosing a different scaling option can sometimes "
    "improve performance for particularly numerically difficult "
    "models.",
    COPT_INTPARAM_SCALING, values_autonoyes_, -1);

  AddSolverOption("tech:threads threads",
    "Number of threads to use;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_THREADS, -1, 128);

  AddSolverOption("tech:barrierthreads barthreads",
      "Number of threads used by the barrier algorithm;\n"
      "default -1 ==> use value in tech:threads.",
    COPT_INTPARAM_BARTHREADS, -1, 128);

  AddSolverOption("tech:crossoverthreads crossoverthreads",
    "Number of threads used by crossover;\n"
    "default -1 ==> use value in tech:threads.",
    COPT_INTPARAM_CROSSOVERTHREADS, -1, 128);

  AddSolverOption("tech:simplexthreads simplexthreads",
    "Number of threads used by dual simplex;\n"
    "default -1 ==> use value in tech:threads.",
    COPT_INTPARAM_SIMPLEXTHREADS, -1, 128);


  AddSolverOption("tech:miptasks miptasks",
    "Number of MIP tasks in parallel;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_MIPTASKS, -1, 255);

  AddSolverOption("lim:time timelim timelimit",
      "limit on solve time (in seconds; default: no limit).",
      COPT_DBLPARAM_TIMELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:soltime soltimelim soltimelimit",
    "Time limit if a primal feasible solution has been found "
    "(in seconds; default: no limit).",
    COPT_DBLPARAM_SOLTIMELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:nodes nodelim nodelimit",
    "Maximum MIP nodes to explore (default: no limit).",
    COPT_INTPARAM_NODELIMIT, 0, INT_MAX);

  AddSolverOption("lp:method method lpmethod",
    "Which algorithm to use for non-MIP problems:\n"
    "\n.. value-table::\n", COPT_INTPARAM_LPMETHOD, 
    lp_values_method, -1);

  AddSolverOption("alg:feastol feastol",
    "Primal feasibility tolerance (default 1e-6).",
    COPT_DBLPARAM_FEASTOL, 1e-9, 1e-4);

  AddSolverOption("alg:dualfeastol dualfeastol",
    "Tolerance for dual solutions and reduced cost (default 1e-6).",
    COPT_DBLPARAM_DUALTOL, 1e-6, 1e-4);

  AddSolverOption("alg:matrixtol matrixtol",
    "nput matrix coefficient tolerance (default 1e-10).",
    COPT_DBLPARAM_MATRIXTOL, 0.0, 1e-7);

  AddSolverOption("alg:iismethod iismethod",
    "Which method to use when finding an IIS (irreducible infeasible "
    "set of constraints, including variable bounds):\n"
    "\n.. value-table::\n",
    COPT_INTPARAM_IISMETHOD, values_iismethod, -1);

  AddSolverOption("bar:crossover crossover",
    "Execute crossover to transform a barrier solution to a basic one (default 1).",
    COPT_INTPARAM_CROSSOVER, values_crossover, 1);

  ReplaceOptionDescription("alg:feasrelax",
    "Whether to modify the problem into a feasibility "
    "relaxation problem:\n"
    "\n"
    "| 0 = No (default)\n"
    "| 1 = Yes, minimizing the weighted sum of violations\n"
    "| 2 = Yes, minimizing the weighted sum of squared violations\n"
    "| 3 = Yes, minimizing the weighted count of violations\n"
    "| 4-6 = Same objective as 1-3, but also optimize the "
    "original objective, subject to the violation "
    "objective being minimized.\n"
    "\n"
    "Weights are given by suffixes .lbpen and .ubpen on variables "
    "and .rhspen on constraints (when nonnegative), else by keywords "
    "alg:lbpen, alg:ubpen, and alg:rhspen, respectively (default values = 1). "
    "Weights < 0 are treated as Infinity, allowing no violation.");

    /*
          AddStoredOption("alg:feasrelax feasrelax",
                      "Whether to modify the problem into a feasibility "
                          "relaxation problem:\n"
                          "\n"
                          "| 0 = No (default)\n"
                          "| 1 = Yes, minimizing the weighted sum of violations\n"
                          "| 2 = Yes, minimizing the weighted sum of squared violations\n"
                          "| 3 = Yes, minimizing the weighted count of violations\n"
                          "| 4-6 = Same objective as 1-3, but also optimize the "
                             "original objective, subject to the violation "
                             "objective being minimized.\n"
                      "\n"
                      "Weights are given by suffixes .lbpen and .ubpen on variables "
                      "and .rhspen on constraints (when nonnegative), else by keywords "
                      "alg:lbpen, alg:ubpen, and alg:rhspen, respectively (default values = 1). "
                      "Weights < 0 are treated as Infinity, allowing no violation.",
          feasrelax().mode_);
      AddStoredOption("alg:lbpen lbpen", "See alg:feasrelax.",
          storedOptions_.lbpen_);
      AddStoredOption("alg:ubpen ubpen", "See alg:feasrelax.",
          storedOptions_.ubpen_);
      AddStoredOption("alg:rhspen rhspen", "See alg:feasrelax.",
          storedOptions_.rhspen_);
          */
}


double CoptBackend::MIPGap() {
  // TODO Check what happens if not solved
  return getDblAttr(COPT_DBLATTR_BESTGAP);
}
double CoptBackend::BestDualBound() {
  // TODO Check what happens if not solved
  return getDblAttr(COPT_DBLATTR_BESTBND);
}

double CoptBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> CoptBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  COPT_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case COPT_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case COPT_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case COPT_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case COPT_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case COPT_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Copt VBasis value: {}", s));
    }
  }
  return vars;
}

ArrayRef<int> CoptBackend::ConStatii() {
  std::vector<int> cons(NumLinCons());
  COPT_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case COPT_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case COPT_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case COPT_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case COPT_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case COPT_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Copt VBasis value: {}", s));
    }
  }
  return cons;
}

void CoptBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = COPT_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = COPT_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = COPT_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = COPT_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = COPT_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!COPT_GetColInfo(lp(), COPT_DBLINFO_LB, 1, index, &lb) && 
        !COPT_GetColInfo(lp(), COPT_DBLINFO_UB, 1, index, &ub))
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
  COPT_SetBasis(lp(), stt.data(), NULL);
}

void CoptBackend::ConStatii(ArrayRef<int> cst) {
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = COPT_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = COPT_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  COPT_SetBasis(lp(), NULL, stt.data());
}

SolutionBasis CoptBackend::GetBasis() {
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

void CoptBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void CoptBackend::ComputeIIS() {
  COPT_CCALL(COPT_ComputeIIS(lp()));
  SetStatus(ConvertCOPTStatus());   // could be new information
}

IIS CoptBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}


std::vector<int> getIIS(copt_prob* prob, int num, 
    int(*getlb)(copt_prob*, int, const int*, int*),
  int(*getub)(copt_prob*, int, const int*, int*)) {
  std::vector<int> iis_lb(num);
  std::vector<int> iis_ub(num);
  COPT_CCALL(getlb(prob, num, NULL, iis_lb.data()));
  COPT_CCALL(getub(prob, num, NULL, iis_ub.data()));
  for (size_t i = iis_lb.size(); i--; ) {
    if (iis_ub[i]) {
      if (iis_lb[i])
        iis_lb[i] = (int)IISStatus::fix;
      else
        iis_lb[i] = (int)IISStatus::upp;
    }
    else {
      if (iis_lb[i])
        iis_lb[i] = (int)IISStatus::low;
      else
        iis_lb[i] = (int)IISStatus::non;
    }
  }
  return iis_lb;
}

static void ConvertIIS2AMPL(std::vector<int>& ai) {
  for (int i = ai.size(); i--; ) {
    ai[i] = int(ai[i] ? IISStatus::mem : IISStatus::non);
  }
}


ArrayRef<int> CoptBackend::VarsIIS() {
  return getIIS(lp(), NumVars(), COPT_GetColLowerIIS, COPT_GetColUpperIIS);
}
pre::ValueMapInt CoptBackend::ConsIIS() {
  auto iis_lincon =  getIIS(lp(), NumLinCons(), COPT_GetRowLowerIIS, COPT_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  COPT_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  COPT_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };

}


ArrayRef<double> CoptBackend::Ray() {
  auto uray_pres = getVarInfo(COPT_DBLINFO_PRIMALRAY);
  auto mv = GetValuePresolver().PostsolveSolution({ uray_pres });
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> CoptBackend::DRay() {
  auto fd = getConInfo(COPT_DBLINFO_DUALFARKAS);
  auto vm = GetValuePresolver().PostsolveSolution({
                                               {},
                                               {{{CG_Linear, std::move(fd)}}}
    });
  return vm.GetConValues().MoveOut();        // need the vector itself
}


void CoptBackend::AddPrimalDualStart(Solution sol0_unpres) {
  if (IsMIP())
    return;
  auto mv = GetValuePresolver().PresolveSolution(
    { sol0_unpres.primal, sol0_unpres.dual });
  auto x0 = mv.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);
  COPT_CCALL(COPT_SetLpSolution(lp(), x0.data(), nullptr, pi0.data(), nullptr));
}

void CoptBackend::AddMIPStart(
		ArrayRef<double> x0_unpres, ArrayRef<int> sparsity_unpres) {
  auto mv = GetValuePresolver().PresolveSolution({ x0_unpres });
  auto ms = GetValuePresolver().PresolveGenericInt({ sparsity_unpres });
  auto x0 = mv.GetVarValues()();
  auto s0 = ms.GetVarValues()();
  std::vector<int> idx;                 // Create sparse vector
  idx.reserve(x0.size());
  std::vector<double> val;
  val.reserve(x0.size());
  for (int i = 0; i < (int)x0.size(); ++i) {
    if (s0[i]) {
      idx.push_back(i);
      val.push_back(x0[i]);
    }
  }
  COPT_CCALL(COPT_AddMipStart(lp(), val.size(), idx.data(), const_cast<double*>(val.data())));
}


} // namespace mp

  // AMPLs
AMPLS_MP_Solver* Open_copt(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::CoptBackend()}, 
    cb);
}

void AMPLSClose_copt(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_copt(AMPLS_MP_Solver* slv) {
  return dynamic_cast<mp::CoptBackend*>(AMPLSGetBackend(slv))->lp();
}

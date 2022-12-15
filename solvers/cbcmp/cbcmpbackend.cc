#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cbcmpbackend.h"

extern "C" {
  #include "cbcmp-ampls-c-api.h"    // Cbcmp AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptCbcmp(void* prob) {
  //return CBCMP_Interrupt((cbcmp_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCbcmpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CbcmpBackend()};
}


namespace mp {

/// Create Cbcmp Model Manager
/// @param gc: the Cbcmp common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return CbcmpModelMgr
std::unique_ptr<BasicModelManager>
CreateCbcmpModelMgr(CbcmpCommon&, Env&, pre::BasicValuePresolver*&);


CbcmpBackend::CbcmpBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateCbcmpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

CbcmpBackend::~CbcmpBackend() {
  CloseSolver();
}

void CbcmpBackend::OpenSolver() {
  set_lp(Cbc_newModel()); // Assign it
}

void CbcmpBackend::CloseSolver() {
  if ( lp() != NULL ) 
    Cbc_deleteModel(lp());
}

const char* CbcmpBackend::GetBackendName()
  { return "CbcmpBackend"; }

std::string CbcmpBackend::GetSolverVersion() {
  return Cbc_getVersion();
}


bool CbcmpBackend::IsMIP() const {
  return Cbc_getNumIntegers(lp())> 0;
}

bool CbcmpBackend::IsQCP() const {
  return false;
}

ArrayRef<double> CbcmpBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  auto sol = Cbc_getColSolution(lp());
  for (int i = 0; i < num_vars; ++i)
    x[i] = sol[i];
  return x;
}

pre::ValueMapDbl CbcmpBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> CbcmpBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
 // int error = CBCMP_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 0;
  if (error)
    pi.clear();
  return pi;
}

double CbcmpBackend::ObjectiveValue() const {
  return Cbc_getObjValue(lp());
}

double CbcmpBackend::NodeCount() const {
    return Cbc_getNodeCount(lp());
}

double CbcmpBackend::SimplexIterations() const {
  // TODO which one is it?
  return Cbc_getIterationCount(lp());
}

int CbcmpBackend::BarrierIterations() const {
  // TODO which one is it?
  return Cbc_getIterationCount(lp());
}

inline bool ends_with(std::string const& value, std::string const& ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}
void CbcmpBackend::DoWriteProblem(const std::string & name) {
  if (ends_with(name, ".lp"))
    Cbc_writeLp(lp(), name.c_str());
  else if (ends_with(name, ".mps"))
    Cbc_writeMps(lp(), name.c_str());
  else
    throw std::runtime_error("Can only export '.lp' or '.mps' files.");
}


void CbcmpBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCbcmp, lp());
  // TODO Check interrupter
  //CBCMP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CbcmpBackend::Solve() {
  if (storedOptions_.timeLimit_ != 0)
    Cbc_setMaximumSeconds(lp(), storedOptions_.timeLimit_);
  CBCMP_CCALL(Cbc_solve(lp()));
  WindupCBCMPSolve();
}

void CbcmpBackend::WindupCBCMPSolve() { }

void CbcmpBackend::ReportResults() {
  ReportCBCMPResults();
  BaseBackend::ReportResults();
}

void CbcmpBackend::ReportCBCMPResults() {
  SetStatus( ConvertCBCMPStatus() );
  AddCBCMPMessages();
  if (need_multiple_solutions())
    ReportCBCMPPool();
}
std::vector<double> CbcmpBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // CBCMP_CCALL(CBCMP_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double CbcmpBackend::getPoolObjective(int i)
{
  double obj;
 // CBCMP_CCALL(CBCMP_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void CbcmpBackend::ReportCBCMPPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(CBCMP_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void CbcmpBackend::AddCBCMPMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> CbcmpBackend::ConvertCBCMPStatus() {
  namespace sol = mp::sol;
  if (Cbc_isProvenOptimal(lp()))
    return { sol::SOLVED, "optimal solution" };
  if (Cbc_isProvenInfeasible(lp()))
    return { sol::INFEASIBLE, "infeasible problem" };
  if (Cbc_isContinuousUnbounded(lp()))
    return { sol::UNBOUNDED, "unbounded problem" };

  switch(Cbc_status(lp())){
  case -1:
    return { sol::UNKNOWN, "unfinished" };
  case 1:
    return { sol::LIMIT, "Hit a limit" };
  case 2:
    return { sol::NUMERIC, "Numeric issues" };
  case 5:
    return { sol::INTERRUPTED, "Interrupted" };
  default:
    return { sol::UNKNOWN, "not solved" };
  }
}


void CbcmpBackend::FinishOptionParsing() {
  int v = Cbc_getLogLevel(lp());
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo lp_values_method[] = {
  { "-1", "Automatic (default)", -1},
  { "1", "Dual simplex", 1},
  { "2", "Barrier", 2},
  { "3", "Crossover", 3},
  { "4", "Concurrent (simplex and barrier simultaneously)", 4}
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

static const mp::OptionValueInfo offon_values_[] = {
{"off", "Turn off", 0},
{"on", "Turn on", 0}
};

static const mp::OptionValueInfo biasLU_values_[] = {
{"UU", "", 0},
{"UX", "", 0},
{"LX", "", 0},
{"LL", "", 0}
};

static const mp::OptionValueInfo bscale_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"off1", "", 0},
{"on1", "", 0},
{"off2", "", 0},
{"on2", "", 0}
};

static const mp::OptionValueInfo cholesky_values_[] = {
{"native", "", 0},
{"dense", "", 0},
{"fudgeLong_dummy", "", 0},
{"wssmp_dummy", "", 0},
{"UniversityOfFlorida_dummy", "", 0},
{"Taucs_dummy", "", 0},
{"Mumps_dummy", "", 0},
{"Pardiso_dummy", "", 0}
};

static const mp::OptionValueInfo OffOnBothBefore_values_[] = {
{"off", "disabled", 0},
{"on", "use every node in the tree", 0},
{"both", "", 0},
{"before", "", 0}
};
static const mp::OptionValueInfo OffOnBothBeforeOften_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"both", "", 0},
{"before", "", 0},
{"often", "", 0}
};

static const mp::OptionValueInfo VndVariableNeighborhoodSearch_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"both", "", 0},
{"before", "", 0},
{"intree", "", 0}
};

static const mp::OptionValueInfo combineSolutions_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"both", "", 0},
{"before", "", 0},
{"onquick", "", 0},
{"bothquick", "", 0},
{"beforequick", "", 0}
};

static const mp::OptionValueInfo constraintfromCutoff_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"variable", "", 0},
{"forcevariable", "", 0},
{"conflict", "", 0}
};

static const mp::OptionValueInfo costStrategy_values_[] = {
{"off", "", 0},
{"priorities", "", 0},
{"columnOrder", "", 0},
{"01first", "", 0},
{"01last", "", 0},
{"length", "", 0},
{"singletons", "", 0},
{"nonzero", "", 0},
{"generalForce", "", 0}
};

static const mp::OptionValueInfo crossover_values_[] = {
{"on", "", 0},
{"off", "", 0},
{"maybe", "", 0},
{"presolve", "", 0}
};

static const mp::OptionValueInfo dualPivot_values_[] = {
{"automatic", "", 0},
{"dantzig", "", 0},
{"partial", "", 0},
{"steepest", "", 0},
{"PEsteepest", "", 0},
{"PEdantzig", "", 0}
};

static const mp::OptionValueInfo primalPivot_values_[] = {
{"auto!matic", "", 0},
{"exact", "", 0},
{"dantzig", "", 0},
{"partial", "", 0},
{"steepest", "", 0},
{"change", "", 0},
{"sprint", "", 0},
{"PEsteepest", "", 0},
{"PEdantzig", "", 0}
};


static const mp::OptionValueInfo factorization_values_[] = {
{"normal", "", 0},
{"dense", "", 0},
{"simple", "", 0},
{"osl", "", 0}
};

static const mp::OptionValueInfo gammadelta_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"gamma", "", 0},
{"delta", "", 0},
{"onstrong", "", 0},
{"gammastrong", "", 0},
{"deltastrong", "", 0}
};




static const mp::OptionValueInfo cutsOnOff_values_[] = {
{"off", "disabled", 0},
{"on", "enabled", 0},
{"root", "enabled only on root node", 0},
{"ifmove", "enabled in the tree if it moves the objective value", 0},
{"forceOn", "enabled at every node", 0}
};


static const mp::OptionValueInfo cutsToOnGlobal[] = {
{"off", "disabled", 0},
{"on", "enabled", 0},
{"root", "enabled only on root node", 0},
{"ifmove", "enabled in the tree if it moves the objective value", 0},
{"forceOn", "enabled at every node", 0},
{"onglobal", "", 0}
};

static const mp::OptionValueInfo reduce2AndSplitCuts_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"root", "", 0},
{"longOn", "", 0},
{"longRoot", "", 0}
};



static const mp::OptionValueInfo knapsackCuts_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"root", "", 0},
{"ifmove", "", 0},
{"forceOn", "", 0},
{"onglobal", "", 0},
{"forceandglobal", "", 0}
};
static const mp::OptionValueInfo gomoryCuts_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"root", "", 0},
{"ifmove", "", 0},
{"forceOn", "", 0},
{"onglobal", "", 0},
{"forceandglobal", "", 0},
{"forceLongOn", "", 0},
{"long", "", 0}
};

static const mp::OptionValueInfo GMICuts_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"root", "", 0},
{"ifmove", "", 0},
{"forceOn", "", 0},
{"endonly", "", 0},
{"long", "", 0},
{"longroot", "", 0},
{"longifmove", "", 0},
{"forceLongOn", "", 0},
{"longendonly", "", 0}
};
static const mp::OptionValueInfo twoMirCuts_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"root", "", 0},
{"ifmove", "", 0},
{"forceOn", "", 0},
{"onglobal", "", 0},
{"forceandglobal", "", 0},
{"forceLongOn", "", 0}
};

static const mp::OptionValueInfo latwomirCuts_values_[] = {
{"off", "", 0},
{"endonlyroot", "", 0},
{"endcleanroot", "", 0},
{"endbothroot", "", 0},
{"endonly", "", 0},
{"endclean", "", 0},
{"endboth", "", 0},
{"onlyaswell", "", 0},
{"cleanaswell", "", 0},
{"bothaswell", "", 0},
{"onlyinstead", "", 0},
{"cleaninstead", "", 0},
{"bothinstead", "", 0}
};


static const mp::OptionValueInfo lagomoryCuts_values_[] = {
{"off", "", 0},
{"endonlyroot", "", 0},
{"endcleanroot", "", 0},
{"root", "", 0},
{"endonly", "", 0},
{"endclean", "", 0},
{"endboth", "", 0},
{"onlyaswell", "", 0},
{"cleanaswell", "", 0},
{"bothaswell", "", 0},
{"onlyinstead", "", 0},
{"cleaninstead", "", 0},
{"bothinstead", "", 0},
{"onlyaswellroot", "", 0},
{"cleanaswellroot", "", 0},
{"bothaswellroot", "", 0}
};


static const mp::OptionValueInfo probingCuts_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"root", "", 0},
{"ifmove", "", 0},
{"forceOn", "", 0},
{"onglobal", "", 0},
{"forceonglobal", "", 0},
{"forceOnBut", "", 0},
{"forceOnStrong", "", 0},
{"forceOnButStrong", "", 0},
{"strongRoot", "", 0}
};

static const mp::OptionValueInfo nodeStrategy_values_[] = {
{"hybrid", "", 0},
{"fewest", "", 0},
{"depth", "", 0},
{"upfewest", "", 0},
{"downfewest", "", 0},
{"updepth", "", 0},
{"downdepth", "", 0}
};


static const mp::OptionValueInfo presolve_values_[] = {
{"on", "", 0},
{"off", "", 0},
{"more", "", 0},
{"file", "", 0}
};


static const mp::OptionValueInfo preprocess_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"save", "", 0},
{"equal", "", 0},
{"sos", "", 0},
{"trysos", "", 0},
{"equalall", "", 0},
{"strategy", "", 0},
{"aggregate", "", 0},
{"forcesos", "", 0},
{"stopaftersaving", "", 0}
};

static const mp::OptionValueInfo proximitySearch_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"both", "", 0},
{"before", "", 0},
{"10", "", 0},
{"100", "", 0},
{"300", "", 0}
};

static const mp::OptionValueInfo Rens_values_[] = {
{"off", "", 0},
{"on", "", 0},
{"both", "", 0},
{"before", "", 0},
{"200", "", 0},
{"1000", "", 0},
{"10000", "", 0},
{"dj", "", 0},
{"djbefore", "", 0},
{"usesolution", "", 0}
};



static const mp::OptionValueInfo scaling_values_[] = {
{"off", "", 0},
{"equilibrium", "", 0},
{"geometric", "", 0},
{"automatic", "", 0},
{"dynamic", "", 0},
{"rowsonly", "", 0}
};


void CbcmpBackend::InitCustomOptions() {

  // Uncomment to print the whole list of CBC parameters - from cbc
  // Useful when updating the list of params in the driver.
  // GetCBCParamsList();

  set_option_header(
      "CBCMP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cbcmp_options``. For example::\n"
      "\n"
      "  ampl: option cbcmp_options 'mipgap=1e-6';\n");


  AddSolverOption("tech:outlev outlev",
    "0*/1: Whether to write log lines (chatter) to stdout and to file.",
    "log", 0, 4);

  _options.addOption("timelimit", Cbc_setMaximumSeconds);

  std::string name;
  std::string desc;
  mp::OptionValueInfo k;
  for (auto p : lp()->cbcData->parameters_)
  {
    //   p.printString();
    //fmt::print(p.name());

    if (p.type() >= CLP_PARAM_DBL_PRIMALTOLERANCE &&
      p.type() <= CBC_PARAM_DBL_DEXTRA5) {
      if (p.type() == CBC_PARAM_DBL_DEXTRA5)
        name = "double:dextra5 dextra5";
      else
        name = fmt::format("double:{} {}", p.name(), p.name());
      desc = fmt::format("{} (default {}).", p.shortHelp(), p.doubleValue());
      AddSolverOption(name.c_str(), desc.c_str(), p.name().c_str(),
        p.lowerDoubleValue(), p.upperDoubleValue());
    }
    else if (p.type() >= CLP_PARAM_INT_SOLVERLOGLEVEL &&
      p.type() <= CBC_PARAM_INT_MOREMOREMIPOPTIONS)
    {
      name = fmt::format("int:{} {}", p.name(), p.name());
      desc = fmt::format("{} (default {}).", p.shortHelp(), p.intValue());
      AddSolverOption(name.c_str(), desc.c_str(), p.name().c_str(),
        p.lowerIntValue(), p.upperIntValue());
    }
    // String options are dealt with separately
  }

  AddStoredOption("lim:time timelim timelimit",
    "Limit on solve time (in seconds; default: no limit).",
    storedOptions_.timeLimit_);
  AddSolverOption("pre:autoScale",
    "Whether to scale objective, rhs and bounds of problem if they look odd:\n"
    "\n.. value-table::\n",
    "autoScale", offon_values_, "off");

   AddSolverOption("pre:scaling scaling",
    "Whether to scale problem\n"
    "\n.. value-table::\n",
    "scaling", scaling_values_, "NULL");

  AddSolverOption("pre:biasLU",
    "Whether factorization biased towards U:\n"
    "\n.. value-table::\n",
    "biasLU", biasLU_values_, "UX");


  // ******************** HEURISTICS ********************
  AddSolverOption("mip:combineSolutions combineSolutions",
    "Whether to use combine solution heuristic\n"
    "\n.. value-table::\n",
    "combineSolutions", combineSolutions_values_, "NULL");

  AddSolverOption("mip:combine2Solutions combine2Solutions",
    "Whether to use crossover solution heuristic\n"
    "\n.. value-table::\n",
    "combine2Solutions", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:dins Dins",
    "Whether to try Distance Induced Neighborhood Search\n"
    "\n.. value-table::\n",
    "Dins", OffOnBothBeforeOften_values_, "NULL");

  AddSolverOption("mip:divingsome DivingSome",
    "Whether to try Diving heuristics\n"
    "\n.. value-table::\n",
    "DivingSome", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:divingcoefficient DivingCoefficient",
    "Whether to try Coefficient diving heuristic\n"
    "\n.. value-table::\n",
    "DivingCoefficient", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:divingfractional DivingFractional",
    "Whether to try Fractional diving heuristic\n"
    "\n.. value-table::\n",
    "DivingFractional", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:divingguided DivingGuided",
    "Whether to try Guided diving heuristic\n"
    "\n.. value-table::\n",
    "DivingGuided", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:divinglinesearch DivingLineSearch",
    "Whether to try Linesearch diving heuristic\n"
    "\n.. value-table::\n",
    "DivingLineSearch", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:divingpseudocost DivingPseudoCost",
    "Whether to try Pseudocost diving heuristic\n"
    "\n.. value-table::\n",
    "DivingPseudoCost", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:divingvectorlength DivingVectorLength",
    "Whether to try Vectorlength diving heuristic\n"
    "\n.. value-table::\n",
    "DivingVectorLength", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:dwHeuristic dwHeuristic",
    "Whether to try Dantzig Wolfe heuristic\n"
    "\n.. value-table::\n",
    "dwHeuristic", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:feasibilitypumpfeasibilityPump",
    "Whether to try the Feasibility Pump heuristic\n"
    "\n.. value-table::\n",
    "feasibilityPump", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:greedyheuristic greedyHeuristic",
    "Whether to use a greedy heuristic\n"
    "\n.. value-table::\n",
    "greedyHeuristic", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:naiveheuristics naiveHeuristics",
    "Whether to try some stupid heuristic\n"
    "\n.. value-table::\n",
    "naiveHeuristics", OffOnBothBefore_values_, "NULL");


  AddSolverOption("mip:pivotandcomplement pivotAndComplement",
    "Whether to try Pivot and Complement heuristic\n"
    "\n.. value-table::\n",
    "pivotAndComplement", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:pivotandfix pivotAndFix",
    "Whether to try Pivot and Fix heuristic\n"
    "\n.. value-table::\n",
    "pivotAndFix", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:randomizedrounding randomizedRounding",
    "Whether to try randomized rounding heuristic\n"
    "\n.. value-table::\n",
    "randomizedRounding", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:rens Rens",
    "Whether to try Relaxation Enforced Neighborhood Search\n"
    "\n.. value-table::\n",
    "Rens", Rens_values_, "NULL");


  AddSolverOption("mip:rins Rins",
    "Whether to try Relaxed Induced Neighborhood Search\n"
    "\n.. value-table::\n",
    "Rins", OffOnBothBeforeOften_values_, "NULL");

  AddSolverOption("mip:roundingheuristic roundingHeuristic",
    "Whether to use simple (but effective) Rounding heuristic\n"
    "\n.. value-table::\n",
    "roundingHeuristic", OffOnBothBefore_values_, "NULL");

  AddSolverOption("mip:vndvariableneighborhoodsearch VndVariableNeighborhoodSearch",
    "Whether to try Variable Neighborhood Search\n"
    "\n.. value-table::\n",
    "VndVariableNeighborhoodSearch", VndVariableNeighborhoodSearch_values_, "NULL");

  AddSolverOption("mip:heuristics heuristicsOnOff",
    "Switches most primal heuristics on or off\n"
    "\n.. value-table::\n",
    "heuristicsOnOff", offon_values_, "NULL");

  // ******************** END HEURISTICS ********************


  AddSolverOption("bar:bscale bscale",
    "Whether to scale in barrier (and ordering speed)\n"
    "\n.. value-table::\n",
    "bscale", bscale_values_, "NULL");

  AddSolverOption("bar:cholesky cholesky",
    "Which cholesky algorithm\n"
    "\n.. value-table::\n",
    "cholesky", cholesky_values_, "NULL");



  AddSolverOption("mip:constraintfromCutoff constraintfromCutoff",
    "Whether to use cutoff as constraint\n"
    "\n.. value-table::\n",
    "constraintfromCutoff", constraintfromCutoff_values_, "NULL");

  AddSolverOption("mip:costStrategy costStrategy",
    "How to use costs for branching priorities\n"
    "\n.. value-table::\n",
    "costStrategy", costStrategy_values_, "NULL");

  AddSolverOption(":cplexUse cplexUse",
    "Whether to use Cplex!\n"
    "\n.. value-table::\n",
    "cplexUse", offon_values_, "NULL");

  AddSolverOption("bar:crash crash",
    "Whether to create basis for problem\n"
    "\n.. value-table::\n",
    "crash", offon_values_, "NULL");

  AddSolverOption("bar:crossover crossover",
    "Whether to get a basic solution with the simplex algorithm after the barrier algorithm finished\n"
    "\n.. value-table::\n",
    "crossover", crossover_values_, "NULL");

  AddSolverOption("cut:cut cutsOnOff",
    "Switches all cut generators on or off\n"
    "\n.. value-table::\n",
    "cutsOnOff", cutsOnOff_values_, "NULL");

  AddSolverOption("cut:probingCuts probingCuts",
    "Whether to use Probing cuts\n"
    "\n.. value-table::\n",
    "probingCuts", probingCuts_values_, "NULL");


  AddSolverOption("lp:dualpivot dualPivot",
    "Dual pivot choice algorithm\n"
    "\n.. value-table::\n",
    "dualPivot", dualPivot_values_, "NULL");
  
  AddSolverOption("lp:primalpivot primalPivot",
    "Primal pivot choice algorithm\n"
    "\n.. value-table::\n",
    "primalPivot", primalPivot_values_, "NULL");


  AddSolverOption("pre:factorization factorization",
    "Which factorization to use\n"
    "\n.. value-table::\n",
    "factorization", factorization_values_, "NULL");

  AddSolverOption("pre:sparsefactor sparseFactor",
    "Whether factorization treated as sparse\n"
    "\n.. value-table::\n",
    "sparseFactor", offon_values_, "NULL");




  // ******************** CUTS ********************
  AddSolverOption("cut:cliqueCuts cliqueCuts",
    "Whether to use Clique cuts\n"
    "\n.. value-table::\n",
    "cliqueCuts", cutsToOnGlobal, "NULL");

  AddSolverOption("cut:flowcovercuts flowCoverCuts",
    "Whether to use Flow Cover cuts\n"
    "\n.. value-table::\n",
    "flowCoverCuts", cutsToOnGlobal, "NULL");

  AddSolverOption("cut:gmicuts GMICuts",
    "Whether to use alternative Gomory cuts\n"
    "\n.. value-table::\n",
    "GMICuts", GMICuts_values_, "NULL");

  AddSolverOption("cut:gomorycuts gomoryCuts",
    "Whether to use Gomory cuts\n"
    "\n.. value-table::\n",
    "gomoryCuts", gomoryCuts_values_, "NULL");

  AddSolverOption("cut:knapsackcuts knapsackCuts",
    "Whether to use Knapsack cuts\n"
    "\n.. value-table::\n",
    "knapsackCuts", knapsackCuts_values_, "NULL");

  AddSolverOption("cut:lagomorycuts lagomoryCuts",
    "Whether to use Lagrangean Gomory cuts\n"
    "\n.. value-table::\n",
    "lagomoryCuts", lagomoryCuts_values_, "NULL");

  AddSolverOption("cut:latwomircuts latwomirCuts",
    "Whether to use Lagrangean TwoMir cuts\n"
    "\n.. value-table::\n",
    "latwomirCuts", latwomirCuts_values_, "NULL");

  AddSolverOption("cut:liftandprojectcuts liftAndProjectCuts",
    "Whether to use Lift and Project cuts\n"
    "\n.. value-table::\n",
    "liftAndProjectCuts", cutsOnOff_values_, "NULL");

  AddSolverOption(":mixedIntegerRoundingCuts mixedIntegerRoundingCuts",
    "Whether to use Mixed Integer Rounding cuts\n"
    "\n.. value-table::\n",
    "mixedIntegerRoundingCuts", cutsToOnGlobal, "NULL");

  AddSolverOption("cut:reduceandsplitcuts reduceAndSplitCuts",
    "Whether to use Reduce-and-Split cuts\n"
    "\n.. value-table::\n",
    "reduceAndSplitCuts", cutsOnOff_values_, "NULL");

  AddSolverOption("cut:reduce2andsplitcts reduce2AndSplitCuts",
    "Whether to use Reduce-and-Split cuts - style 2\n"
    "\n.. value-table::\n",
    "reduce2AndSplitCuts", reduce2AndSplitCuts_values_, "NULL");

  AddSolverOption("cut:residualcapacitycuts residualCapacityCuts",
    "Whether to use Residual Capacity cuts\n"
    "\n.. value-table::\n",
    "residualCapacityCuts", cutsOnOff_values_, "NULL");

  AddSolverOption("cut:twomircuts twoMirCuts",
    "Whether to use Two phase Mixed Integer Rounding cuts\n"
    "\n.. value-table::\n",
    "twoMirCuts", twoMirCuts_values_, "NULL");

  AddSolverOption("cut:zeroHalfCuts zeroHalfCuts",
    "Whether to use zero half cuts\n"
    "\n.. value-table::\n",
    "zeroHalfCuts", cutsToOnGlobal, "NULL");

  AddSolverOption("bar:gammadelta gamma(Delta)",
    "Whether to regularize barrier\n"
    "\n.. value-table::\n",
    "gamma(Delta)", gammadelta_values_, "NULL");






  AddSolverOption("bar:kkt KKT",
    "Whether to use KKT factorization in barrier\n"
    "\n.. value-table::\n",
    "KKT", offon_values_, "NULL");


  AddSolverOption("mip:localtreesearch localTreeSearch",
    "Whether to use local tree search when a solution is found\n"
    "\n.. value-table::\n",
    "localTreeSearch", offon_values_, "NULL");


  AddSolverOption("tech:messages messages",
    "Controls if Clpnnnn is printed\n"
    "\n.. value-table::\n",
    "messages", offon_values_, "NULL");



  AddSolverOption("mip:nodestrategy nodeStrategy",
    "What strategy to use to select the next node from the branch and cut tree\n"
    "\n.. value-table::\n",
    "nodeStrategy", nodeStrategy_values_, "NULL");


  AddSolverOption("alg:perturbation perturbation",
    "Whether to perturb the problem\n"
    "\n.. value-table::\n",
    "perturbation", offon_values_, "NULL");

  AddSolverOption("lp:pfi PFI",
    "Whether to use Product Form of Inverse in simplex\n"
    "\n.. value-table::\n",
    "PFI", offon_values_, "NULL");


  AddSolverOption("pre:presolve presolve",
    "Whether to presolve problem\n"
    "\n.. value-table::\n",
    "presolve", presolve_values_, "NULL");


  AddSolverOption("pre:preprocess preprocess",
    "Whether to use integer preprocessing\n"
    "\n.. value-table::\n",
    "preprocess", preprocess_values_, "NULL");
  

  AddSolverOption("mip:proximitysearch proximitySearch",
    "Whether to do proximity search heuristic\n"
    "\n.. value-table::\n",
    "proximitySearch", proximitySearch_values_, "NULL");
 }


double CbcmpBackend::MIPGap() {
  if (!(ObjectiveValue() == 0))
    return MIPGapAbs() / std::fabs(ObjectiveValue());
  return Infinity();
}
double CbcmpBackend::BestDualBound() {
  if (!IsMIP())
    return ObjectiveValue();
  // Todo: Check for objective constant
  return Cbc_getBestPossibleObjValue(lp());
}

double CbcmpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}

void OsiToMP(std::vector<int>& b) {
  int map[] =
  {
    (int)BasicStatus::sup,
    (int)BasicStatus::bas,
    (int)BasicStatus::low,
    (int)BasicStatus::upp
  };
  for (auto& s : b) 
    s = map[s];
}
SolutionBasis CbcmpBackend::GetBasis() {
  std::vector<int> varstt(NumVars());
  std::vector<int> constt(NumLinCons());
  lp()->solver_->getBasisStatus(varstt.data(), constt.data());
  OsiToMP(varstt);
  OsiToMP(constt);
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

void MPToOsi(std::vector<int>& b) {
  int map[] =
  {
    0,  // none
    1,  // bas
    0,  // sup
    2,  // low
    3,  // upp
    3, // equ
    2  // btw
  };
  for (auto& s : b)
    s = map[s];
}
void CbcmpBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  MPToOsi(varstt);
  MPToOsi(constt);
  lp()->solver_->setBasisStatus(varstt.data(), constt.data());
}


void CbcmpBackend::AddMIPStart(ArrayRef<double> x0) {
  //CBCMP_CCALL(CBCMP_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));

}


} // namespace mp

  // AMPLs
AMPLS_MP_Solver* AMPLSOpenCbcmp(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::CbcmpBackend()},
    slv_opt, cb);
}

void AMPLSCloseCbcmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void*  GetCbcmpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CbcmpBackend*>(AMPLSGetBackend(slv))->lp();
}

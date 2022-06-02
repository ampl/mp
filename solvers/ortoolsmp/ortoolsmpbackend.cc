#include <vector>
#include <climits>
#include <cfloat>
#include <iostream> // for file output
#include <fstream>
#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "ortoolsmpbackend.h"

extern "C" {
  #include "ortoolsmp-ampls-c-api.h"    // Ortools AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptOrtools(void* prob) { 
  return static_cast < operations_research::MPSolver*>(prob)->InterruptSolve();
}
}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateOrtoolsBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::OrtoolsBackend()};
}


namespace mp {

/// Create Ortools Model Manager
/// @param gc: the Ortools common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return OrtoolsModelMgr
std::unique_ptr<BasicModelManager>
CreateOrtoolsModelMgr(OrtoolsCommon&, Env&, pre::BasicValuePresolver*&);


OrtoolsBackend::OrtoolsBackend() : 
  status_(operations_research::MPSolver::NOT_SOLVED) {
  onlyRecordOptions_ = true;
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateOrtoolsModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

OrtoolsBackend::~OrtoolsBackend() {
  CloseSolver();
}

  void OrtoolsBackend::OpenSolver() {
  auto lp = operations_research::MPSolver::CreateSolver(storedOptions_.solver_);
  if (!lp)
    throw std::runtime_error(fmt::format(
      "Failed to create solver {}", storedOptions_.solver_));
  set_lp(lp); // Assign it to the underling "common" object
}

void OrtoolsBackend::CloseSolver() {
  delete lp();
}

const char* OrtoolsBackend::GetBackendName()
  { return "OrtoolsBackend"; }

std::string OrtoolsBackend::GetSolverVersion() {
  return "20220601"; // TODO: version is initialized in backend 
  // before the solver is instantiated - and actually before the 
  // specific solver is known.
  //return lp()->SolverVersion();
}


bool OrtoolsBackend::IsMIP() const {
  return lp()->IsMIP();
}

bool OrtoolsBackend::IsQCP() const {
  return false;
}

Solution OrtoolsBackend::GetSolution() {
  auto mv = GetValuePresolver().PostsolveSolution(
        { PrimalSolution(), DualSolution() } );
  return { mv.GetVarValues()(), mv.GetConValues()(),
    GetObjectiveValues() };   
}

ArrayRef<double> OrtoolsBackend::PrimalSolution() {
  std::vector<double> x(NumVars());
  for (std::size_t i = 0; i < NumVars(); i++)
    x[i] = lp()->variable(i)->solution_value();
  return x;
}

pre::ValueMapDbl OrtoolsBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> OrtoolsBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  for (std::size_t i = 0; i < num_cons; i++)
    pi[i] = lp()->constraint(i)->dual_value();
  return pi;
}

double OrtoolsBackend::ObjectiveValue() const {
  return lp()->Objective().Value();
}

double OrtoolsBackend::NodeCount() const { 
  return lp()->nodes(); }

double OrtoolsBackend::SimplexIterations() const { 
  return lp()->iterations(); }

int OrtoolsBackend::BarrierIterations() const {
  // TODO Check if i can access this through the MPSolver interface
  return 0;
}

bool endsWith(std::string str, std::string const& suffix) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  if (str.length() < suffix.length()) {
    return false;
  }
  return str.compare(str.length() - suffix.length(), suffix.length(), suffix) ==
         0;
}

void OrtoolsBackend::ExportModel(const std::string &file) {
  std::string output;
  if (endsWith(file, "mps"))
    lp()->ExportModelAsMpsFormat(false, false, &output);
  else if (endsWith(file, "lp"))
    lp()->ExportModelAsLpFormat(false, &output);
  else
    throw std::runtime_error("Output file should be either .mps or .lp");
  // Write file
  std::ofstream outputFile;
  outputFile.open(file);
  outputFile << output;
}


void OrtoolsBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptOrtools, lp());
}

void OrtoolsBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  status_ = lp()->Solve(params_);
  WindupORTOOLSSolve();
}

void OrtoolsBackend::WindupORTOOLSSolve() { }

void OrtoolsBackend::ReportResults() {
  ReportORTOOLSResults();
  BaseBackend::ReportResults();
}

void OrtoolsBackend::ReportORTOOLSResults() {
  SetStatus( ConvertORTOOLSStatus() );
  AddORTOOLSMessages();
  if (need_multiple_solutions())
    ReportORTOOLSPool();
}
void OrtoolsBackend::ReportORTOOLSPool() {
  if (!IsMIP())
    return;
  do {
    ReportIntermediateSolution({PrimalSolution(), {}, {ObjectiveValue()}});
  } while (lp()->NextSolution());
}


void OrtoolsBackend::AddORTOOLSMessages() {
  if (auto ni = SimplexIterations())
    AddToSolverMessage(
          fmt::format("{} simplex iterations\n", ni));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> OrtoolsBackend::ConvertORTOOLSStatus() {
  namespace sol = mp::sol;
  namespace ort = operations_research;
  switch (status_)
  {
  case ort::MPSolver::OPTIMAL:
    return { sol::SOLVED, "optimal solution" };
  case ort::MPSolver::FEASIBLE:
    return { sol::INTERRUPTED, "interrupted" };
  case ort::MPSolver::INFEASIBLE:
    return { sol::INFEASIBLE, "infeasible problem" };
      /// proven unbounded.
  case ort::MPSolver::UNBOUNDED:
    return { sol::UNBOUNDED, "unbounded problem" };
  case ort::MPSolver::ABNORMAL:
  case ort::MPSolver::MODEL_INVALID:
  case ort::MPSolver::NOT_SOLVED:
    return { sol::FAILURE, "failed" };
  default:
    return { sol::UNKNOWN, "not solved" };
  }
}


void OrtoolsBackend::FinishOptionParsing() {
  OpenSolver();
  ReplaySolverOptions();
  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
  set_verbose_mode(storedOptions_.outlev_ > 0);
  // Set stored options
  if (storedOptions_.outlev_ > 0) lp()->EnableOutput();
  if (storedOptions_.timelimit_ > 0)
    lp()->SetTimeLimit(absl::Seconds(storedOptions_.timelimit_));
  if (storedOptions_.threads_ > 0) lp()->SetNumThreads(storedOptions_.threads_);
}

// LP Basis

SolutionBasis OrtoolsBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
        {std::move(varstt), {{{CG_Linear, std::move(constt)}}}});
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
    assert(constt.size());
  }
  return {std::move(varstt), std::move(constt)};
}


orr::MPSolver::BasisStatus invbasismap[] = {
    orr::MPSolver::FREE,            // what to do
    orr::MPSolver::BASIC,           // basic
    orr::MPSolver::FREE,            // superbasic
    orr::MPSolver::AT_LOWER_BOUND,  // low
    orr::MPSolver::AT_UPPER_BOUND,  // upp
    orr::MPSolver::FIXED_VALUE,     //. equ
    orr::MPSolver::AT_UPPER_BOUND   // btw
};

int basismap[] = {(int)BasicStatus::sup, (int)BasicStatus::low,
                  (int)BasicStatus::upp, (int)BasicStatus::equ,
                  (int)BasicStatus::bas};

void OrtoolsBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis({basis.varstt, basis.constt});
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  std::vector<orr::MPSolver::BasisStatus> varstatus(varstt.size());
  std::vector<orr::MPSolver::BasisStatus> constatus(constt.size());
  for (auto i = 0; i < varstt.size(); i++) 
    varstatus[i] = invbasismap[varstt[i]];
  for (auto i = 0; i < constt.size(); i++) 
    constatus[i] = invbasismap[constt[i]];
  lp()->SetStartingLpBasis(varstatus, constatus);
}

ArrayRef<int> OrtoolsBackend::VarStatii() {
  std::vector<int> stt(NumVars());
  for (int i = 0; i < NumVars(); i++) {
    stt[i] = basismap[lp()->variable(i)->basis_status()];
  }
  return stt;
}

ArrayRef<int> OrtoolsBackend::ConStatii() {
  std::vector<int> stt(NumLinCons());
  for (int i = 0; i < NumLinCons(); i++) {
    stt[i] = basismap[lp()->constraint(i)->basis_status()];
  }
  return stt;
}




////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo values_method[] = {
  {"10", "Dual simplex", 0}, 
  {"11", "Primal simplex", 1}, 
  {"12", "Barrier", 2}};



void OrtoolsBackend::InitCustomOptions() {

  set_option_header(
      "ORTOOLS Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``ortools_options``. For example::\n"
      "\n"
      "  ampl: option ortools_options 'solver=scip';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  AddStoredOption("tech:solver solver",
      "Specifies the name of the solver to be used. Available options are: "
      "glop (default), cbc, scip",
      storedOptions_.solver_);
  
  AddStoredOption("tech:outlev outlev",
                  "0*/1: Whether to write log lines (chatter) to stdout",
                  storedOptions_.outlev_, 0, 1);

  AddStoredOption("tech:threads threads",
                    "How many threads to use when using the barrier algorithm "
                    "or solving MIP problems; default 0 ==> automatic choice.",
                    storedOptions_.threads_, 0,
                    std::numeric_limits<int>::max());

  AddStoredOption("lim:time timelim timelimit",
                    "Limit on solve time (in seconds; default: no limit).",
                      storedOptions_.timelimit_, 0.0,
                      std::numeric_limits<double>::max());
  
  AddSolverOption("mip:gap mipgap",
                  "Max. relative MIP optimality gap (default 1e-4).",
                  "RELATIVE_MIP_GAP", 1e-4, std::numeric_limits<double>::max());

  AddSolverOption("alg:feastol feastol",
                  "Primal feasibility tolerance (default 1e-6).",
                  "PRIMAL_TOLERANCE", 0.0, DBL_MAX);

  AddSolverOption("alg:dualfeastol dualfeastol",
                  "Dual feasibility tolerance (default 1e-6).",
                  "DUAL_TOLERANCE", 0.0, DBL_MAX);

  AddSolverOption("pre:solve presolve",
                  "Whether to use presolve:\n"
                  "\n.. value-table::\n",
                  "PRESOLVE", values_01_noyes_1default_, 0);

  AddSolverOption("alg:method method lpmethod simplex",
                  "Which algorithm to use for non-MIP problems or for the root "
                  "node of MIP problems:\n"
                  "\n.. value-table::\n",
                  "LP_ALGORITHM", values_method, 10);

  AddSolverOption(
        "pre:scale scale",
        "Whether to scale the problem:\n"
        "\n.. value-table::\n"
        "Scaling typically reduces solution times, but it may lead "
        "to larger constraint violations in the original, unscaled "
        "model. Choosing a different scaling option can sometimes "
        "improve performance for particularly numerically difficult "
        "models.", "SCALING", values_01_noyes_1default_, 1);
}

} // namespace mp


// AMPLs

AMPLS_MP_Solver* AMPLSOpenOrtools(
  const char* slv_opt) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::OrtoolsBackend()},
    slv_opt);
}

void AMPLSCloseOrtools(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetOrtoolsmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::OrtoolsBackend*>(AMPLSGetBackend(slv))->lp();
}

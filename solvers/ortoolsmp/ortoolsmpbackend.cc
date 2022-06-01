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
  int status = 0;
  //getenv("ortools_solver")
  auto lp = operations_research::MPSolver::CreateSolver(storedOptions_.solver_);
  if (!lp)
    throw std::runtime_error(fmt::format(
      "Failed to create solver {}", storedOptions_.solver_));
  set_lp(lp); // Assign it
}

void OrtoolsBackend::CloseSolver() {
  delete lp();
}

const char* OrtoolsBackend::GetBackendName()
  { return "OrtoolsBackend"; }

std::string OrtoolsBackend::GetSolverVersion() {
  return "20220601";
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
  status_ = lp()->Solve();
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
std::vector<double> OrtoolsBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // ORTOOLS_CCALL(ORTOOLS_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double OrtoolsBackend::getPoolObjective(int i)
{
  double obj;
 // ORTOOLS_CCALL(ORTOOLS_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void OrtoolsBackend::ReportORTOOLSPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(ORTOOLS_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
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


////////////////////////////// OPTIONS /////////////////////////////////


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

  AddStoredOption(
      "tech:solver solver",
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

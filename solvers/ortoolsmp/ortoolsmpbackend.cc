#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "ortoolsmpbackend.h"

extern "C" {
  #include "ortoolsmp-ampls-c-api.h"    // Ortools AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptOrtools(void* prob) {
  //return ORTOOLS_Interrupt((ortools_prob*)prob);
  return true;
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


OrtoolsBackend::OrtoolsBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateOrtoolsModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

OrtoolsBackend::~OrtoolsBackend() {
  CloseSolver();
}

void OrtoolsBackend::OpenSolver() {
  int status = 0;
  std::string solver = "SCIP";
  auto lp = operations_research::MPSolver::CreateSolver(solver);
  if (!lp)
    throw std::runtime_error(fmt::format(
      "Failed to create solver {}", solver));
  set_lp(lp); // Assign it
}

void OrtoolsBackend::CloseSolver() {
  delete lp();
}

const char* OrtoolsBackend::GetBackendName()
  { return "OrtoolsBackend"; }

std::string OrtoolsBackend::GetSolverVersion() {
  return lp()->SolverVersion();
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
  return 0;
}

double OrtoolsBackend::SimplexIterations() const {
  return 0;
}

int OrtoolsBackend::BarrierIterations() const {
  return 0;
}

void OrtoolsBackend::ExportModel(const std::string &file) {
  // TODO export proper by file extension
  //ORTOOLS_CCALL(ORTOOLS_WriteLp(lp(), file.data()));
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
  int v=-1;
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

void OrtoolsBackend::InitCustomOptions() {

  set_option_header(
      "ORTOOLS Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``ortools_options``. For example::\n"
      "\n"
      "  ampl: option ortools_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

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

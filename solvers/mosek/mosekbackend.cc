#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mosekbackend.h"

extern "C" {
  #include "mosek-ampls-c-api.h"    // Mosek AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptMosek(void* prob) {
  // TODO 
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateMosekBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::MosekBackend()};
}


namespace mp {

/// Create Mosek Model Manager
/// @param gc: the Mosek common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return MosekModelMgr
std::unique_ptr<BasicModelManager>
CreateMosekModelMgr(MosekCommon&, Env&, pre::BasicValuePresolver*&);


MosekBackend::MosekBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateMosekModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

MosekBackend::~MosekBackend() {
  CloseSolver();
}

static void MSKAPI printstr(void* handle,
  const char* str)
{
  fmt::print("{}", str);
  fflush(stdout);
}

void MosekBackend::OpenSolver() {
  int status = MSK_RES_OK;
  MSKtask_t task;
  status = MSK_maketask(NULL, 0, 0, &task);
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create task, error code {}.", status ) );
  set_lp(task); // Assign it
  /// Turn off verbosity by default
  MOSEK_CCALL(MSK_putintparam(task, MSK_IPAR_LOG, 0));

  // Register callback for console logging (controlled by the outlev param
  // in all AMPL solver drivers
  MSK_linkfunctotaskstream(lp(), MSK_STREAM_LOG, NULL, printstr);

}

void MosekBackend::CloseSolver() {
  if (lp()) MSK_deletetask(&lp_ref());
}

const char* MosekBackend::GetBackendName()
  { return "MosekBackend"; }

std::string MosekBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", MSK_VERSION_MAJOR,
    MSK_VERSION_MINOR, MSK_VERSION_REVISION);
}


bool MosekBackend::IsMIP() const {
  return getIntAttr(MSK_IINF_ANA_PRO_NUM_VAR_BIN)+
    getIntAttr(MSK_IINF_ANA_PRO_NUM_VAR_INT);
}

bool MosekBackend::IsQCP() const {
  // TODO
// return getIntAttr(MOSEK_INTATTR_QELEMS) > 0;
  return false;
}

Solution MosekBackend::GetSolution() {
  auto mv = GetValuePresolver().PostsolveSolution(
        { PrimalSolution(), DualSolution() } );
  return { mv.GetVarValues()(), mv.GetConValues()(),
    GetObjectiveValues() };   
}

ArrayRef<double> MosekBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  // TODO get appropriate solution
  MSK_getxx(lp(), MSK_SOL_BAS, x.data());
  return x;
}

pre::ValueMapDbl MosekBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> MosekBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  // TODO get appropriate solution
  MSKrescodee error = MSK_gety(lp(), MSK_SOL_BAS, pi.data());
  if (error != MSK_RESPONSE_OK)
    pi.clear();
  return pi;
}

double MosekBackend::ObjectiveValue() const {
  double v;
  // TODO get appropriate solution
  MOSEK_CCALL(MSK_getprimalobj(lp(), MSK_SOL_BAS, &v));
  return v;
}

double MosekBackend::NodeCount() const {
  return getIntAttr(MSK_IINF_MIO_NUM_ACTIVE_NODES); // TODO check
}

double MosekBackend::SimplexIterations() const {
  return getIntAttr(MSK_IINF_SIM_PRIMAL_ITER);
}

int MosekBackend::BarrierIterations() const {
  return getIntAttr(MSK_IINF_INTPNT_ITER);
}

void MosekBackend::ExportModel(const std::string &file) {
  MOSEK_CCALL(MSK_writedata(lp(), file.data()));
}


void MosekBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptMosek, lp());
  // TODO Check interrupter
}

void MosekBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  MOSEK_CCALL(MSK_optimizetrm(lp(), &termcode));
  MSK_getsolsta(lp(), MSK_SOL_BAS, &solsta); // TODO appropriate sol
  WindupMOSEKSolve();
}

void MosekBackend::WindupMOSEKSolve() { }

void MosekBackend::ReportResults() {
  ReportMOSEKResults();
  BaseBackend::ReportResults();
}

void MosekBackend::ReportMOSEKResults() {
  SetStatus( ConvertMOSEKStatus() );
  AddMOSEKMessages();
  if (need_multiple_solutions())
    ReportMOSEKPool();
}
std::vector<double> MosekBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
  // TODO get solutions in the pool 
  return vars;
}
double MosekBackend::getPoolObjective(int i)
{
  double obj;
  // TODO get objective value of solution i
  return obj;
}
void MosekBackend::ReportMOSEKPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  // TODO
  /*
  while (++iPoolSolution < getIntAttr(MOSEK_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void MosekBackend::AddMOSEKMessages() {
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

std::pair<int, std::string> MosekBackend::ConvertMOSEKStatus() {
  namespace sol = mp::sol;
  // TODO check the logic here
  switch (termcode) {
  case MSK_RES_OK:
    switch (solsta) {
    case MSK_SOL_STA_OPTIMAL:
    case MSK_SOL_STA_INTEGER_OPTIMAL:
      return { sol::SOLVED, "optimal" };
    case MSK_SOL_STA_PRIM_FEAS:
      return { sol::UNCERTAIN, "feasible primal" };
    case MSK_SOL_STA_DUAL_FEAS:
      return { sol::UNCERTAIN, "feasible dual" };
    case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      return { sol::UNCERTAIN, "feasible solution" };
    case MSK_SOL_STA_PRIM_INFEAS_CER:
    case MSK_SOL_STA_DUAL_INFEAS_CER:
      return { sol::INFEASIBLE, "unfeasible solution" };
    case MSK_SOL_STA_UNKNOWN:
    case MSK_SOL_STA_PRIM_ILLPOSED_CER:
    case MSK_SOL_STA_DUAL_ILLPOSED_CER:
    default:
      return { sol::UNKNOWN, "unknown" };
    }
  case MSK_RES_TRM_MAX_ITERATIONS:
    return { sol::LIMIT, "max number of iterations reached" };
  case MSK_RES_TRM_MAX_TIME:
    return { sol::LIMIT, "maximum allowed time reached" };
  case MSK_RES_TRM_MIO_NUM_RELAXS:
  case MSK_RES_TRM_MIO_NUM_BRANCHES:
  case MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS:
    return { sol::LIMIT, "limit hit" };
  default:
    return { sol::UNKNOWN, "unfinished" };
  }
}

void MosekBackend::FinishOptionParsing() {
  int v=-1;
  GetSolverOption(MSK_IPAR_LOG, v);
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

void MosekBackend::InitCustomOptions() {

  set_option_header(
      "MOSEK Optimizer Options for AMPL\n"
      "--------------------------------\n\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``mosek_options``. For example::\n"
      "\n"
      "  ampl: option mosek_options 'threads=3';\n");

  // Example of stored option, to be acted upon in the driver code
  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``..task``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  // Example of direct solver option (set directly by the framework)
  AddSolverOption("tech:threads threads",
    "Controls the number of threads employed by the optimizer. "
    "Default 0 ==> number of threads used will be equal to the number "
    "of cores detected on the machine.",
    MSK_IPAR_NUM_THREADS, 0, INT_MAX);

  AddSolverOption("tech:outlev outlev",
    "0*/1: Whether to write mosek log lines to stdout.",
    MSK_IPAR_LOG, 0, 1);
}

double MosekBackend::MIPGap() {
  return getDblAttr(MSK_DINF_MIO_OBJ_ABS_GAP);
}
double MosekBackend::BestDualBound() {
  // TODO
  return 0;
  //return getDblAttr(MOSEK_DBLATTR_BESTBND);
}

double MosekBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> MosekBackend::VarStatii() {
  // TODO 
  std::vector<int> vars(NumVars());
  /*
  MOSEK_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case MOSEK_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case MOSEK_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case MOSEK_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case MOSEK_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case MOSEK_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Mosek VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> MosekBackend::ConStatii() {
  // TODO
  std::vector<int> cons(NumLinCons());
  /*
  MOSEK_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case MOSEK_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case MOSEK_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case MOSEK_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case MOSEK_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case MOSEK_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Mosek VBasis value: {}", s));
    }
  }*/
  return cons;
}

void MosekBackend::VarStatii(ArrayRef<int> vst) {
  // TODO
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = MOSEK_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = MOSEK_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = MOSEK_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = MOSEK_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = MOSEK_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!MOSEK_GetColInfo(lp(), MOSEK_DBLINFO_LB, 1, index, &lb) && 
        !MOSEK_GetColInfo(lp(), MOSEK_DBLINFO_UB, 1, index, &ub))
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
  MOSEK_SetBasis(lp(), stt.data(), NULL);
  */
}

void MosekBackend::ConStatii(ArrayRef<int> cst) {
  // TODO
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = MOSEK_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = MOSEK_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  MOSEK_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis MosekBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
    assert(constt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void MosekBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

} // namespace mp


// AMPLs

AMPLS_MP_Solver* AMPLSOpenMosek(
  const char* slv_opt) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::MosekBackend()},
    slv_opt);
}

void AMPLSCloseMosek(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

MSKtask_t GetMosekmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::MosekBackend*>(AMPLSGetBackend(slv))->lp();
}

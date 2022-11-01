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

void CbcmpBackend::ExportModel(const std::string &file) {
  // TODO export proper by file extension
  //CBCMP_CCALL(CBCMP_WriteLp(lp(), file.data()));
}


void CbcmpBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCbcmp, lp());
  // TODO Check interrupter
  //CBCMP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CbcmpBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
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

void CbcmpBackend::InitCustomOptions() {

  
  set_option_header(
      "CBCMP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cbcmp_options``. For example::\n"
      "\n"
      "  ampl: option cbcmp_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);


  AddSolverOption("tech:outlev outlev",
    "0*/1: Whether to write log lines (chatter) to stdout and to file.",
    "log", 0, 4);

  _options.addOption("timelimit", Cbc_setMaximumSeconds);

}


double CbcmpBackend::MIPGap() {
  return 0;
//  return getDblAttr(CBCMP_DBLATTR_BESTGAP);
}
double CbcmpBackend::BestDualBound() {
  return 0;
  //return getDblAttr(CBCMP_DBLATTR_BESTBND);
}

double CbcmpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> CbcmpBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  CBCMP_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case CBCMP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case CBCMP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case CBCMP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case CBCMP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case CBCMP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Cbcmp VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> CbcmpBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  CBCMP_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case CBCMP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case CBCMP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case CBCMP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case CBCMP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case CBCMP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Cbcmp VBasis value: {}", s));
    }
  }*/
  return cons;
}

void CbcmpBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = CBCMP_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = CBCMP_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = CBCMP_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = CBCMP_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = CBCMP_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!CBCMP_GetColInfo(lp(), CBCMP_DBLINFO_LB, 1, index, &lb) && 
        !CBCMP_GetColInfo(lp(), CBCMP_DBLINFO_UB, 1, index, &ub))
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
  CBCMP_SetBasis(lp(), stt.data(), NULL);
  */
}

void CbcmpBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = CBCMP_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = CBCMP_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  CBCMP_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis CbcmpBackend::GetBasis() {
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

void CbcmpBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void CbcmpBackend::ComputeIIS() {
  //CBCMP_CCALL(CBCMP_ComputeIIS(lp()));
  SetStatus(ConvertCBCMPStatus());   // could be new information
}

IIS CbcmpBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> CbcmpBackend::VarsIIS() {
  return ArrayRef<int>();
//  return getIIS(lp(), NumVars(), CBCMP_GetColLowerIIS, CBCMP_GetColUpperIIS);
}
pre::ValueMapInt CbcmpBackend::ConsIIS() {
  /*auto iis_lincon = getIIS(lp(), NumLinCons(), CBCMP_GetRowLowerIIS, CBCMP_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  CBCMP_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  CBCMP_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };
      */
  return { {{ 0, std::vector<int>()}} };
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

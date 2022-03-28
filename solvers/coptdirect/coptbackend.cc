#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/backend_model_api_base.h"
#include "coptbackend.h"


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
CreateCoptModelMgr(CoptCommon&, Env&, pre::BasicPresolver*&);


CoptBackend::CoptBackend() {
  OpenSolver();

  pre::BasicPresolver* pPre;
  auto data = CreateCoptModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetPresolver(pPre);

  copy_handlers_to_other_copt();
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
  return getIntAttr(COPT_INTATTR_QELEMS) > 0;
}

Solution CoptBackend::GetSolution() {
  auto mv = GetPresolver().PostsolveSolution(
        { PrimalSolution(), DualSolution() } );
  return { mv.GetVarValues()(), mv.GetConValues()(),
    GetObjectiveValues() };   // TODO postsolve obj values
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
  int error = COPT_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  if (error)
    pi.clear();
  return pi;
}

double CoptBackend::ObjectiveValue() const {
  if (IsMIP())
    return getDblAttr(COPT_DBLATTR_BESTOBJ);
  else
    return getDblAttr(COPT_DBLATTR_LPOBJVAL);
  // TODO Why failsafe below ?
  //double objval = -1e308;
  //CPXgetobjval (env(), lp(), &objval );   // failsafe
  //return objval;
}

double CoptBackend::NodeCount() const {
  return getIntAttr(COPT_INTATTR_NODECNT);
}

double CoptBackend::SimplexIterations() const {
  // TODO check
  return getIntAttr(COPT_INTATTR_SIMPLEXITER);
//  return std::max(
  //      CPXgetmipitcnt (env(), lp()), CPXgetitcnt (env(), lp()));
}

int CoptBackend::BarrierIterations() const {
  return getIntAttr(COPT_INTATTR_BARRIERITER);
}

void CoptBackend::ExportModel(const std::string &file) {
  // TODO export proper by file extension
  COPT_CCALL(COPT_WriteLp(lp(), file.data()));
}


void CoptBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCopt, lp());
  // TODO Check interrupter
  //COPT_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CoptBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
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
}

void CoptBackend::AddCOPTMessages() {
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

std::pair<int, std::string> CoptBackend::ConvertCOPTStatus() {
  namespace sol = mp::sol;
  if (IsMIP())
  {
    int optstatus = getIntAttr(COPT_INTATTR_MIPSTATUS);
    switch (optstatus) {
    case COPT_MIPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case COPT_MIPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case COPT_MIPSTATUS_INF_OR_UNB:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case COPT_MIPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case COPT_MIPSTATUS_TIMEOUT:
    case COPT_MIPSTATUS_NODELIMIT:
    case COPT_MIPSTATUS_INTERRUPTED:
      return { sol::INTERRUPTED, "interrupted" };
    }
  }
  else {
    int optstatus = getIntAttr(COPT_INTATTR_LPSTATUS);
    switch (optstatus) {
    case COPT_LPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case COPT_LPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case COPT_LPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case COPT_LPSTATUS_TIMEOUT:
      return { sol::INTERRUPTED, "interrupted" };
    default:
      return { sol::UNKNOWN, "unfinished" };
    }
  }
    /*
  int optimstatus = CPXgetstat(env(), lp());
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter()->Stop()) {
      return { sol::INTERRUPTED, "interrupted" };
    }
    int solcount;
    solcount = CPXgetsolnpoolnumsolns (env(), lp());  // Can we use it without CPXpopulate?
    if (solcount>0) {
      return { sol::UNCERTAIN, "feasible solution" };
    }
    return { sol::UNKNOWN, "unknown solution status" };
  case CPX_STAT_OPTIMAL:
  case CPXMIP_OPTIMAL:
  case CPX_STAT_MULTIOBJ_OPTIMAL:
    return { sol::SOLVED, "optimal solution" };
  case CPX_STAT_INFEASIBLE:
  case CPXMIP_INFEASIBLE:
  case CPX_STAT_MULTIOBJ_INFEASIBLE:
    return { sol::INFEASIBLE, "infeasible problem" };
  case CPX_STAT_INForUNBD:
  case CPXMIP_INForUNBD:
  case CPX_STAT_MULTIOBJ_INForUNBD:
    return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
  case CPX_STAT_UNBOUNDED:
  case CPXMIP_UNBOUNDED:
  case CPX_STAT_MULTIOBJ_UNBOUNDED:
    return { sol::UNBOUNDED, "unbounded problem" };
  case CPX_STAT_FEASIBLE_RELAXED_INF:
  case CPX_STAT_FEASIBLE_RELAXED_QUAD:
  case CPX_STAT_FEASIBLE_RELAXED_SUM:
  case CPX_STAT_NUM_BEST:
  case CPX_STAT_OPTIMAL_INFEAS:
  case CPX_STAT_OPTIMAL_RELAXED_INF:
  case CPX_STAT_OPTIMAL_RELAXED_QUAD:
  case CPX_STAT_OPTIMAL_RELAXED_SUM:
    return { sol::UNCERTAIN, "feasible or optimal but numeric issue" };
  }
  */
}


void CoptBackend::FinishOptionParsing() {
  int v=-1;
  GetSolverOption(COPT_INTPARAM_LOGGING, v);
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////

void CoptBackend::InitCustomOptions() {

  set_option_header(
      "IBM ILOG COPT Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``copt_options``. For example::\n"
      "\n"
      "  ampl: option copt_options 'mipgap=1e-6';\n");

  AddSolverOption("tech:outlev outlev",
      "0-1: output logging verbosity. "
      "Default = 0 (no logging).",
    COPT_INTPARAM_LOGTOCONSOLE, 0, 1);
  SetSolverOption(COPT_INTPARAM_LOGTOCONSOLE, 0);

  AddStoredOption("tech:exportfile writeprob",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);
  
  AddSolverOption("mip:gap mipgap",
      "Relative optimality gap |bestbound-bestinteger|/(1e-10+|bestinteger|).",
        COPT_DBLPARAM_RELGAP, 0.0, 1.0);

  AddSolverOption("tech:threads threads",
    "Number of threads to use;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_BARTHREADS, -1, 128);

  AddSolverOption("tech:barrierthreads barthreads",
      "Number of threads used by the barrier algorithm;\n"
      "default -1 ==> see use value in tech:threads.",
    COPT_INTPARAM_BARTHREADS, -1, 128);

  AddSolverOption("tech:crossoverthreads crossoverthreads",
    "Number of threads used by crossover;\n"
    "default -1 ==> see use value in tech:threads.",
    COPT_INTPARAM_CROSSOVERTHREADS, -1, 128);

  AddSolverOption("tech:simplexthreads simplexthreads",
    "Number of threads used by dual simplex;\n"
    "default -1 ==> see use value in tech:threads.",
    COPT_INTPARAM_SIMPLEXTHREADS, -1, 128);


  AddSolverOption("tech:miptasks miptasks",
    "Number of MIP tasks in parallel;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_MIPTASKS, -1, 255);


  AddSolverOption("lim:time timelim timelimit",
      "limit on solve time (in seconds; default: no limit).",
      COPT_DBLPARAM_TIMELIMIT, 0.0, DBL_MAX);
      
}


} // namespace mp

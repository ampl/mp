#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "gcgbackend.h"

extern "C" {
  #include "gcg-ampls-c-api.h"    // Gcg AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptGcg(void* prob) {
  //return GCG_Interrupt((gcg_prob*)prob);
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

void GcgBackend::OpenSolver() {
  int status = 0;
  // TODO Typically this function creates an instance of the solver environment
  // and an empty model
  void* env_p;
  // Typically try the registered function first;
  // if not available call the solver's API function directly
  /*
  const auto& create_fn = GetCallbacks().cb_initsolver_;
  if (create_fn)
    set_env((GRBenv*)create_fn());
  else
    status = createEnv(&env_p);
    */
  // set_env(env_p);

  /* Todo catch errors
  if ( env() == NULL ) {
    // char  errmsg[CPXMESSAGEBUFSIZE];
    // CPXgeterrorstring (env(), status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open GCG environment.\n{}", status) );
  }
  */

  /* TODO Create problem instance
  gcg_prob* prob;
  status = GCG_CreateProb(env_p, &prob);
 */
  Solver::SolverModel* prob = Solver::CreateSolverModel();
  set_lp(prob); // Assign it
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  /* TODO Typically check call */
  /// Turn off verbosity by default
  // GCG_CCALL(GCG_SetIntParam(prob, "Logging", 0));

}

void GcgBackend::CloseSolver() {
  /* TODO Cleanup: close problem and environment
  if ( lp() != NULL ) {
    GCG_CCALL(GCG_DeleteProb(&lp_) );
  }
  if ( env() != NULL ) {
    GCG_CCALL(GCG_DeleteEnv(&env_) );
  }
  */
}

const char* GcgBackend::GetBackendName()
  { return "GcgBackend"; }

std::string GcgBackend::GetSolverVersion() {
  // TODO Return version from solver API
  return "0.0.0";
  //return fmt::format("{}.{}.{}", GCG_VERSION_MAJOR, 
  //  GCG_VERSION_MINOR, GCG_VERSION_TECHNICAL);
}


bool GcgBackend::IsMIP() const {
  // TODO
  return getIntAttr(Solver::VARS_INT) > 0;
  //return getIntAttr(GCG_INTATTR_ISMIP);
}

bool GcgBackend::IsQCP() const {
  return getIntAttr(Solver::CONS_QUAD) > 0;
// return getIntAttr(GCG_INTATTR_QELEMS) > 0;
}

ArrayRef<double> GcgBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  /*
  if (IsMIP()) 
    error = GCG_GetSolution(lp(), x.data());
  else
    error = GCG_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);
  if (error)
    x.clear();
    */
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
 /* if (IsMIP())
    return getDblAttr(GCG_DBLATTR_BESTOBJ);
  else
    return getDblAttr(GCG_DBLATTR_LPOBJVAL);
    */
  return 0;
}

double GcgBackend::NodeCount() const {
  return 0;
//  return getIntAttr(GCG_INTATTR_NODECNT);
}

double GcgBackend::SimplexIterations() const {
  return 0;
//  return getIntAttr(GCG_INTATTR_SIMPLEXITER);
}

int GcgBackend::BarrierIterations() const {
  return 0;
//  return getIntAttr(GCG_INTATTR_BARRIERITER);
}

void GcgBackend::ExportModel(const std::string &file) {
  // TODO export proper by file extension
  //GCG_CCALL(GCG_WriteLp(lp(), file.data()));
}


void GcgBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptGcg, lp());
  // TODO Check interrupter
  //GCG_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void GcgBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  //GCG_CCALL(GCG_Solve(lp()));
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
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> GcgBackend::ConvertGCGStatus() {
  namespace sol = mp::sol;
  if (IsMIP())
  {
    /*
    int optstatus = getIntAttr(GCG_INTATTR_MIPSTATUS);
    switch (optstatus) {
    case GCG_MIPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case GCG_MIPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case GCG_MIPSTATUS_INF_OR_UNB:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case GCG_MIPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case GCG_MIPSTATUS_TIMEOUT:
    case GCG_MIPSTATUS_NODELIMIT:
    case GCG_MIPSTATUS_INTERRUPTED:
      return { sol::INTERRUPTED, "interrupted" };
    }
    */
  }
  else {
    /*
    int optstatus = getIntAttr(GCG_INTATTR_LPSTATUS);
    switch (optstatus) {
    case GCG_LPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case GCG_LPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case GCG_LPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case GCG_LPSTATUS_TIMEOUT:
      return { sol::INTERRUPTED, "interrupted" };
    default:
      return { sol::UNKNOWN, "unfinished" };
    }
    */
  }
  return { sol::UNKNOWN, "not solved" };
}


void GcgBackend::FinishOptionParsing() {
  int v=-1;
 // GetSolverOption(GCG_INTPARAM_LOGGING, v);
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

void GcgBackend::InitCustomOptions() {

  set_option_header(
      "GCG Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``gcg_options``. For example::\n"
      "\n"
      "  ampl: option gcg_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

}


double GcgBackend::MIPGap() {
  return 0;
//  return getDblAttr(GCG_DBLATTR_BESTGAP);
}
double GcgBackend::BestDualBound() {
  return 0;
  //return getDblAttr(GCG_DBLATTR_BESTBND);
}

double GcgBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
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

void GcgBackend::AddMIPStart(ArrayRef<double> x0) {
  //GCG_CCALL(GCG_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));

}


} // namespace mp


// AMPLs
void* AMPLSOpenGcg(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::GcgBackend()},
    slv_opt, cb);
}

void AMPLSCloseGcg(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetGcgmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::GcgBackend*>(AMPLSGetBackend(slv))->lp();
}

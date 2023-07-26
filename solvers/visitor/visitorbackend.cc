#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "visitorbackend.h"

extern "C" {
  #include "visitor-ampls-c-api.h"    // Visitor AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptVisitor(void* prob) {
  //return VISITOR_Interrupt((visitor_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateVisitorBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::VisitorBackend()};
}


namespace mp {

/// Create Visitor Model Manager
/// @param gc: the Visitor common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return VisitorModelMgr
std::unique_ptr<BasicModelManager>
CreateVisitorModelMgr(VisitorCommon&, Env&, pre::BasicValuePresolver*&);


VisitorBackend::VisitorBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateVisitorModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

VisitorBackend::~VisitorBackend() {
  CloseSolver();
}

void VisitorBackend::OpenSolver() {
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
       fmt::format("Could not open VISITOR environment.\n{}", status) );
  }
  */

  /* TODO Create problem instance
  visitor_prob* prob;
  status = VISITOR_CreateProb(env_p, &prob);
 */
  Solver::SolverModel* prob = Solver::CreateSolverModel();
  set_lp(prob); // Assign it
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  /* TODO Typically check call */
  /// Turn off verbosity by default
  // VISITOR_CCALL(VISITOR_SetIntParam(prob, "Logging", 0));

}

void VisitorBackend::CloseSolver() {
  /* TODO Cleanup: close problem and environment
  if ( lp() != NULL ) {
    VISITOR_CCALL(VISITOR_DeleteProb(&lp_) );
  }
  if ( env() != NULL ) {
    VISITOR_CCALL(VISITOR_DeleteEnv(&env_) );
  }
  */
}

const char* VisitorBackend::GetBackendName()
  { return "VisitorBackend"; }

std::string VisitorBackend::GetSolverVersion() {
  // TODO Return version from solver API
  return "0.0.0";
  //return fmt::format("{}.{}.{}", VISITOR_VERSION_MAJOR, 
  //  VISITOR_VERSION_MINOR, VISITOR_VERSION_TECHNICAL);
}


bool VisitorBackend::IsMIP() const {
  // TODO
  return getIntAttr(Solver::VARS_INT) > 0;
  //return getIntAttr(VISITOR_INTATTR_ISMIP);
}

bool VisitorBackend::IsQCP() const {
  return getIntAttr(Solver::CONS_QUAD) > 0;
// return getIntAttr(VISITOR_INTATTR_QELEMS) > 0;
}

ArrayRef<double> VisitorBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  /*
  if (IsMIP()) 
    error = VISITOR_GetSolution(lp(), x.data());
  else
    error = VISITOR_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);
  if (error)
    x.clear();
    */
  return x;
}

pre::ValueMapDbl VisitorBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> VisitorBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
 // int error = VISITOR_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 0;
  if (error)
    pi.clear();
  return pi;
}

double VisitorBackend::ObjectiveValue() const {
 /* if (IsMIP())
    return getDblAttr(VISITOR_DBLATTR_BESTOBJ);
  else
    return getDblAttr(VISITOR_DBLATTR_LPOBJVAL);
    */
  return 0;
}

double VisitorBackend::NodeCount() const {
  return 0;
//  return getIntAttr(VISITOR_INTATTR_NODECNT);
}

double VisitorBackend::SimplexIterations() const {
  return 0;
//  return getIntAttr(VISITOR_INTATTR_SIMPLEXITER);
}

int VisitorBackend::BarrierIterations() const {
  return 0;
//  return getIntAttr(VISITOR_INTATTR_BARRIERITER);
}

void VisitorBackend::ExportModel(const std::string &file) {
  // TODO export proper by file extension
  //VISITOR_CCALL(VISITOR_WriteLp(lp(), file.data()));
}


void VisitorBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptVisitor, lp());
  // TODO Check interrupter
  //VISITOR_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void VisitorBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  //VISITOR_CCALL(VISITOR_Solve(lp()));
  WindupVISITORSolve();
}

void VisitorBackend::WindupVISITORSolve() { }

void VisitorBackend::ReportResults() {
  ReportVISITORResults();
  BaseBackend::ReportResults();
}

void VisitorBackend::ReportVISITORResults() {
  SetStatus( ConvertVISITORStatus() );
  AddVISITORMessages();
  if (need_multiple_solutions())
    ReportVISITORPool();
}
std::vector<double> VisitorBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // VISITOR_CCALL(VISITOR_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double VisitorBackend::getPoolObjective(int i)
{
  double obj;
 // VISITOR_CCALL(VISITOR_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void VisitorBackend::ReportVISITORPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(VISITOR_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void VisitorBackend::AddVISITORMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> VisitorBackend::ConvertVISITORStatus() {
  namespace sol = mp::sol;
  if (IsMIP())
  {
    /*
    int optstatus = getIntAttr(VISITOR_INTATTR_MIPSTATUS);
    switch (optstatus) {
    case VISITOR_MIPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case VISITOR_MIPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case VISITOR_MIPSTATUS_INF_OR_UNB:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case VISITOR_MIPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case VISITOR_MIPSTATUS_TIMEOUT:
    case VISITOR_MIPSTATUS_NODELIMIT:
    case VISITOR_MIPSTATUS_INTERRUPTED:
      return { sol::INTERRUPTED, "interrupted" };
    }
    */
  }
  else {
    /*
    int optstatus = getIntAttr(VISITOR_INTATTR_LPSTATUS);
    switch (optstatus) {
    case VISITOR_LPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case VISITOR_LPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case VISITOR_LPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case VISITOR_LPSTATUS_TIMEOUT:
      return { sol::INTERRUPTED, "interrupted" };
    default:
      return { sol::UNKNOWN, "unfinished" };
    }
    */
  }
  return { sol::UNKNOWN, "not solved" };
}


void VisitorBackend::FinishOptionParsing() {
  int v=-1;
 // GetSolverOption(VISITOR_INTPARAM_LOGGING, v);
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

void VisitorBackend::InitCustomOptions() {

  set_option_header(
      "VISITOR Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``visitor_options``. For example::\n"
      "\n"
      "  ampl: option visitor_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

}


double VisitorBackend::MIPGap() {
  return 0;
//  return getDblAttr(VISITOR_DBLATTR_BESTGAP);
}
double VisitorBackend::BestDualBound() {
  return 0;
  //return getDblAttr(VISITOR_DBLATTR_BESTBND);
}

double VisitorBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> VisitorBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  VISITOR_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case VISITOR_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case VISITOR_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case VISITOR_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case VISITOR_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case VISITOR_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Visitor VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> VisitorBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  VISITOR_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case VISITOR_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case VISITOR_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case VISITOR_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case VISITOR_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case VISITOR_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Visitor VBasis value: {}", s));
    }
  }*/
  return cons;
}

void VisitorBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = VISITOR_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = VISITOR_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = VISITOR_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = VISITOR_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = VISITOR_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!VISITOR_GetColInfo(lp(), VISITOR_DBLINFO_LB, 1, index, &lb) && 
        !VISITOR_GetColInfo(lp(), VISITOR_DBLINFO_UB, 1, index, &ub))
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
  VISITOR_SetBasis(lp(), stt.data(), NULL);
  */
}

void VisitorBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = VISITOR_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = VISITOR_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  VISITOR_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis VisitorBackend::GetBasis() {
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

void VisitorBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void VisitorBackend::ComputeIIS() {
  //VISITOR_CCALL(VISITOR_ComputeIIS(lp()));
  SetStatus(ConvertVISITORStatus());   // could be new information
}

IIS VisitorBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> VisitorBackend::VarsIIS() {
  return ArrayRef<int>();
//  return getIIS(lp(), NumVars(), VISITOR_GetColLowerIIS, VISITOR_GetColUpperIIS);
}
pre::ValueMapInt VisitorBackend::ConsIIS() {
  /*auto iis_lincon = getIIS(lp(), NumLinCons(), VISITOR_GetRowLowerIIS, VISITOR_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  VISITOR_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  VISITOR_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };
      */
  return { {{ 0, std::vector<int>()}} };
}

void VisitorBackend::AddMIPStart(
    ArrayRef<double> x0, ArrayRef<int> sparsity) {
  //VISITOR_CCALL(VISITOR_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));
}


} // namespace mp


// AMPLs
void* AMPLSOpenVisitor(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::VisitorBackend()},
        cb);
}

void AMPLSCloseVisitor(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetVisitormodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::VisitorBackend*>(AMPLSGetBackend(slv))->lp();
}

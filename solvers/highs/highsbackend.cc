#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "highsbackend.h"

extern "C" {
  #include "highs-ampls-c-api.h"    // Highs AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptHighs(void* prob) {
  //return HIGHS_Interrupt((highs_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateHighsBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()};
}


namespace mp {

/// Create Highs Model Manager
/// @param gc: the Highs common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return HighsModelMgr
std::unique_ptr<BasicModelManager>
CreateHighsModelMgr(HighsCommon&, Env&, pre::BasicPresolver*&);


HighsBackend::HighsBackend() {
  OpenSolver();

  pre::BasicPresolver* pPre;
  auto data = CreateHighsModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetPresolver(pPre);

  copy_handlers_to_other_highs();
}

HighsBackend::~HighsBackend() {
  CloseSolver();
}

const char* HighsBackend::GetBackendName()
  { return "HighsBackend"; }

std::string HighsBackend::GetSolverVersion() {
  // TODO Return version from solver API
  return "0.0.0";
  //return fmt::format("{}.{}.{}", HIGHS_VERSION_MAJOR, 
  //  HIGHS_VERSION_MINOR, HIGHS_VERSION_TECHNICAL);
}


bool HighsBackend::IsMIP() const {
  // TODO
  return false;
}

bool HighsBackend::IsQCP() const {
  return false; // TODO
}

Solution HighsBackend::GetSolution() {
  auto mv = GetPresolver().PostsolveSolution(
        { PrimalSolution(), DualSolution() } );
  return { mv.GetVarValues()(), mv.GetConValues()(),
    GetObjectiveValues() };   // TODO postsolve obj values
}

ArrayRef<double> HighsBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  Highs_getSolution(lp(), x.data(), NULL, NULL, NULL);
  return x;
}

pre::ValueMapDbl HighsBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> HighsBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  Highs_getSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 0;
  if (error)
    pi.clear();
  return pi;
}

double HighsBackend::ObjectiveValue() const {
  return Highs_getObjectiveValue(lp());
}

double HighsBackend::NodeCount() const {
  return getInt64Attr("mip_node_count");
}

double HighsBackend::SimplexIterations() const {
  return getIntAttr("simplex_iteration_count");
}

int HighsBackend::BarrierIterations() const {
  return getIntAttr("ipm_iteration_count");
}

void HighsBackend::ExportModel(const std::string &file) {
  HIGHS_CCALL(Highs_writeModel(lp(), file.data()));
}


void HighsBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptHighs, lp());
  // TODO Check interrupter
  //HIGHS_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void HighsBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  HIGHS_CCALL(Highs_run(lp()));
  WindupHIGHSSolve();
}

void HighsBackend::WindupHIGHSSolve() { }

void HighsBackend::ReportResults() {
  ReportHIGHSResults();
  BaseBackend::ReportResults();
}

void HighsBackend::ReportHIGHSResults() {
  SetStatus( ConvertHIGHSStatus() );
  AddHIGHSMessages();
  if (need_multiple_solutions())
    ReportHIGHSPool();
}
std::vector<double> HighsBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // HIGHS_CCALL(HIGHS_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double HighsBackend::getPoolObjective(int i)
{
  double obj;
 // HIGHS_CCALL(HIGHS_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void HighsBackend::ReportHIGHSPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(HIGHS_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


SolutionBasis HighsBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetPresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
    assert(constt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void HighsBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetPresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

ArrayRef<int> HighsBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  HIGHS_CCALL(Highs_getBasis(lp(), vars.data(), NULL));
  for (auto& s : vars) {
    switch (s) {
    case kHighsBasisStatusBasic:
      s = (int)BasicStatus::bas;
      break;
    case kHighsBasisStatusLower:
      s = (int)BasicStatus::low;
      break;
    case kHighsBasisStatusUpper:
      s = (int)BasicStatus::upp;
      break;
    case kHighsBasisStatusNonbasic:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Highs VBasis value: {}", s));
    }
  }
  return vars;
}

ArrayRef<int> HighsBackend::ConStatii() {
  std::vector<int> cons(NumLinCons());

  HIGHS_CCALL(Highs_getBasis(lp(), NULL, cons.data()));
  for (auto& s : cons) {
    switch (s) {
    case kHighsBasisStatusBasic:
      s = (int)BasicStatus::bas;
      break;
    case kHighsBasisStatusLower:
      s = (int)BasicStatus::low;
      break;
    case kHighsBasisStatusUpper:
      s = (int)BasicStatus::upp;
      break;
    case kHighsBasisStatusNonbasic:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Highs VBasis value: {}", s));
    }
  }
  return cons;
}

void HighsBackend::VarStatii(ArrayRef<int> vst) {
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  std::vector<int> indicesOfMissing;
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = kHighsBasisStatusBasic;
      break;
    case BasicStatus::low:
    case BasicStatus::equ:
      s = kHighsBasisStatusLower;
      break;
    case BasicStatus::upp:
      s = kHighsBasisStatusUpper;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = kHighsBasisStatusNonbasic;
      break;
    case BasicStatus::none:
      indicesOfMissing.push_back(j);
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }

  if (indicesOfMissing.size() > 0)
  {
    /// 'none' is assigned to new variables. Compute low/upp/sup:
    /// Depending on where 0.0 is between bounds
    std::vector<double> lb(indicesOfMissing.size());
    std::vector<double> ub(indicesOfMissing.size());
    Highs_getColsBySet(lp(), indicesOfMissing.size(), indicesOfMissing.data(),
      NULL, NULL, lb.data(), ub.data(), NULL, NULL, NULL, NULL);
    for (int i = 0; i < indicesOfMissing.size(); i++) {
     
        if (lb[i] >= -1e-6)
          stt[indicesOfMissing[i]] = kHighsBasisStatusLower;
        else if (ub[i] <= 1e-6)
          stt[indicesOfMissing[i]] = kHighsBasisStatusUpper;
        else
          stt[indicesOfMissing[i]] = kHighsBasisStatusNonbasic;
      }
  }
 HIGHS_CCALL(Highs_setBasis(lp(), stt.data(), NULL));
}

void HighsBackend::ConStatii(ArrayRef<int> cst) {
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = kHighsBasisStatusBasic;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status,
    case BasicStatus::low:    // as Gurobi 9.5 does not accept partial basis.
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = kHighsBasisStatusNonbasic;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  HIGHS_CCALL(Highs_setBasis(lp(), NULL, stt.data()));
}

void HighsBackend::AddHIGHSMessages() {
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

std::pair<int, std::string> HighsBackend::ConvertHIGHSStatus() {
  namespace sol = mp::sol;
 // TODO Complete
  int optstatus = Highs_getModelStatus(lp());
    switch (optstatus) {
    case kHighsModelStatusOptimal:
      return { sol::SOLVED, "optimal solution" };
    case kHighsModelStatusInfeasible:
      return { sol::INFEASIBLE, "infeasible problem" };
    case kHighsModelStatusUnbounded:
      return { sol::UNBOUNDED, "unbounded problem" };
    case kHighsModelStatusTimeLimit:
    case kHighsModelStatusIterationLimit:
      return { sol::INTERRUPTED, "interrupted" };
    default:
      return { sol::UNKNOWN, "unfinished" };
    }
  return { sol::UNKNOWN, "not solved" };
}


void HighsBackend::FinishOptionParsing() {
  int v=-1;
 // GetSolverOption(HIGHS_INTPARAM_LOGGING, v);
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo lp_values_method[] = {
  { "choiche", "Automatic (default)", -1},
  { "simplex", "Simplex", 1},
  { "ipm", "Interior point method", 2},
};
/*

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
*/
void HighsBackend::InitCustomOptions() {

  set_option_header(
      "HIGHS Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``highs_options``. For example::\n"
      "\n"
      "  ampl: option highs_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  std::string c = "choice";
  AddSolverOption("alg:method method lpmethod",
    "Which algorithm to use :\n"
    "\n.. value-table::\n", "solver", lp_values_method, c);


  AddSolverOption("tech:threads threads",
    "How many threads to use when using the barrier algorithm "
    "or solving MIP problems; default 0 ==> automatic choice.",
    "threads", 0, 128);

   
}


double HighsBackend::MIPGap() {
  return 0;
//  return getDblAttr(HIGHS_DBLATTR_BESTGAP);
}
double HighsBackend::BestDualBound() {
  return 0;
  //return getDblAttr(HIGHS_DBLATTR_BESTBND);
}

double HighsBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}

void HighsBackend::AddMIPStart(ArrayRef<double> x0) {
  //HIGHS_CCALL(HIGHS_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));

}


} // namespace mp


// AMPLs

AMPLS_MP_Solver* AMPLSOpenHighs(
  const char* slv_opt) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()},
    slv_opt);
}

void AMPLSCloseHighs(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetHighsmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::HighsBackend*>(AMPLSGetBackend(slv))->lp();
}

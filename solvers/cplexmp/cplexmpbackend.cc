#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cplexmpbackend.h"

extern "C" {
#include "cplexmp-ampls-c-api.h"    // CPLEX AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {

volatile int terminate_flag = 0;
bool InterruptCplex(void *) {
  terminate_flag = 1;
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCplexBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CplexBackend()};
}


namespace mp {

  /// Create Cplex Model Manager
  /// @param gc: the Cplex Backend
  /// @param e: environment
  /// @param pre: presolver to be returned,
  /// need it to convert solution data
  /// @return GurobiModelMgr
  std::unique_ptr<BasicModelManager>
    CreateCplexModelMgr(CplexCommon&, Env&, pre::BasicValuePresolver*&);


  void CplexBackend::InitOptionParsing() {
    OpenSolver();
  }
  CplexBackend::CplexBackend() {

    pre::BasicValuePresolver* pPre;
    auto data = CreateCplexModelMgr(*this, *this, pPre);
    SetMM(std::move(data));
    SetValuePresolver(pPre);

    /// Copy env/lp to ModelAPI
    copy_common_info_to_other();
  }

  CplexBackend::~CplexBackend() {
    CloseSolver();
  }

  void CplexBackend::OpenSolver() {
    int status = 0;
    // Typically try the registered function first;
    // if not available call the solver's API function directly
    const auto create_fn = GetCallbacks().init;
    if (create_fn)
      set_env((CPXENVptr)create_fn());
    else
      set_env(CPXopenCPLEX(&status));
    if (env() == NULL) {
      char  errmsg[CPXMESSAGEBUFSIZE];
      CPXgeterrorstring(env(), status, errmsg);
      throw std::runtime_error(
        fmt::format("Could not open CPLEX environment.\n{}", errmsg));
    }
    /* Avoid most error messages on screen */
    CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ScreenOutput, CPX_OFF));
    /* Do not echo the params twice */
    CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ParamDisplay, 0)); 
    /* defaults */
    CPXsetintparam(env(), CPX_PARAM_SIMDISPLAY, 0);
    CPXsetintparam(env(), CPX_PARAM_MIPDISPLAY, 0);
    CPXsetintparam(env(), CPX_PARAM_BARDISPLAY, 0);

    /* Create an empty model */
    set_lp(CPXcreateprob(env(), &status, "amplcplex"));
    if (status)
      throw std::runtime_error(fmt::format(
        "Failed to create problem, error code {}.", status));
    /* Copy handlers to ModelAPI */
    copy_common_info_to_other();
  }

  void CplexBackend::CloseSolver() {
    if (lp() != NULL) {
      CPLEX_CALL(CPXfreeprob(env(), &lp_ref()));
    }
    /* Free up the CPLEX environment, if necessary */
    if (env() != NULL) {
      CPLEX_CALL(CPXcloseCPLEX(&env_ref()));
    }
  }


  const char* CplexBackend::GetBackendName()
  {
    return "CplexBackend";
  }

  std::string CplexBackend::GetSolverVersion() {
    return fmt::format("{}.{}.{}", CPX_VERSION_VERSION,
      CPX_VERSION_RELEASE, CPX_VERSION_MODIFICATION);
  }


  bool CplexBackend::IsMIP() const {
    int probtype = CPXgetprobtype(env(), lp());
    return
      CPXPROB_MILP == probtype ||
      CPXPROB_MIQP == probtype ||
      CPXPROB_MIQCP == probtype;
  }

  bool CplexBackend::IsQCP() const {
    int probtype = CPXgetprobtype(env(), lp());
    return probtype >= 5;
  }
#define getAndReturnDblParam(function)\
  double value;\
  CPLEX_CALL(function(env(), lp(), &value));\
  return value;

#define getDblParam(function, var)\
  double var;\
  CPLEX_CALL(function(env(), lp(), &var));

  double  CplexBackend::MIPGap() {
    getAndReturnDblParam(CPXgetmiprelgap);
  }

  double  CplexBackend::MIPGapAbs() {
    double obj;
    int status = CPXgetobjval(env(), lp(), &obj);
    if (status)
      return AMPLInf(); // no solution found
    return std::abs(obj - BestDualBound());
  }
  double  CplexBackend::BestDualBound() {
    getDblParam(CPXgetbestobjval, bobj);
    if (bobj == Infinity())
      return AMPLInf();
    if (bobj == -Infinity())
      return -AMPLInf();
    return bobj;
  }

  ArrayRef<int> CplexBackend::VarStatii() {
    std::vector<int> vars(NumVars());
    CPLEX_CALL(CPXgetbase(env(), lp(), vars.data(), nullptr));
    for (auto& s : vars) {
      switch (s) {
      case CPX_BASIC:
        s = (int)BasicStatus::bas;
        break;
      case CPX_AT_LOWER:
        s = (int)BasicStatus::low;
        break;
      case CPX_AT_UPPER:
        s = (int)BasicStatus::upp;
        break;
      case CPX_FREE_SUPER:
        s = (int)BasicStatus::sup;
        break;
      default:
        MP_RAISE(fmt::format("Unknown CPLEX cstat value: {}", s));
      }
    }
    return vars;
  }

  ArrayRef<int> CplexBackend::ConStatii() {
    std::vector<int> vars(NumLinCons());
    for (auto& s : vars) {
      switch (s) {
      case CPX_BASIC:
        s = (int)BasicStatus::bas;
        break;
      case CPX_AT_LOWER:
        s = (int)BasicStatus::low;
        break;
      case CPX_AT_UPPER: // just for range constraints
        s = (int)BasicStatus::upp;
        break;
      default:
        MP_RAISE(fmt::format("Unknown CPLEX rstat value: {}", s));
      }
    }
    return vars;
  }

  void CplexBackend::VarConStatii(ArrayRef<int> vstt, ArrayRef<int> cstt) {
    std::vector<int> vst= std::vector<int>(vstt.data(), vstt.data() + vstt.size());
    std::vector<int> cst = std::vector<int>(cstt.data(), cstt.data() + cstt.size());
    for (auto j = vst.size(); j--; ) {
      auto& s = vst[j];
      switch ((BasicStatus)s) {
      case BasicStatus::bas:
        s = (int)CPX_BASIC;
        break;
      case BasicStatus::low:
      case BasicStatus::equ:
      case BasicStatus::none:
        s = CPX_AT_LOWER;
        break;
      case BasicStatus::upp:
        s = CPX_AT_UPPER;
        break;
      case BasicStatus::sup:
      case BasicStatus::btw:
        s = CPX_FREE_SUPER;
        break;
      default:
        MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
      }
    }
      for (auto j = cst.size(); j--; ) {
        auto& s = cst[j];
        switch ((BasicStatus)s) {
        case BasicStatus::bas:
          s = (int)CPX_BASIC;
          break;
        case BasicStatus::low:
        case BasicStatus::equ:
        case BasicStatus::none:
        case BasicStatus::sup:
        case BasicStatus::btw:
          s = CPX_AT_LOWER;
          break;
        case BasicStatus::upp:
          s = CPX_AT_UPPER;
          break;
        default:
          MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
        }
    }
      CPLEX_CALL(CPXcopybase(env(), lp(), vst.data(), cst.data()));
  }


  SolutionBasis CplexBackend::GetBasis() {
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
  void CplexBackend::SetBasis(SolutionBasis basis) {
    auto mv = GetValuePresolver().PresolveBasis(
      { basis.varstt, basis.constt });
    auto varstt = mv.GetVarValues()();
    auto constt = mv.GetConValues()(CG_Linear);
    assert(varstt.size());
    assert(constt.size());
    VarConStatii(varstt, constt);
  }

ArrayRef<double> CplexBackend::PrimalSolution() {
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  int error = CPXgetx (env(), lp(), x.data(), 0, num_vars-1);
  if (error)
    x.clear();
  return x;
}

pre::ValueMapDbl CplexBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> CplexBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  if (IsMIP())
    return pi; // when implementing fixed model, get rid of this clause
  int error = CPXgetpi (env(), lp(), pi.data(), 0, num_cons-1);
  if (error)
    pi.clear();
  return pi;
}

double CplexBackend::ObjectiveValue() const {
  double objval = -Infinity();
  CPXgetobjval (env(), lp(), &objval );
  
  return objval;
}

double CplexBackend::NodeCount() const {
  return CPXgetnodecnt (env(), lp());
}

double CplexBackend::SimplexIterations() const {
  if (IsMIP())
    return CPXgetmipitcnt(env(), lp());
  else
    return CPXgetitcnt(env(), lp());
}

int CplexBackend::BarrierIterations() const {
  return CPXgetbaritcnt (env(), lp());
}

void CplexBackend::DoWriteProblem(const std::string &file) {
  CPLEX_CALL( CPXwriteprob (env(), lp(), file.c_str(), NULL) );
}


void CplexBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCplex, nullptr);
  CPLEX_CALL( CPXsetterminate (env(), &terminate_flag) );
}

void CplexBackend::Solve() {
  if(IsMIP())
    CPLEX_CALL(CPXmipopt(env(), lp()));
  else
    CPLEX_CALL(CPXlpopt(env(), lp()));

  WindupCPLEXSolve();
}

void CplexBackend::WindupCPLEXSolve() { }

void CplexBackend::ReportResults() {
  ReportCPLEXResults();
  BaseBackend::ReportResults();
}

void CplexBackend::ReportCPLEXResults() {
  SetStatus( ConvertCPLEXStatus() );
  AddCPLEXMessages();
  if (need_multiple_solutions())
    ReportCPLEXPool();
}
void CplexBackend::ReportCPLEXPool() {
  if (!IsMIP())
    return;
  if(storedOptions_.populate_==1) CPLEX_CALL(CPXpopulate(env(), lp()));

  int nsols = CPXgetsolnpoolnumsolns(env(), lp());

  double inttol = GetCPLEXDblParam(CPX_PARAM_EPINT);
  double feastol = GetCPLEXDblParam(CPX_PARAM_EPRHS);
  if (inttol > feastol) SetCPLEXParam(CPX_PARAM_EPRHS, inttol);

  int iPoolSolution = -1;
  double currentObj;
  int NumVarsm1 = NumVars();
  std::vector<double> x(NumVars());
  NumVarsm1--; // To use as limit when getting the solution
  while (++iPoolSolution < nsols) {
    
    CPLEX_CALL(CPXgetsolnpoolobjval(env(), lp(), iPoolSolution, &currentObj));
    CPLEX_CALL(CPXgetsolnpoolx(env(), lp(), iPoolSolution, x.data(), 0, NumVarsm1));

    auto mv = GetValuePresolver().PostsolveSolution(  // only single-obj with pool
      { { x },
        {},                                       // no duals
        std::vector<double>{ currentObj } });
    ReportIntermediateSolution(
      { mv.GetVarValues()(), mv.GetConValues()(),
        mv.GetObjValues()() });

    // Restore feasibility tolerance
    if (inttol > feastol) SetCPLEXParam(CPX_PARAM_EPRHS, feastol);
  }
}
void CplexBackend::AddCPLEXMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> CplexBackend::ConvertCPLEXStatus() {
  namespace sol = mp::sol;
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
}



void CplexBackend::ComputeIIS() {
  int status;
  int cs;
  int nr, nc, nr2, nc2;
  static char* abort_reason[7] = {
  "contradiction", "time limit", "iteration limit", "node limit",
  "objective limit", "memory limit", "user request" };
  
  int nGroups = NumQPCons() + NumSOSCons() + NumIndicatorCons();
  //if (nGroups > 0) // use ext routines if we have non linear constraints
  if(true)
  {
    
    // Check if we have bounds on variables, in which case, consider them as possible
    // conflicts, otherwise don't.
    std::vector<double> lbs(NumVars()), ubs(NumVars());
    CPXgetlb(env(), lp(), lbs.data(), 0, NumVars() - 1);
    CPXgetub(env(), lp(), ubs.data(), 0, NumVars() - 1);
    for (int i = 0; i < NumVars(); i++) {
      if (lbs[i] > MinusInfinity()) nGroups++;
      if (ubs[i] < Infinity()) nGroups++;
    }
    // Add all linear constraints
    nGroups += NumLinCons();

    std::vector<double> grpPref(nGroups, 1);
    std::vector<int> grpBeg(nGroups), grpInd(nGroups);
    std::vector<char> grpType(nGroups);
    for (int i = 0; i < nGroups; i++) grpBeg[i] = i; // each entity has its own group
    int j = 0;
    // First variables with bounds
    for (int i = 0; i < NumVars(); i++)
    {
      if (lbs[i] > MinusInfinity()) {
        grpInd[j] = i;
        grpType[j++] = CPX_CON_LOWER_BOUND;
      }
      if (ubs[i] < Infinity()) {
        grpInd[j] = i;
        grpType[j++] = CPX_CON_UPPER_BOUND;
      }
    }
    // Then all linear constraints
    for (int i = 0; i < NumLinCons(); i++) {
      grpInd[j] = i;
      grpType[j++] = CPX_CON_LINEAR;
    }
    // Then all quadratic constraints
    for (int i = 0; i < NumQPCons(); i++) {
      grpInd[j] = i;
      grpType[j++] = CPX_CON_QUADRATIC;
    }
    // Then all indicator constraints
    for (int i = 0; i < NumIndicatorCons(); i++) {
      grpInd[j] = i;
      grpType[j++] = CPX_CON_INDICATOR;
    }
    // Then all SOS  constraints
    for (int i = 0; i < NumSOSCons(); i++) {
      grpInd[j] = i;
      grpType[j++] = CPX_CON_SOS;
    }
    
    // Calculate information on the above defined groups
    CPXrefineconflictext(env(), lp(), nGroups, nGroups,
      grpPref.data(), grpBeg.data(), grpInd.data(), grpType.data());
    
    // Get calculated information
    std::vector<int> grpStat(nGroups);
    CPXgetconflictext(env(), lp(), grpStat.data(), 0, nGroups - 1);
    for (int i = j = 0; i < NumVars(); i++)
    {
      if (lbs[i] > MinusInfinity()) {
        iisColIndices.push_back(i);
        iisColValues.push_back(grpStat[j++]);
      }
      if (ubs[i] < Infinity()) {
        iisColIndices.push_back(i);
        iisColValues.push_back(grpStat[j++]);
      }
    }
    for (int i = 0; i < NumLinCons(); i++) {
      iisRowIndices.push_back(i);
      iisRowValues.push_back(grpStat[j++]);
    }
    for (int i = 0; i < NumQPCons(); i++) {
      iisRowIndices.push_back(i);
      iisRowValues.push_back(grpStat[j++]);
    }
    for (int i = 0; i < NumIndicatorCons(); i++) {
      iisRowIndices.push_back(i);
      iisRowValues.push_back(grpStat[j++]);
    }
    for (int i = 0; i < NumSOSCons(); i++) {
      iisRowIndices.push_back(i);
      iisRowValues.push_back(grpStat[j++]);
    }
  }
  else {
    CPLEX_CALL(CPXrefineconflict(env(), lp(), &nr, &nc));
    iisRowIndices.resize(nr);
    iisRowValues.resize(nr);
    iisColIndices.resize(nc);
    iisColValues.resize(nc);
    status = CPXgetconflict(env(), lp(), &cs,
      iisRowIndices.data(), iisRowValues.data(), &nr2,
      iisColIndices.data(), iisColValues.data(), &nc2);

    if (cs == CPX_STAT_CONFLICT_FEASIBLE) {
      fmt::print("No IIS after all: problem is feasible!");
      return;
    }
    if (cs != CPX_STAT_CONFLICT_MINIMAL) {
      if ((cs - CPX_STAT_CONFLICT_ABORT_CONTRADICTION) < 0
        || cs > CPX_STAT_CONFLICT_ABORT_CONTRADICTION + 6)
        fmt::print("Surprise conflict status = {} from CPXgetconflict\n", cs);
      else
        fmt::print("Search for conflicts aborted because of {}",
          abort_reason[cs - CPX_STAT_CONFLICT_ABORT_CONTRADICTION]);
      return;
    }
    if (nr2 > nr || nc2 > nc) {
      fmt::print("Surprise confnumrows = {} (should be <= {}), "
        "\nconfnumcols = {} (should be <= {}) from CPXgetconflict.",
        nr2, nr, nc2, nc);
      return;
    }
  }
  
}

int IISCplexToAMPL(int i) {
  static int stmap[7] = { 0, 5, 6, 7, 4, 1, 3 };
  i++;
  if ((i < 0) || (i > 6))
    return 8;
  return stmap[i];
}
IIS CplexBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}
ArrayRef<int> CplexBackend::VarsIIS() {

  std::vector<int> iis(NumVars(), 0);
  for (int i = 0; i < iisColIndices.size(); i++)
    iis[iisColIndices[i]] = (int)IISCplexToAMPL(iisColValues[i]);
  return iis;
}
pre::ValueMapInt CplexBackend::ConsIIS() {
  std::vector<int> iis_lincon(NumLinCons(), 0), iis_qc(NumQPCons(), 0),
    iis_indcon(NumIndicatorCons(), 0), iis_soscon(NumSOSCons(), 0);

  for (int i = 0; i < iisRowIndices.size(); i++)
    iis_lincon[iisRowIndices[i]] = (int)IISCplexToAMPL(iisRowValues[i]);
  int j = NumLinCons();
  for (int i = 0; i < NumQPCons(); i++)
    iis_qc[i] = (int)IISCplexToAMPL(iisRowValues[j++]);
  for (int i = 0; i < NumIndicatorCons(); i++)
    iis_indcon[i] = (int)IISCplexToAMPL(iisRowValues[j++]);
  for (int i = 0; i < NumSOSCons(); i++)
    iis_soscon[i] = (int)IISCplexToAMPL(iisRowValues[j++]);
  return { {{ CG_Linear, iis_lincon },
        { CG_Quadratic, iis_qc },
      { CG_SOS, iis_soscon },
      { CG_General, iis_indcon }} };
}

void CplexBackend::FinishOptionParsing() {
  // Set output on screen
  int lp, mip, bar;
  GetSolverOption(CPX_PARAM_SIMDISPLAY, lp);
  GetSolverOption(CPX_PARAM_MIPDISPLAY, mip);
  GetSolverOption(CPX_PARAM_BARDISPLAY, bar);
  assert(storedOptions_.outlev_ <= 2);
  int olp[] = { 0, 1, 2 };
  int omip[] = { 0, 3, 5 };
  lp = lp ? lp : olp[storedOptions_.outlev_];
  mip = mip ? mip : omip[storedOptions_.outlev_];
  bar = bar ? bar : olp[storedOptions_.outlev_];
  SetSolverOption(CPX_PARAM_SIMDISPLAY, lp);
  SetSolverOption(CPX_PARAM_MIPDISPLAY, mip);
  SetSolverOption(CPX_PARAM_BARDISPLAY, bar);
  if (!storedOptions_.logFile_.empty())
  {
    if (lp < 1) SetSolverOption(CPX_PARAM_SIMDISPLAY, 1);
    if (mip < 1) SetSolverOption(CPX_PARAM_MIPDISPLAY, 1);
    CPLEX_CALL(CPXsetlogfilename(env(), storedOptions_.logFile_.data(), "w"));
  }
  set_verbose_mode(storedOptions_.outlev_ > 0);

  // Set behaviour for solultion pool related options
  if (!need_multiple_solutions()) {
    storedOptions_.populate_ = 0;
    storedOptions_.poolIntensity_ = 0;
  }
  else {
    int poolIntensity = 0, populate = 0;
    switch (storedOptions_.nPoolMode_) {
    case 0: poolIntensity = 0; populate = 0; break;
    case 1: poolIntensity = 2; populate = 1; break;
    case 2: poolIntensity = 4; populate = 1; break;
    }
    // Override the below only if not set
    if (storedOptions_.populate_ < 0) storedOptions_.populate_ = populate;
    if (storedOptions_.poolIntensity_ < 0) storedOptions_.poolIntensity_ = poolIntensity;
  }
  SetSolverOption(CPX_PARAM_SOLNPOOLINTENSITY, storedOptions_.poolIntensity_);

}

static const mp::OptionValueInfo lpmethod_values_[] = {
  { "choose", "Automatic (default)", -1},
  { "simplex", "Simplex", 1},
  { "ipm", "Interior point method", 2},
};
static const mp::OptionValueInfo bardisplay_values_[] = {
  { "0", "no information (default)", 0},
  { "1", "balanced setup and iteration information", 1},
  { "2", "diagnostic information", 2}
};
static const mp::OptionValueInfo display_values_[] = {
  { "0", "never (default)", 0},
  { "1", "each factorization", 1},
  { "2", "each iteration", 2}
};
static const mp::OptionValueInfo mipdisplay_values_[] = {
  { "0", "no node log displayed (default)", 0},
  { "1", "each integer feasible solution", 1},
  { "2", "every \"mipinterval\" nodes", 2},
  { "3", "same as 2 plus cutting planes info and info about new incumbents found through MIP starts", 3},
  { "4", "same as 3, plus LP root relaxation info (according to \"display\")"},
  { "5", "same as 4, plus LP subproblems (according to \"display\")"}
};
static const mp::OptionValueInfo mipinterval_values_[] = {
  { "0", "automatic (default)", 0},
  { "n > 0", "every n nodes and every incumbent", 1},
  { "n < 0", "new incumbents and less info the more negative n is", 2}
};
static const mp::OptionValueInfo outlev_values_[] = {
  { "0", "no output (default)", 0},
  { "1", "equivalent to \"bardisplay\"=1, \"display\"=1, \"mipdisplay\"=3", 1},
  { "2", "equivalent to \"bardisplay\"=2, \"display\"=2, \"mipdisplay\"=5", 2}
};

static const mp::OptionValueInfo values_pool_mode[] = {
  { "0", "Just collect solutions during normal solve", 0},
  { "1", "Make some effort at finding additional solutions => poolintensity=2, populate=1" , 1},
  { "2", "Seek \"sol:poollimit\" best solutions (default) => poolintensity=2, populate=1", 2}
};

static const mp::OptionValueInfo values_populate[] = {
  { "0", "no; just keep solutions found during the initial solve", 0},
  { "1", "run \"populate\" after finding a MIP solution" , 1},
};

static const mp::OptionValueInfo values_poolintensity[] = {
  { "0", "Treated as 1 if poolstub is specified without populate, or 2 if populate is specified (default)", 0},
  { "3", "More additions to the solution pool" , 1},
  { "4", "Tries to generate all MIP solutions and keep the best \"sol:poollimit\" ones.", 2}
};

static const mp::OptionValueInfo values_poolreplace[] = {
  { "0", "FIFO (first-in, first-out); default", 0},
  { "1", "Keep best solutions" , 1},
  { "2", "Keep most diverse solutions", 2}
};
////////////////////////////// OPTIONS /////////////////////////////////

void CplexBackend::InitCustomOptions() {

  set_option_header(
      "IBM ILOG CPLEX Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cplex_options``. For example::\n"
      "\n"
      "  ampl: option cplex_options 'mipgap=1e-6';\n");
  
  // Solution pool controls
  AddSolverOption("sol:poolgap ams_eps poolgap",
    "Relative tolerance for reporting alternate MIP solutions "
    "(default: 1e75).", CPX_PARAM_SOLNPOOLGAP, 0.0, DBL_MAX);
  AddSolverOption("sol:poolgapabs ams_epsabs poolagap",
    "Absolute tolerance for reporting alternate MIP solutions "
    "(default: 1e75).",CPX_PARAM_SOLNPOOLAGAP, 0.0, DBL_MAX);

  AddStoredOption("sol:poolpopulate populate",
    "Whether to run CPLEX's \"populate\" algorithm in an "
    "attempt to add more solutions to the MIP solution pool:\n"
    "\n.. value-table::\n",
    storedOptions_.populate_, values_populate);

  AddStoredOption("sol:poolintensity poolintensity",
    "How hard to try adding MIP solutions to the solution\n\
		pool.  Useful only if poolstub is specified.\n"
    "\n.. value-table::\n",
    storedOptions_.poolIntensity_, values_poolintensity);

  AddStoredOption("sol:poolmode ams_mode poolmode",
    "Search mode for MIP solutions when sol:stub/sol:count are specified "
    "to request finding several alternative solutions. Overriden by sol:populate and"
    "sol:poolintensity. Values:\n"
    "\n.. value-table::\n",
    storedOptions_.nPoolMode_, values_pool_mode);
  AddOptionSynonyms_Inline_Front("ams_stub", "sol:stub");

  AddSolverOption("sol:poollimit ams_limit poolcapacity poollimit solnlimit",
    "Limit on the number of alternate MIP solutions written. Default: 2100000000.",
    CPX_PARAM_SOLNPOOLCAPACITY, 1, 2100000000);

  AddSolverOption("sol:poolpopulatelim populatelim",
    "Limit on number of solutions added to the solution pool by the populate algorithm. "
    "Default: 20.",
    CPX_PARAM_POPULATELIM, 1, 20);

  AddSolverOption("sol:poolreplace poolreplace",
    "Policy for replacing solutions in the solution pool if "
    "more than poolcapacity solutions are generated:\n"
    "\n.. value-table::\n",
    CPX_PARAM_SOLNPOOLREPLACE, values_poolreplace, 0);

  ReplaceOptionDescription("sol:stub",
    "Stub for alternative MIP solutions, written to files with "
    "names obtained by appending \"1.sol\", \"2.sol\", etc., to "
    "<solutionstub>.  The number of such files written is affected "
    "by the keywords poolgap, poolgapabs, poollimit, poolpopulatelim, "
    "poolpopulate, poolintensity and poolmode. "
    "The number of alternative MIP solution files written is "
    "returned in suffix .nsol on the problem.");
  ReplaceOptionDescription("sol:count",
    "0*/1: Whether to count the number of solutions "
    "and return it in the ``.nsol`` problem suffix. "
    "The number and kind of solutions are controlled by the "
    "sol:pool... parameters. Value 1 implied by sol:stub.");
  // end solution pool controls



  AddStoredOption("tech:outlev outlev",
    "Whether to write CPLEX log lines (chatter) to stdout,"
    "for granular control see \"tech:lpdisplay\", \"tech:mipdisplay\", \"tech:bardisplay\"."
    "Values:\n"
    "\n.. value-table::\n",
    storedOptions_.outlev_, outlev_values_);

  AddSolverOption("tech:bardisplay bardisplay",
    "Specifies how much the barrier algorithm chatters:\n"
    "\n.. value-table::\n",
    CPX_PARAM_BARDISPLAY, bardisplay_values_, 0);

  AddSolverOption("tech:lpdisplay display lpdisplay",
    "Frequency of displaying LP progress information:\n"
    "\n.. value-table::\n",
    CPX_PARAM_SIMDISPLAY, display_values_, 0);

  AddSolverOption("tech:mipdisplay mipdisplay",
    "Frequency of displaying branch-and-bound information:\n"
    "\n.. value-table::\n",
    CPX_PARAM_MIPDISPLAY, mipdisplay_values_, 0);

    AddSolverOption("tech:mipinterval mipinterval",
      "Frequency of node logging for \"tech::mipdisplay\" >=2:\n"
      "\n.. value-table::\n",
      CPX_PARAM_MIPINTERVAL, mipinterval_values_, 0);

  AddStoredOption("tech:logfile logfile",
    "Log file name.", storedOptions_.logFile_);

  AddSolverOption("mip:gap mipgap",
      "Relative optimality gap |bestbound-bestinteger|/(1e-10+|bestinteger|).",
      CPXPARAM_MIP_Tolerances_MIPGap, 0.0, 1.0);

  AddSolverOption("tech:threads threads",
      "How many threads to use when using the barrier algorithm\n"
      "or solving MIP problems; default 0 ==> automatic choice.",
      CPXPARAM_Threads, 0, INT_MAX);

  AddSolverOption("lim:time timelim timelimit",
      "limit on solve time (in seconds; default: no limit).",
      CPXPARAM_TimeLimit, 0.0, DBL_MAX);

}


} // namespace mp


AMPLS_MP_Solver* Open_cplexmp(const char* slv_opt, CCallbacks cb = {}) {
  AMPLS_MP_Solver* slv = 
    AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::CplexBackend()},
      slv_opt, cb);
  return slv;
}

void AMPLSClose_cplexmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

CPXLPptr AMPLSGetModel_cplexmp(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CplexBackend*>(AMPLSGetBackend(slv))->lp();
}

CPXENVptr AMPLSGetEnv_cplexmp(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CplexBackend*>(AMPLSGetBackend(slv))->env();
}

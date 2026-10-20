#include <vector>
#include <climits>
#include <cfloat>
#include <stdio.h> // for fflush

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "xpressbackend.h"

extern "C" {
  #include "xpress-ampls-c-api.h"    // Xpressmp AMPLS C API
  #include "xprs_mse_defaulthandler.h"
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptXpressmp(void* prob) {
  XPRSinterrupt((XPRSprob)prob, XPRS_STOP_CTRLC);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateXpressmpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::XpressmpBackend()};
}


namespace mp {
  int XpressmpBackend::outlev_ = 3;
/// Create Xpressmp Model Manager
/// @param gc: the Xpressmp common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return XpressmpModelMgr
std::unique_ptr<BasicModelManager>
CreateXpressmpModelMgr(XpressmpCommon&, Env&, pre::BasicValuePresolver*&);



XpressmpBackend::XpressmpBackend() : msp_(NULL), mse_(NULL) {
  pre::BasicValuePresolver* pPre;
  auto data = CreateXpressmpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

XpressmpBackend::~XpressmpBackend() {
  CloseSolver();
}

void XpressmpBackend::Init(char** argv)  {
  OpenSolver();
  BackendWithModelManager::Init(argv);
}


void XpressmpBackend::OpenSolver() {
  char message[512];
  // Try the registered function first; if not available
  // call the solver's API function directly
  const auto& create_fn = GetCallbacks().init;
  if (create_fn)
    create_fn();
  else 
    int status  = XPRSinit(NULL);
  int status = XPRScreateprob(lp_ref());

  if (status) {
    const auto diag = GetCallbacks().diagnostics;
    if (diag)
      diag();
    XPRSgetlicerrmsg(message, 512);
    MP_RAISE(fmt::format("XPRSinit failed with message:\n{}",
      message));

  }
  model_fixed_ = lp();
  XPRESSMP_CCALL(XPRS_msp_create(&msp_));
  XPRESSMP_CCALL(XPRS_mse_create(&mse_));
  copy_common_info_to_other();
}

void XpressmpBackend::ApplyNativeOptions() {
  for (const auto& entry: storedOptions_.nativeOptions_) {
    auto sep = entry.find_first_of("= ");
    if (std::string::npos == sep)
      MP_RAISE(fmt::format(
          "tech:optionnative: expected \"NAME=value\", got \"{}\"", entry));
    auto value = entry.substr(sep+1);
    auto vstart = value.find_first_not_of(" ");
    ApplyNativeOption(entry.substr(0, sep),
                      std::string::npos == vstart ? "" : value.substr(vstart));
  }
}

void XpressmpBackend::ApplyNativeOption(
    const std::string& name, const std::string& value) {
  int id = 0, type = XPRS_TYPE_NOTDEFINED;
  if (XPRSgetcontrolinfo(lp(), name.c_str(), &id, &type)
      || XPRS_TYPE_NOTDEFINED == type)
    MP_RAISE(fmt::format("Unknown Xpress control \"{}\"", name));
  bool query = ("?" == value);
  auto bad_value = [&]() {
    MP_RAISE(fmt::format("Invalid value \"{}\" for Xpress control \"{}\"",
                         value, name));
  };
  switch (type) {
  case XPRS_TYPE_INT:
  case XPRS_TYPE_INT64: {
    int v = 0;
    if (query) GetSolverOption(id, v);
    else {
      std::size_t pos = 0;
      try { v = std::stoi(value, &pos); } catch (const std::exception&) { bad_value(); }
      if (pos != value.size()) bad_value();
      SetSolverOption(id, v);
    }
    if (query) Print("{}={}\n", name, v);
    break;
  }
  case XPRS_TYPE_DOUBLE: {
    double v = 0.0;
    if (query) GetSolverOption(id, v);
    else {
      std::size_t pos = 0;
      try { v = std::stod(value, &pos); } catch (const std::exception&) { bad_value(); }
      if (pos != value.size()) bad_value();
      SetSolverOption(id, v);
    }
    if (query) Print("{}={}\n", name, v);
    break;
  }
  case XPRS_TYPE_STRING: {
    std::string v;
    if (query) GetSolverOption(id, v);
    else SetSolverOption(id, v = value);
    if (query) Print("{}={}\n", name, v);
    break;
  }
  default:
    MP_RAISE(fmt::format("Xpress control \"{}\" has unsupported type {}",
                         name, type));
  }
}

void XpressmpBackend::ReadBasisFile() {
  const auto& file = storedOptions_.readBasis_;
  if (XPRSreadbasis(lp(), file.c_str(), ""))
    MP_RAISE(fmt::format("Could not read the starting basis from \"{}\":\n{}",
                         file, getErr()));
}

void XpressmpBackend::FinishOptionParsing() {
  ApplyNativeOptions();
  if (storedOptions_.nonconvex_ >= 0)        // deprecated qp:nonconvex
    SetSolverOption(XPRS_IFCHECKCONVEXITY,
                    2 == storedOptions_.nonconvex_ ? 0 : storedOptions_.nonconvex_);
  bool doLog = outlev_ > 0 && outlev_ < 5;
  set_verbose_mode(doLog);
  if (doLog && !msgCallbackSet_)
  {
    msgCallbackSet_ = true;
    XPRSaddcbmessage(lp(), xpdisplay, NULL, 0);
  }

  if (storedOptions_.logFile_.size())
  {
    XPRESSMP_CCALL(XPRSsetlogfile(lp(), storedOptions_.logFile_.c_str()));
    XPRSsetintcontrol(lp(), XPRS_OUTPUTLOG, 1);        /* Chat mode */
  }
}

void XpressmpBackend::CloseSolver() {
  if (mse_) {
    XPRS_mse_destroy(mse_);
    mse_ = nullptr;
  }
  if (msp_) {
    XPRS_msp_destroy(msp_);
    msp_ = nullptr;
  }

  /* Free the fixed model */
  if (lp() != nullptr && lp() != model_fixed_) {
    assert(model_fixed_);
    XPRSdestroyprob(model_fixed_);
  }
  model_fixed_ = nullptr;

  if(lp()!=nullptr)
    XPRESSMP_CCALL(XPRSdestroyprob(lp()));
  *lp_ref() = nullptr;
}

const char* XpressmpBackend::GetBackendName()
  { return "XpressmpBackend"; }

std::string XpressmpBackend::GetSolverVersion() {
  char vbuf[16], rbuf[16];
  int vlen;
  vlen = 0;
  if (XPRSgetstringattrib(lp(), XPRS_XPRESSVERSION, rbuf, (int)sizeof(rbuf), &vlen))
    *rbuf = 0;
  XPRSgetversion(vbuf);
  return fmt::format("{} ({})", rbuf, vbuf);
}

int XpressmpBackend::GetGlobalFlag(const SolverOption& ) const
{ return global_flag_; }

void XpressmpBackend::SetGlobalFlag(const SolverOption& , int v) {
  if (!v) {
    global_flag_ = 0;
    this->ParseSolverOptions(
        nullptr, 0, 0,
        "acc:abs=0 acc:acos=0 acc:acosh=0 acc:asin=0 "
        "acc:asinh=0 acc:atan=0 acc:atanh=0 acc:cos=0 "
        "acc:cosh=0 acc:sin=0 acc:sinh=0 acc:tan=0 acc:tanh=0 "
        "acc:div=0 acc:exp=0 acc:expA=0 "
        "acc:log=0 acc:logA=0 acc:powconstexp=0 acc:pow=0"        );
  }
}

bool XpressmpBackend::IsQCP() const {
  return numQuadCons() > 0;
}

ArrayRef<double> XpressmpBackend::PrimalSolution() {
  if (fixed_model_failed_) return mip_primal_before_fixed_;
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
    int solst;
    // After XPRSoptimize():
  error = XPRSgetsolution(lp(), &solst, x.data(), 0, num_vars-1);
  if (error
      || XPRS_SOLAVAILABLE_NOTFOUND==solst
          //   Infeasible - should return:
          // || XPRS_SOLSTATUS_INFEASIBLE==solst
      )    // @todo keep up2date
    x.clear();
  return x;
}

pre::ValueMapDbl XpressmpBackend::DualSolution() {
  return {{ { CG_Algebraic, DualSolution_LP() } }};
}

ArrayRef<double> XpressmpBackend::DualSolution_LP() {
  if (fixed_model_failed_) return {};
  std::vector<double> pi(NumLinCons());
  int solst;
  int opttype;
  XPRSgetintattrib(model_fixed_, XPRS_OPTIMIZETYPEUSED, &opttype);
  // LP only? Use XPRSgetduals() after XPRSoptimize()
  int error
      = XPRSgetduals(model_fixed_,
                     &solst, pi.data(), 0, NumLinCons()-1);
  if (error)
    pi.clear();
  return pi;
}

double XpressmpBackend::ObjectiveValue() const {
  if (fixed_model_failed_ && !mip_objectives_before_fixed_.empty())
    return mip_objectives_before_fixed_.front();
  return getDblAttr(XPRS_OBJVAL);
}

int XpressmpBackend::NodeCount() const {
  return getIntAttr(XPRS_NODES);
}

int XpressmpBackend::SimplexIterations() const {
  return getIntAttr(XPRS_SIMPLEXITER);
}

int XpressmpBackend::BarrierIterations() const {
  return getIntAttr(XPRS_BARITER);
}

int XpressmpBackend::CrossoverIterations() const {
  return getIntAttr(XPRS_CROSSOVERITER);
}

int XpressmpBackend::SLPIterations() const {
  return getIntAttr(XPRS_SLPITER);
}

std::map<std::string, std::variant<int, double, std::string>>
XpressmpBackend::SolutionStats() {
  std::map<std::string, std::variant<int, double, std::string>> stats;
  stats["simplex_iterations"] = SimplexIterations();
  stats["barrier_iterations"] = BarrierIterations();
  stats["crossover_iterations"] = CrossoverIterations();
  stats["SLP_iterations"] = SLPIterations();
  stats["node_count"] = NodeCount();
  return stats;
}

void XpressmpBackend::DoWriteProblem(const std::string& name) {
  char const* wpflags = NULL;
  if (ends_with(name, ".lp"))
    wpflags = "l";
  else if (ends_with(name, ".mps"))
    wpflags = "";
  else
    throw std::runtime_error(
        fmt::format(
            "Expected \"writeprob=...\" to specify a filename ending in \".lp\"\n"
      "or \".mps\"; got \"{}\".\n", name));
  XPRESSMP_CCALL(XPRSwriteprob(lp(), name.c_str(), wpflags));
}

void XpressmpBackend::DoWriteSolution(const std::string& name) {
  char const* wpflags = "v";  // do not add extension
  auto CheckEnds = [&](std::string e) {
    return ends_with(name, e);
  };
  if (CheckEnds("sol") || CheckEnds("prt"))
    XPRESSMP_CCALL(XPRSwriteprtsol(lp(), name.c_str(), wpflags));
  else if (CheckEnds("iis"))
    XPRESSMP_CCALL(XPRSiiswrite(lp(), 1, name.c_str(), 0, "lv"));
  else if (CheckEnds("bas") || CheckEnds("bss"))
    XPRESSMP_CCALL(XPRSwritebasis(lp(), name.c_str(), wpflags));
}


void XpressmpBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptXpressmp, lp());
}

void XpressmpBackend::DoXPRESSTune() {
  SetSolverOption(XPRS_TUNEROUTPUTPATH, tunebase().data());
  if (tunename().size())
    SetSolverOption(XPRS_TUNERSESSIONNAME, tunename().data());

  if (!tunermethodfilewrite().empty())
      XPRESSMP_CCALL(XPRStunerwritemethod(lp(), tunermethodfilewrite().c_str()));

  if (!tunermethodfileread().empty())
      XPRESSMP_CCALL(XPRStunerreadmethod(lp(), tunermethodfileread().c_str()));
  XPRESSMP_CCALL(XPRStune(lp(), ""));
}

void XpressmpBackend::GetSolverOption(MspControl key, int& value) const {
  XPRESSMP_CCALL(XPRS_msp_getintcontrol(msp_, int(key), &value));
}
void XpressmpBackend::SetSolverOption(MspControl key, int value) {
  XPRESSMP_CCALL(XPRS_msp_setintcontrol(msp_, int(key), value));
}
void XpressmpBackend::GetSolverOption(MspControl key, double& value) const {
  XPRESSMP_CCALL(XPRS_msp_getdblcontrol(msp_, int(key), &value));
}
void XpressmpBackend::SetSolverOption(MspControl key, double value) {
  XPRESSMP_CCALL(XPRS_msp_setdblcontrol(msp_, int(key), value));
}
void XpressmpBackend::GetSolverOption(MseControl key, int& value) const {
  XPRESSMP_CCALL(XPRS_mse_getintcontrol(mse_, int(key), &value));
}
void XpressmpBackend::SetSolverOption(MseControl key, int value) {
  XPRESSMP_CCALL(XPRS_mse_setintcontrol(mse_, int(key), value));
}
void XpressmpBackend::GetSolverOption(MseControl key, double& value) const {
  XPRESSMP_CCALL(XPRS_mse_getdblcontrol(mse_, int(key), &value));
}
void XpressmpBackend::SetSolverOption(MseControl key, double value) {
  XPRESSMP_CCALL(XPRS_mse_setdblcontrol(mse_, int(key), value));
}

void XpressmpBackend::Solve() {
  fixed_model_failed_ = false;
  int nsols = 10;
  if (tunebase().size())
    DoXPRESSTune();
  
  auto flags = XPRESSSolveFlags();
  int solvest=0, solst=0;
  if (IsMIP()) {
    if (need_multiple_solutions() || storedOptions_.nbest_ > 1) {
      if (storedOptions_.nbest_ == 0) {
        storedOptions_.nbest_ = 20;
      }
      XPRESSMP_CCALL(XPRS_mse_opt(mse_, lp(),
        msp_, XPRS_mse_defaulthandler, 0,
        &storedOptions_.nbest_));
    }
    else
      XPRESSMP_CCALL(XPRSoptimize(lp(), flags.c_str(), &solvest, &solst));
  }
  else
      XPRESSMP_CCALL(XPRSoptimize(lp(), flags.c_str(), &solvest, &solst));

  WindupXPRESSMPSolve();
}

std::string XpressmpBackend::XPRESSSolveFlags() {
  if (bool(storedOptions_.fBarrier_)
      + bool(storedOptions_.fPrimal_)
      + bool(storedOptions_.fDual_)
      + bool(storedOptions_.fNetwork_) >= 2) {
    AddWarning("Ambiguous LP method",
               "Only one of barrier/primal/dual/network should be specified.");
  }
  if (storedOptions_.fBarrier_)
    return "b";
  if (storedOptions_.fPrimal_)
    return "p";
  if (storedOptions_.fDual_)
    return "d";
  if (storedOptions_.fNetwork_)
    return "n";
  return "";
}

void XpressmpBackend::WindupXPRESSMPSolve() { }

void XpressmpBackend::ReportResults() {
  ReportXPRESSMPResults();
  BaseBackend::ReportResults();
}


void XpressmpBackend::ReportXPRESSMPResults() {
  SetStatus( GetSolveResult() );
  AddXPRESSMPMessages();
  if (need_multiple_solutions())
    ReportXPRESSMPPool();
  if (need_fixed_MIP())
    ConsiderXpressFixedModel();
}

void XpressmpBackend::ConsiderXpressFixedModel() {
  if (!IsMIP())
    return;
  // The native operation modifies the main model. Retain the incumbent for
  // reporting if the supplementary solve reaches a limit.
  mip_primal_before_fixed_ = PrimalSolution();
  mip_objectives_before_fixed_ = GetObjectiveValues();
  std::string msg;
  try {
    if (XPRSfixmipentities(model_fixed_, 1))
      msg = "XPRSfixmipentities()";
    else {
      if (outlev_)
        Print("\n\n     ======= Solving the fixed MIP =======\n\n");
      msg = DoXpressFixedModel();
    }
  } catch (const std::exception& e) {
    msg = e.what();
  }
  if (!msg.empty()) {
    fixed_model_failed_ = true;
    AddToSolverMessage(msg + " failed in DoXpressFixedModel().\n");
  }
}
std::string XpressmpBackend::DoXpressFixedModel()
{
  double old_time_limit;
  int old_iteration_limit, old_barrier_limit;
  XPRESSMP_CCALL(XPRSgetdblcontrol(model_fixed_, XPRS_TIMELIMIT, &old_time_limit));
  XPRESSMP_CCALL(XPRSgetintcontrol(model_fixed_, XPRS_LPITERLIMIT, &old_iteration_limit));
  XPRESSMP_CCALL(XPRSgetintcontrol(model_fixed_, XPRS_BARITERLIMIT, &old_barrier_limit));
  auto restore_limits = [&](void*) {
    XPRSsetdblcontrol(model_fixed_, XPRS_TIMELIMIT, old_time_limit);
    XPRSsetintcontrol(model_fixed_, XPRS_LPITERLIMIT, old_iteration_limit);
    XPRSsetintcontrol(model_fixed_, XPRS_BARITERLIMIT, old_barrier_limit);
  };
  std::unique_ptr<void, decltype(restore_limits)> limits_guard(model_fixed_, restore_limits);
  XPRESSMP_CCALL(XPRSsetdblcontrol(model_fixed_, XPRS_TIMELIMIT,
      fixed_model_time_limit() < 0 ? 1e20 : fixed_model_time_limit()));
  XPRESSMP_CCALL(XPRSsetintcontrol(model_fixed_, XPRS_LPITERLIMIT,
      fixed_model_iter_limit() < 0 ? INT_MAX : fixed_model_iter_limit()));
  XPRESSMP_CCALL(XPRSsetintcontrol(model_fixed_, XPRS_BARITERLIMIT,
      fixed_model_iter_limit() < 0 ? INT_MAX : fixed_model_iter_limit()));
  if (XPRSlpoptimize(model_fixed_, ""))
    return "optimize()";
  int status;
  if (XPRSgetintattrib(model_fixed_, XPRS_LPSTATUS, &status))
    return "getStatus()";
  static const char* statusName[] = {
    "Infeasible",
    "Cut off",
    "Unfinished",
    "Unbounded",
    "Cutoff in dual",
    "Unsolved",
    "Nonconvex" };
  if (status != XPRS_LP_OPTIMAL) {
    if (status >= XPRS_LP_INFEAS && status <= XPRS_LP_NONCONVEX)
      return fmt::format(
        "Fixed model status: {}. XPRSlpoptimize",
        statusName[status - XPRS_LP_INFEAS]);
    else
      return fmt::format(
        "Surprise status {} after XPRSlpoptimize", status);
  }
  int f;
  if (!XPRSgetintattrib(model_fixed_, XPRS_SIMPLEXITER, &f)) {
    if (f)
      AddToSolverMessage(fmt::format(
        "Fixed MIP for mip:basis: {} simplex iteration{}\n",
        f, "s"[f == 1.]));
  }
  return {};

}
  std::vector<double> XpressmpBackend::getPoolSolution(int id)
  {
    std::vector<double> vars(NumVars());
    int j, nx;
    XPRS_msp_getsol(msp_, id, &j, vars.data(), 0, NumVars() - 1, &nx);
    return vars;
  }
  double XpressmpBackend::getPoolObjective(int id)
  {
    double obj;
    int j;
    XPRS_msp_getdblattribprobsol(msp_, lp(), id, &j, XPRS_MSP_SOLPRB_OBJ, &obj);
    return obj;
  }
  void XpressmpBackend::ReportXPRESSMPPool() {
    if (!IsMIP())
      return;
    int iPoolSolution = -1;
    int nPool;
    XPRS_msp_getintattrib(msp_, XPRS_MSP_SOLUTIONS, &nPool);
    auto sid = std::vector<int>(nPool);
    int nret, nsols;
    XPRS_msp_getsollist(msp_, 0, 0, 0, 1, nPool, sid.data(),
      &nret, &nsols);
    nPool = nret;
    int id;
    for (int i = 0; i < nPool; i++) {
      id = sid[i];
      auto mv = GetValuePresolver().PostsolveSolution({
            getPoolSolution(id) });
      ReportIntermediateSolution(
        { mv.GetVarValues()(), mv.GetConValues()(),
          { getPoolObjective(id) } });   // not when multiobj
    }
  }


  void XpressmpBackend::AddXPRESSMPMessages() {
    auto iter = SimplexIterations();    // Always
    AddToSolverMessage(
        fmt::format("{} simplex iterations\n",iter));
    if (auto nbi = BarrierIterations())
      AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
    if (auto nci = CrossoverIterations())
      AddToSolverMessage(
          fmt::format("{} crossover iterations\n", nci));
    if (auto nslpi = SLPIterations())
      AddToSolverMessage(
          fmt::format("{} SLP iterations\n", nslpi));
    if (auto nnd = NodeCount())
      AddToSolverMessage(
        fmt::format("{} branching nodes\n", nnd));
  }

  std::pair<int, std::string> XpressmpBackend::GetSolveResult() {
    auto opttype = getIntAttr(XPRS_OPTIMIZETYPEUSED);

    // Catch some specific codes
    switch (opttype) {
    case XPRS_OPTIMIZETYPE_NONE:
      return GetSolveResult_General();
    case XPRS_OPTIMIZETYPE_LP:
      return GetSolveResult_LP();
    case XPRS_OPTIMIZETYPE_MIP:
      return GetSolveResult_MIP();
    case XPRS_OPTIMIZETYPE_LOCAL:
    case XPRS_OPTIMIZETYPE_GLOBAL:
      return GetSolveResult_NLP();    // Same codes?
    default:
      return { sol::FAILURE, "failure, unknown OPTIMIZETYPE code" };
    }

    return { sol::UNKNOWN, "not solved" };
  }

  std::pair<int, std::string> XpressmpBackend::GetSolveResult_LP() {
    namespace sol = mp::sol;
    auto solstatus = getIntAttr(XPRS_SOLSTATUS);
    bool fFeasible = (XPRS_SOLSTATUS_FEASIBLE==solstatus);

    {
      auto status = getIntAttr(XPRS_LPSTATUS);
      switch (status) {
      case XPRS_LP_OPTIMAL:
        return { sol::SOLVED, "optimal solution" };
      case XPRS_LP_INFEAS:
        return { sol::INFEASIBLE, "infeasible problem" };
      case XPRS_LP_UNBOUNDED:
        if (fFeasible)
          return { sol::UNBOUNDED_FEAS,
                  "unbounded problem, feasible solution returned" };
        return { sol::LIMIT_INF_UNB,
                "unbounded or infeasible, "
                "switch off dual presolve or run IIS Finder" };
      case XPRS_LP_CUTOFF:
        return { sol::UNCERTAIN, "objective cutoff" };
      case XPRS_LP_CUTOFF_IN_DUAL:
        return { sol::UNCERTAIN, "objective cutoff in dual" };
      case XPRS_LP_UNSOLVED:
        return { sol::NUMERIC, "numerical issues" };
      case XPRS_LP_UNSTARTED:
        return { sol::UNKNOWN, "LP unstarted" };
      case XPRS_LP_NONCONVEX:
        return { sol::FAILURE, "problem is nonconvex, consider using FICO Xpress Global" };
      case XPRS_LP_UNFINISHED:
      default:
        return GetSolveResult_General();
      }
    }
    return { sol::UNKNOWN, "unknown LP status" };
  }

  std::pair<int, std::string> XpressmpBackend::GetSolveResult_MIP() {
    namespace sol = mp::sol;
    auto solstatus = getIntAttr(XPRS_SOLSTATUS);
    bool fFeasible = (XPRS_SOLSTATUS_FEASIBLE==solstatus);

    {
      auto status = getIntAttr(XPRS_MIPSTATUS);
      switch (status) {
      case XPRS_MIP_OPTIMAL:
        return { sol::SOLVED, "optimal solution" };
      case XPRS_MIP_SOLUTION:
        return { sol::LIMIT_FEAS, "limit, feasible solution" };
      case XPRS_MIP_INFEAS:
        return { sol::INFEASIBLE, "infeasible problem" };
      case XPRS_MIP_UNBOUNDED:
        if (fFeasible)
          return { sol::UNBOUNDED_FEAS,
                  "unbounded problem, feasible solution returned" };
        return { sol::LIMIT_INF_UNB,
                "unbounded or infeasible, "
                "switch off dual presolve or run IIS Finder" };
      case XPRS_MIP_LP_NOT_OPTIMAL:
      case XPRS_MIP_NO_SOL_FOUND:
      case XPRS_MIP_LP_OPTIMAL:
      default:
        return GetSolveResult_General();
      }
    }
    return { sol::UNKNOWN, "unknown MIP status" };
  }

  std::pair<int, std::string> XpressmpBackend::GetSolveResult_NLP() {
    namespace sol = mp::sol;
    // auto solstatus = getIntAttr(XPRS_SOLSTATUS);
    // bool fFeasible = (XPRS_SOLSTATUS_FEASIBLE==solstatus);

    {
      auto status = getIntAttr(XPRS_NLPSTATUS);
      auto solstatus = getIntAttr(XPRS_NLPSOLSTATUS);
      // should work for UNBOUNDED (@todo also other problem types):
      bool fFeasible = (XPRS_NLPSOLSTATUS_NONE!=solstatus);
      std::string msg_duals =
          (XPRS_NLPSOLSTATUS_SOLUTION_NODUALS==solstatus ||
                        XPRS_NLPSOLSTATUS_GLOBALLYOPTIMAL_NODUALS==solstatus)
                                  ? ", no duals" : ", with duals";

      switch (status) {
      case XPRS_NLPSTATUS_LOCALLY_OPTIMAL:
        return { sol::OPTIMAL_LOCALLY,
                "locally optimal solution" + msg_duals };
      case XPRS_NLPSTATUS_OPTIMAL:
        return { sol::SOLVED, "optimal solution" + msg_duals };
      case XPRS_NLPSTATUS_LOCALLY_INFEASIBLE:
        return { sol::INFEASIBLE_LOCALLY, "locally infeasible" };
      case XPRS_NLPSTATUS_INFEASIBLE:
        return { sol::INFEASIBLE, "infeasible problem" };
      case XPRS_NLPSTATUS_UNBOUNDED:
        if (fFeasible)
          return { sol::UNBOUNDED_FEAS,
              "unbounded problem, feasible solution returned" +
                  msg_duals
          };
        return { sol::LIMIT_INF_UNB,
                "unbounded or infeasible, "
                "switch off dual presolve or run IIS Finder" };
      case XPRS_NLPSTATUS_UNFINISHED:
      case XSLP_NLPSTATUS_UNSOLVED:
      default:
        return GetSolveResult_General();
      }
    }
    return GetSolveResult_General();
  }

  std::pair<int, std::string> XpressmpBackend::GetSolveResult_General() {
    namespace sol = mp::sol;
    auto solvestatus = getIntAttr(XPRS_SOLVESTATUS);
    auto solstatus = getIntAttr(XPRS_SOLSTATUS);
    auto stopstatus = getIntAttr(XPRS_STOPSTATUS);
    bool fFeasible = (XPRS_SOLSTATUS_FEASIBLE==solstatus);

    switch (solstatus) {
    case XPRS_SOLSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case XPRS_SOLSTATUS_FEASIBLE:
      switch (stopstatus) {
      case XPRS_STOP_TIMELIMIT:
        return { sol::LIMIT_FEAS_TIME, "time limit, feasible solution" };
      case XPRS_STOP_WORKLIMIT:
        return { sol::LIMIT_FEAS_WORK, "work limit, feasible solution" };
      case XPRS_STOP_CTRLC:
      case XPRS_STOP_USER:
        return { sol::LIMIT_FEAS_INTERRUPT, "interrupted, feasible solution" };
      case XPRS_STOP_NODELIMIT:
        return { sol::LIMIT_FEAS_NODES, "node limit, feasible solution" };
      case XPRS_STOP_ITERLIMIT:
        return{ sol::LIMIT_FEAS_ITER, "iteration limit, feasible solution" };
      case XPRS_STOP_MIPGAP:
        return{ sol::LIMIT_FEAS_GAP, "MIP gap reached, feasible solution" };
      case XPRS_STOP_SOLLIMIT:
        return { sol::LIMIT_FEAS_NUMSOLS, "solution limit" };
      case XPRS_STOP_MEMORYERROR:
        return { sol::LIMIT_FEAS_SOFTMEM, "memory limit, feasible solution" };
      case XPRS_STOP_NUMERICALERROR:
        return { sol::UNCERTAIN, "numerical error, solution candidate returned" };
      default:
        return { sol::LIMIT_FEAS, "limit, feasible solution" };
      }
    case XPRS_SOLSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case XPRS_SOLSTATUS_UNBOUNDED:
      if (fFeasible)
        return { sol::UNBOUNDED_FEAS,
                "unbounded problem, feasible solution returned" };
      return { sol::LIMIT_INF_UNB,
              "unbounded or infeasible, "
              "switch off dual presolve or run IIS Finder" };
    case XPRS_SOLSTATUS_NOTFOUND:
      if (XPRS_SOLVESTATUS_FAILED==solvestatus)
        return { sol::FAILURE, "failure, no solution" };
      switch (stopstatus) {
      case XPRS_STOP_TIMELIMIT:
        return { sol::LIMIT_NO_FEAS_TIME, "time limit, without a feasible solution" };
      case XPRS_STOP_WORKLIMIT:
        return { sol::LIMIT_NO_FEAS_WORK, "work limit, without a feasible solution" };
      case XPRS_STOP_CTRLC:
      case XPRS_STOP_USER:
        return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted, without a feasible solution" };
      case XPRS_STOP_NODELIMIT:
        return { sol::LIMIT_NO_FEAS_NODES, "node limit, without a feasible solution" };
      case XPRS_STOP_ITERLIMIT:
        return{ sol::LIMIT_NO_FEAS_ITER, "iteration limit, without a feasible solution" };
      case XPRS_STOP_MEMORYERROR:
        return { sol::LIMIT_NO_FEAS_SOFTMEM, "memory limit, without a feasible solution" };
      case XPRS_STOP_NUMERICALERROR:
        return { sol::UNCERTAIN, "numerical error" };
      default:
        return { sol::LIMIT_NO_FEAS, "limit, without a feasible solution" };
      }
    }
    return { sol::UNKNOWN, "unknown status" };
  }


  void XpressmpBackend::CreateSolutionPoolEnvironment() {
    XPRESSMP_CCALL(XPRS_msp_probattach(msp_, lp()));
    if (outlev_ > 0)
      XPRESSMP_CCALL(XPRS_mse_addcbmsghandler(mse_, xp_mse_display, NULL, 0));

    //SetSolverOption(XPRS_HEURSTRATEGY, 0);
    SetSolverOption(XPRS_MIPDUALREDUCTIONS, 2);

    if (storedOptions_.pooldualred_ != 2 || storedOptions_.pooldupcol_ != 2)
    {
      // Adapt presolve operations to solution pool specific controls
      int i, j = 0, k = 0;
      GetSolverOption(XPRS_PRESOLVEOPS, i);
      if (storedOptions_.pooldualred_ != 2) {
        j = 1 << 3;
        if (storedOptions_.pooldualred_ == 1)
          k = j;
      }
      if (storedOptions_.pooldupcol_ != 2) {
        j |= 1 << 5;
        if (storedOptions_.pooldupcol_ == 1)
          k |= 1 << 5;
      }
      i = (i & ~j) | k;
      SetSolverOption(XPRS_PRESOLVEOPS, i);
    }
  }

  void XpressmpBackend::InputExtras() {
    BaseBackend::InputExtras();
    InputXPRESSExtras();
  }
  void XpressmpBackend::InputXPRESSExtras() {
    if (need_multiple_solutions())
      CreateSolutionPoolEnvironment();
    if (!storedOptions_.readBasis_.empty())
      ReadBasisFile();
  }

  ArrayRef<double> XpressmpBackend::GetObjectiveValues()
  {
    if (fixed_model_failed_) return mip_objectives_before_fixed_;
    int n = NumObjs();
    if (n==0)
      return { };
    if (n==1)
      return std::vector<double>{ObjectiveValue()};
    std::vector<double> objs(n,
      std::numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < n; i++)
      if (XPRScalcobjn(lp(), i, NULL, &objs[i]))
        return {};
    return objs;
  }
  
  void XpressmpBackend::ObjPriorities(ArrayRef<int> priority) {
    for (int i = 0; i < (int)priority.size(); ++i) 
      XPRESSMP_CCALL(XPRSsetobjintcontrol(lp(), i, 
        XPRS_OBJECTIVE_PRIORITY, priority[i]));
  }

  void XpressmpBackend::ObjWeights(ArrayRef<double> val) {
    for (int i = 0; i < (int)val.size(); ++i) 
      XPRESSMP_CCALL(XPRSsetobjdblcontrol(lp(), i, 
        XPRS_OBJECTIVE_WEIGHT, val[i]));
    }

  void XpressmpBackend::ObjAbsTol(ArrayRef<double> val) {
    for (int i = 0; i < (int)val.size(); ++i)
      XPRESSMP_CCALL(XPRSsetobjdblcontrol(lp(), i,
        XPRS_OBJECTIVE_ABSTOL, val[i]));
  }

  void XpressmpBackend::ObjRelTol(ArrayRef<double> val) {
    for (int i = 0; i < (int)val.size(); ++i)
      XPRESSMP_CCALL(XPRSsetobjdblcontrol(lp(), i,
        XPRS_OBJECTIVE_RELTOL, val[i]));
  }

  ////////////////////////////// OPTIONS /////////////////////////////////

  static const mp::OptionValueInfo pool_values_[] = {
    { "0", "Yes (default, can be expensive)", 0},
    { "1", "No", 1},
    { "2", "Honor presolveops bit 3 (2^3 = 8)", 2}
  };

  static const mp::OptionValueInfo nlp_presolve_values_[] = {
      { "0", "No", 0},
      { "1", "Yes (default)", 1},
      { "2", "Low memory presolve. Original problem is not restored "
            "by postsolve and dual solution may not be completely "
            "postsolved.", 2}
  };

  static const mp::OptionValueInfo nlp_solver_values_[] = {
      { "-1", "If the license allows and there are no user functions "
             "or multistart jobs, FICO Xpress Global will be called, "
             "otherwise a local solver (default)", -1},
      { "1", "The algorithm selected by alg:localsolver will be used to "
            "find a locally optimal solution", 1},
      { "2", "FICO Xpress Global will be used to find a globally "
            "optimal solution.", 2}
  };

  static const mp::OptionValueInfo nlp_localsolver_values_[] = {
      { "-1", "Automatic selection, based on model characteristics "
             "and solver availability (default)", -1 },
      { "0", "Use Xpress-SLP (always available)", 0},
      { "1", "Use Knitro if available", 1},
      { "2", "Use Xpress-Optimizer if possible (convex quadratic "
            "problems only).", 2}
  };

  static const mp::OptionValueInfo presolveops_values_[] = {
    { "1 = 2^0", "Remove singleton columns", XPRS_PRESOLVEOPS_SINGLETONCOLUMNREMOVAL},
    { "2 = 2^1", "Remove singleton constraints (rows)", XPRS_PRESOLVEOPS_SINGLETONROWREMOVAL},
    { "4 = 2^2", "Forcing row removal", XPRS_PRESOLVEOPS_FORCINGROWREMOVAL},
    { "8 = 2^3", "Dual reductions", XPRS_PRESOLVEOPS_DUALREDUCTIONS},
    { "16 = 2^4", "Redundant row removal", XPRS_PRESOLVEOPS_REDUNDANTROWREMOVAL},
    { "32 = 2^5", "Duplicate column removal", XPRS_PRESOLVEOPS_DUPLICATECOLUMNREMOVAL},
    { "64 = 2^6", "Duplicate row removal", XPRS_PRESOLVEOPS_DUPLICATEROWREMOVAL},
    { "128 = 2^7", "Strong dual reductions", XPRS_PRESOLVEOPS_STRONGDUALREDUCTIONS},
    { "256 = 2^8", "Variable eliminations", XPRS_PRESOLVEOPS_VARIABLEELIMINATIONS},
    { "512 = 2^9", "No IP reductions", XPRS_PRESOLVEOPS_NOIPREDUCTIONS},
    { "1024 = 2^10", "No semi-continuous variable detection", XPRS_PRESOLVEOPS_NOGLOBALDOMAINCHANGE},
    { "2048 = 2^11", "No advanced IP reductions", XPRS_PRESOLVEOPS_NOADVANCEDIPREDUCTIONS},
    { "4096 = 2^12", "No eliminations on integers", 12},
    { "16384 = 2^14", "Linearly dependant row removal", XPRS_PRESOLVEOPS_LINEARLYDEPENDANTROWREMOVAL},
      { "32768 = 2^15", "No integer variable and SOS detection", XPRS_PRESOLVEOPS_NOINTEGERVARIABLEANDSOSDETECTION},
      { "65536 = 2^16", "No implied bounds", XPRS_PRESOLVEOPS_NOIMPLIEDBOUNDS},
      { "131072 = 2^17", "No clique presolve", XPRS_PRESOLVEOPS_NOCLIQUEPRESOLVE},
      { "262144 = 2^18", "No mod2 presolve", XPRS_PRESOLVEOPS_NOMOD2REDUCTIONS},
    { "536870912 = 2^29", "No dual reduction on globals", XPRS_PRESOLVEOPS_NODUALREDONGLOBALS},
  };

  static const mp::OptionValueInfo pooldups_values_[] = {
    { "0", "Retain all duplicates", 0},
    { "1", "Discard exact matches", 1},
    { "2", "Discard exact matches of continuous variables "
          "and matches of rounded values of discrete variables", 2},
    { "3", "Discard matches of rounded values of discrete "
    " variables (default)", 3}
  };





















  




#define xstr(a) stringify(a)
#define stringify(a) #a




















  static const mp::OptionValueInfo values_miqcpalg[]{
    {"-1", "automatic (default)", -1},
    {"0", "barrier", 0},
    {"1", "outer approximations", 1}
  };
  static const mp::OptionValueInfo values_netstalllimit[]{
  {"-1", "automatic (default)", -1},
  {"0", "no limit", 0},
  {"n > 0", "limit to n network simplex iterations", 1}
  };

  static const mp::OptionValueInfo values_nodeselection[]{
    {"1", "local first: choose between descendant and sibling nodes if available; choose from all outstanding nodes otherwise", 1},
    {"2", "best first: choose from all outstanding nodes", 2},
    {"3", "local depth first: choose between descendant and sibling nodes if available; choose from the deepest nodes otherwise", 3},
    {"4", "best first, then local first: best first is used for the first BREADTHFIRST nodes, after which local first is used", 4},
    {"5", "pure depth first: choose from the deepest outstanding nodes", 5}
  };
  static const mp::OptionValueInfo values_numericalemphasis[]{
        {"-1", "automatic (default)", -1},
        {"0", "epmhasize speed", 0},
        {"1", "mild emphasis on numerical stability", 1},
        {"2", "medium emphasis on numerical stability", 2},
        {"3", "strong emphasis on numerical stability", 3},
  };

  static const mp::OptionValueInfo values_outlev[]{
  {"0","none", 0},
  {"1", "all", 1},
  {"2", "information", 2},
  {"3", "warnings & errors only (default)", 3},
  {"4", "errors", 4},
  {"5","none", 5}
  };
  static const mp::OptionValueInfo values_precoefelim[]{
  {"0","disabled", 0},
  {"1", "remove as many coefficients as possible", 1},
  {"2", "cautious eliminations", 2}
};
  static const mp::OptionValueInfo values_preconvertseparable[]{
 {"-1","automatic (default)", -1},
    {"0", "disable", 0},
 {"1", "enable reformulation to diagonal quadratic constraints.", 1},
 {"2", "1, plus reduction to second-order cones", 2},
 {"3", "2, plus the objective function is converted to a constraint and treated as a quadratic constraint", 3}
 };

  static const mp::OptionValueInfo values_predomcol[]{
    {"-1","automatic (default)", -1},
    {"0", "disable", 0},
    {"1", "cautious", 1},
    {"2", "aggressive: all candidate will be checked", 2},
    {"3", "includes 1 and 2 but also looks for more generic column domination", 3}
  };


  static const mp::OptionValueInfo values_preduprow[]{
    {"-1","automatic (default)", -1},
    {"0", "disable", 0},
    {"1", "eliminate only rows that are identical in all variables", 1},
    {"2", "1 plus eliminate duplicate rows with simple penalty variable expressions", 2},
    {"3", "2 plus eliminate duplicate rows with more complex penalty variable expressions", 3}
  };
  
  static const mp::OptionValueInfo values_prepermute[]{
    {"1", "permute rows", 1},
    {"2", "permute columns", 2},
    {"4", "permute global information (for MIP)", 4}
  };

  static const mp::OptionValueInfo values_pricingalg[]{
   {"-1", "partial pricing", 1},
   {"0", "automatic choice (default)", 0},
   {"1", "devex pricing", 1},
   {"2", "steepest edge", 2},
    {"3", "steepest edge with initial weights", 3}
  };












  const mp::OptionValueInfo values_sleeponthreadwait[] = {
   {"-1", "automatically determined", -1},
   {"0", "no (busy-wait)", 0},
   {">0", "yes (sleep, might add overhead)", 1}
  };


 
  

/// @note The AddSolverOptions()
///   are complemented by those from
///   optprm.h and nlprm.h,
///   obtained from optimizerparams.xml and nonlinearparams.xml
///   using parse-xpress-options
///   (see BUILD_XPRESS_OPTION_PARSER).
///
/// @note Might add prefix XPRS_ to the params of optimizerparams.xml
///   as the last argument of parse_xpress_options.
///
/// @note Run the parser as follows:
///   bin/parse-xpress-options ../../solver-libraries/xpress/46.01.01/xml/nonlinearparams.xml ../mp/solvers/xpress/nlprm.h Nonlinear
///   bin/parse-xpress-options ../../solver-libraries/xpress/46.01.01/xml/optimizerparams.xml ../mp/solvers/xpress/optprm.h Optimizer XPRS
///
/// @note Option names are synonymized,
///   descriptions are replaced.
void XpressmpBackend::InitCustomOptions() {

  set_option_header(
      "XPRESS Optimizer Options for AMPL\n"
      "---------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``xpress_options``. For example::\n"
      "\n"
      "  ampl: option xpress_options 'mipgap=1e-6';\n");

  AddToOptionDescription(
      "tech:writesolution",
      "Supported name extensions: "
      "``.sol/.prt``, ``.iss``, ``.bas/.bss``.");


  // ****************************
  // General
  // ****************************
  AddStoredOption("tech:outlev outlev",
    "Whether to write xpress log lines (chatter) to stdout and to file:\n"
    "\n.. value-table::\n",
    outlev_, values_outlev);

  AddSolverOption("tech:cputime",
    "",
    XPRS_CPUTIME, INT_MIN, INT_MAX);

  AddSolverOption("tech:threads",
    "",
    XPRS_THREADS, -1, INT_MAX);

  AddSolverOption("tech:backgroundthreads backgroundthreads",
    "",
    XPRS_BACKGROUNDMAXTHREADS, -1, INT_MAX);

  AddSolverOption("tech:backgroundselect",
    "",
    XPRS_BACKGROUNDSELECT, -1, 1);
  
  AddSolverOption("lim:time timelim",
    "",
    XPRS_TIMELIMIT, 0.0, Infinity());


  AddSolverOption("lim:prerootwork prerootworklim",
    "",
    XPRS_PREROOTWORKLIMIT, -DBL_MAX, DBL_MAX);
  
    AddSolverOption("lim:prerooteffort mip:prerooteffort",
      "",
      XPRS_PREROOTEFFORT, -DBL_MAX, DBL_MAX);

  AddSolverOption("lim:work worklim",
    "",
    XPRS_WORKLIMIT, 0.0, DBL_MAX);


  AddSolverOption("lim:soltime soltimelim",
    "",
    XPRS_SOLTIMELIMIT, 0.0, Infinity());

  AddSolverOption("tech:sleeponthreadwait sleeponthreadwait",
    "Deprecated: ignored by the current Xpress library.",
    XPRS_SLEEPONTHREADWAIT, values_sleeponthreadwait,  -1);

  AddSolverOption("lim:iter iterlim",
    "",
    XPRS_LPITERLIMIT, 2147483645, INT_MAX);

    AddSolverOption("lim:lprefineiter",
      "",
      XPRS_LPREFINEITERLIMIT, -1, INT_MAX);

    AddSolverOption("lim:maxcuttime",
      "",
      XPRS_MAXCUTTIME, 0.0, Infinity());

  AddSolverOption("lim:maxtreefilesize globalfilemax tech:globalfilemax",
    "",
    XPRS_MAXTREEFILESIZE, INT_MIN, INT_MAX);

  AddSolverOption("mip:treefileloginterval globalfileloginterval tech:globalfileloginterval",
    "",
    XPRS_TREEFILELOGINTERVAL, INT_MIN, INT_MAX);

  AddStoredOption("tech:logfile logfile",
    "Log file name; default=no file", storedOptions_.logFile_);

  AddSolverOption("lim:mem memlimit",
    "",
    XPRS_MAXMEMORYHARD, 0, INT_MAX);

    AddSolverOption("lim:softmem softmemlimit",
      "",
      XPRS_MAXMEMORYSOFT, 0, INT_MAX);
  // ****************************
  // Solution pool params
  // ****************************
  AddStoredOption("sol:pooldualred pooldualred",
		"Whether to suppress removal of dominated solutions (via "
      "\"dual reductions\") when poolstub is specified:\n"
    "\n.. value-table::\n",
    storedOptions_.pooldualred_, pool_values_);

  AddStoredOption("sol:pooldupcol pooldupcol",
    "Whether to suppress duplicate variable removal when "
    "poolstub is specified:\n"
    "\n.. value-table::\n",
    storedOptions_.pooldupcol_, pool_values_);

  AddSolverOption("sol:pooldups pooldups",
    "How the solution pool handles duplicate solutions:\n"
    "\n.. value-table::\n",
    MspControl{XPRS_MSP_DUPLICATESOLUTIONSPOLICY}, pooldups_values_, 3);

  AddOptionSynonyms_Inline_Front("ams_stub", "sol:stub");

	AddStoredOption("sol:poollimit poollimit poolnbest",
    "When poollimit = n > 1, the "
    "solution pool (see sol:stub) is allowed to keep n best solutions. "
                  "Default 10.",
    storedOptions_.nbest_);

    AddSolverOption("lim:mipsol",
      "",
      XPRS_MAXMIPSOL, 0, INT_MAX);

    AddSolverOption("lim:nodes nodelim nodelimit",
      "",
      XPRS_MAXNODE, 2147483647, INT_MAX);

    AddSolverOption("lim:stalltime",
      "",
      XPRS_MAXSTALLTIME, 0.0, Infinity());
  // ****************************
  // Generic algorithm controls
  // ****************************
  AddSolverOption("alg:method method lpmethod",
    "",
    XPRS_DEFAULTALG, INT_MIN, INT_MAX);

  AddStoredOption("alg:barrier barrier",
                  "Solve (MIP node) LPs by barrier method.",
                  storedOptions_.fBarrier_);

	AddStoredOption("alg:primal alg:primalopt primal primalopt",
                  "Solve (MIP node) LPs by primal simplex method.",
                  storedOptions_.fPrimal_);

	AddStoredOption("alg:dual alg:dualopt dual dualopt",
                  "Solve (MIP node) LPs by dual simplex method.",
                  storedOptions_.fDual_);

  AddStoredOption("alg:network network",
                  "Solve (substructure of) (MIP node) LPs "
                  "by network simplex method.",
                  storedOptions_.fNetwork_);

  AddSolverOption("alg:clamping",
    "",
    XPRS_CLAMPING, INT_MIN, INT_MAX);

  AddSolverOption("alg:feastol",
    "",
    XPRS_FEASTOL, 0.0, DBL_MAX);

  AddSolverOption("alg:feastolperturb",
    "",
    XPRS_FEASTOLPERTURB, 0.0, DBL_MAX);

  AddSolverOption("alg:feastoltarget",
    "",
    XPRS_FEASTOLTARGET, 0.0, DBL_MAX);

    AddSolverOption("pre:indlinsmallbigm alg:indlinbigm",
      "",
      XPRS_INDLINBIGM, 0.0, DBL_MAX);

    AddSolverOption("lp:folding alg:lpfolding",
      "",
      XPRS_LPFOLDING, INT_MIN, INT_MAX);

    AddSolverOption("iis:max iismax",
      "",
      XPRS_MAXIIS, -1, INT_MAX);
    AddSolverOption("iis:log",
      "",
      XPRS_IISLOG, -1, INT_MAX);
    AddSolverOption("iis:ops",
      "",
      XPRS_IISOPS, -1, INT_MAX);

    AddSolverOption("alg:zerotol",
      "",
      XPRS_MATRIXTOL, 0.0, DBL_MAX);

    AddSolverOption("lp:pivtol pivtol",
      "",
      XPRS_MARKOWITZTOL, 0.0, DBL_MAX);

    AddSolverOption("alg:cutoff cutoff",
      "",
      XPRS_MIPABSCUTOFF, MinusInfinity(), Infinity());

    AddSolverOption("alg:addcutoff addcutoff",
      "",
      XPRS_MIPADDCUTOFF, -1e-10, DBL_MAX);

    AddSolverOption("alg:relcutoff relcutoff",
      "",
      XPRS_MIPRELCUTOFF, 1e-4, Infinity());

    AddSolverOption("tech:seed seed alg:randomseed",
      "",
      XPRS_RANDOMSEED, -INT_MAX, INT_MAX);

    AddSolverOption("alg:refactor",
      "",
      XPRS_REFACTOR, INT_MIN, INT_MAX);

    AddSolverOption("alg:refineops",
      "",
      XPRS_REFINEOPS, INT_MIN, INT_MAX);

    AddSolverOption("alg:resourcestrategy",
      "",
      XPRS_RESOURCESTRATEGY, INT_MIN, INT_MAX);

    AddIntOption("alg:global global",
                 "0/1*: Allow global solving. Passing 0 should linearize "
                 "all expressions requiring Xpress Global.",
                 &XpressmpBackend::GetGlobalFlag, &XpressmpBackend::SetGlobalFlag);

    AddSolverOption("alg:nlpsolver nlpsolver nlp:solver",
                    "Controls whether to call FICO Xpress Global or one "
                    "of the local solvers:\n"
                    "\n.. value-table::\n",
                    XPRS_NLPSOLVER, nlp_solver_values_, -1);

    AddSolverOption("alg:localsolver localsolver nlp:localsolver",
                    "Selects the library to use for local solves:\n"
                    "\n.. value-table::\n"
                    "\n"
                    "When solving problems to global optimality (see alg:nlpsolver), "
                    "alg:localsolver is used to decide which local solver to call "
                    "for reoptimizing NLP-infeasible solutions heuristically.",
                    XPRS_LOCALSOLVER, nlp_localsolver_values_, -1);

    //endalg
    // ****************************
    // PRESOLVER
    // ****************************
    AddSolverOption("pre:solve",
      "",
      XPRS_PRESOLVE, INT_MIN, INT_MAX);

    AddSolverOption("pre:solve_nlp presolve_nlp presolve_slp",
      "",
      XSLP_PRESOLVE, nlp_presolve_values_, 1);

    AddSolverOption("pre:ops",
      "",
      XPRS_PRESOLVEOPS, presolveops_values_, 511);

    AddSolverOption("pre:maxscalefactor",
      "",
      XPRS_MAXSCALEFACTOR, 0, 64);

    AddSolverOption("pre:configuration",
      "",
      XPRS_PRECONFIGURATION, INT_MIN, INT_MAX);

    AddSolverOption("pre:elimfillin",
      "",
      XPRS_ELIMFILLIN, 0, INT_MAX);

    AddSolverOption("pre:elimtol",
      "",
      XPRS_ELIMTOL, 0.0, DBL_MAX);

    AddSolverOption("pre:genconsdualreductions",
      "",
      XPRS_GENCONSDUALREDUCTIONS, INT_MIN, INT_MAX);

    AddSolverOption("pre:indlinbigbigm",
      "",
      XPRS_INDPRELINBIGM, 0.0, DBL_MAX);

    AddSolverOption("pre:maximpliedbound",
      "",
      XPRS_MAXIMPLIEDBOUND, 0.0, DBL_MAX);


    AddSolverOption("pre:objscalefactor",
      "",
      XPRS_OBJSCALEFACTOR, 0, INT_MAX);

    AddSolverOption("pre:basisred",
      "",
      XPRS_PREBASISRED, INT_MIN, INT_MAX);

    AddSolverOption("pre:bndredcone",
      "",
      XPRS_PREBNDREDCONE, INT_MIN, INT_MAX);

    AddSolverOption("pre:bndredquad",
      "",
      XPRS_PREBNDREDQUAD, INT_MIN, INT_MAX);

    AddSolverOption("pre:cliquestrategy",
      "",
      XPRS_PRECLIQUESTRATEGY, -1, INT_MAX);

    AddSolverOption("pre:coefelim",
      "",
      XPRS_PRECOEFELIM, INT_MIN, INT_MAX);
    
    AddSolverOption("pre:components",
      "",
      XPRS_PRECOMPONENTS, INT_MIN, INT_MAX);

    AddSolverOption("pre:componentseffort",
      "",
      XPRS_PRECOMPONENTSEFFORT, 0.0, DBL_MAX);

    AddSolverOption("pre:convertseparable",
      "",
      XPRS_PRECONVERTSEPARABLE, INT_MIN, INT_MAX);

    AddSolverOption("pre:domcol",
      "",
      XPRS_PREDOMCOL, INT_MIN, INT_MAX);

    AddSolverOption("pre:domrow",
      "",
      XPRS_PREDOMROW, INT_MIN, INT_MAX);

    AddSolverOption("pre:duprow",
      "",
      XPRS_PREDUPROW, INT_MIN, INT_MAX);

    AddSolverOption("pre:elimquad",
      "",
      XPRS_PREELIMQUAD, INT_MIN, INT_MAX);

    AddSolverOption("pre:folding",
      "",
      XPRS_PREFOLDING, INT_MIN, INT_MAX);

    AddSolverOption("pre:implications",
      "",
      XPRS_PREIMPLICATIONS, INT_MIN, INT_MAX);

    AddSolverOption("pre:lindep",
      "",
      XPRS_PRELINDEP, INT_MIN, INT_MAX);

    AddSolverOption("pre:objcutdetect",
      "",
      XPRS_PREOBJCUTDETECT, INT_MIN, INT_MAX);

    AddSolverOption("pre:permute",
      "",
      XPRS_PREPERMUTE, INT_MIN, INT_MAX);

    AddSolverOption("tech:permseed permseed pre:permuteseed permuteseed",
      "",
      XPRS_PREPERMUTESEED, -INT_MAX, INT_MAX);

    AddSolverOption("pre:probing",
      "",
      XPRS_PREPROBING, INT_MIN, INT_MAX);

    AddSolverOption("pre:protectdual",
      "",
      XPRS_PREPROTECTDUAL, INT_MIN, INT_MAX);

    AddSolverOption("pre:maxgrow",
      "",
      XPRS_PRESOLVEMAXGROW, 0.0, DBL_MAX);

    AddSolverOption("pre:passes",
      "",
      XPRS_PRESOLVEPASSES, 0, INT_MAX);

    AddSolverOption("pre:pwldualreductions",
      "",
      XPRS_PWLDUALREDUCTIONS, INT_MIN, INT_MAX);

    AddSolverOption("pre:pwlnonconvextransformation",
      "",
      XPRS_PWLNONCONVEXTRANSFORMATION, INT_MIN, INT_MAX);

    AddSolverOption("pre:rootpresolve",
      "",
      XPRS_ROOTPRESOLVE, INT_MIN, INT_MAX);

    AddSolverOption("pre:scaling",
      "",
      XPRS_SCALING, INT_MIN, INT_MAX);

    AddSolverOption("pre:trace",
      "",
      XPRS_TRACE, INT_MIN, INT_MAX);

    AddSolverOption("pre:sosreftol",
      "",
      XPRS_SOSREFTOL, 1e-6, Infinity());

    //endpre
  // ****************************
  // BARRIER ALGORITHM CONTROLS
  // ****************************
  AddSolverOption("bar:alg",
    "",
    XPRS_BARALG, INT_MIN, INT_MAX);

  AddSolverOption("bar:cachesize cachesize",
    "Deprecated: ignored by the current Xpress library.",
    XPRS_CACHESIZE, -1, INT_MAX);

    AddSolverOption("bar:l1cache l1cache",
    "Deprecated: ignored by the current Xpress library.",
                    XPRS_L1CACHE, -1, INT_MAX);

  AddSolverOption("bar:corespercpu",
    "",
    XPRS_CORESPERCPU, -1, INT_MAX);

  AddSolverOption("bar:cores",
    "",
    XPRS_BARCORES, -1, INT_MAX);

  AddSolverOption("bar:choleskyalg",
    "",
    XPRS_CHOLESKYALG, INT_MIN, INT_MAX);

  AddSolverOption("bar:choleskytol",
    "",
    XPRS_CHOLESKYTOL, 1e-15, DBL_MAX);

  AddSolverOption("bar:cpuplatform",
    "",
    XPRS_CPUPLATFORM, INT_MIN, INT_MAX);

  AddSolverOption("bar:crash",
    "",
    XPRS_BARCRASH, INT_MIN, INT_MAX);

  AddSolverOption("bar:crossover",
    "",
    XPRS_CROSSOVER, INT_MIN, INT_MAX);

  AddSolverOption("bar:crossoverops",
    "",
    XPRS_CROSSOVEROPS, INT_MIN, INT_MAX);

  AddSolverOption("bar:crossoverthreads",
    "",
    XPRS_CROSSOVERTHREADS, -1, INT_MAX);

  AddSolverOption("lim:crossoveriter bar:crossoveriterlim crossoveriterlim crossoveritlim",
    "",
    XPRS_CROSSOVERITERLIMIT, 1, INT_MAX);

  AddSolverOption("bar:crossovertol crossovertol",
    "",
    XPRS_CROSSOVERACCURACYTOL, 0.0, DBL_MAX);

  AddSolverOption("bar:densecollimit",
    "",
    XPRS_DENSECOLLIMIT, 0, INT_MAX);

  AddSolverOption("bar:dualstop",
    "",
    XPRS_BARDUALSTOP, 0.0, DBL_MAX);

  AddSolverOption("bar:hgextrapolate",
    "",
    XPRS_BARHGEXTRAPOLATE, 0.0, 1.0);

  AddSolverOption("bar:hgmaxrestarts",
    "",
    XPRS_BARHGMAXRESTARTS, 0, INT_MAX);

  AddSolverOption("bar:hgops",
    "",
    XPRS_BARHGOPS, INT_MIN, INT_MAX);

  AddSolverOption("bar:gapstop",
    "",
    XPRS_BARGAPSTOP, 0.0, DBL_MAX);

  AddSolverOption("bar:gap",
    "",
    XPRS_BARGAPTARGET, 0.0, DBL_MAX);

  AddSolverOption("bar:indeflimit",
    "",
    XPRS_BARINDEFLIMIT, 0, INT_MAX);

  AddSolverOption("lim:bariter bar:iterlim bariterlim",
    "",
    XPRS_BARITERLIMIT, 1, INT_MAX);

  AddSolverOption("bar:start",
    "",
    XPRS_BARSTART, INT_MIN, INT_MAX);

  AddSolverOption("bar:kernel",
    "",
    XPRS_BARKERNEL, -DBL_MAX, DBL_MAX);

  AddSolverOption("bar:objperturb",
    "",
    XPRS_BAROBJPERTURB, -DBL_MAX, DBL_MAX);

  AddSolverOption("bar:objscale",
    "",
    XPRS_BAROBJSCALE, -DBL_MAX, DBL_MAX);

  AddSolverOption("bar:order",
    "",
    XPRS_BARORDER, INT_MIN, INT_MAX);

  AddSolverOption("bar:orderthreads",
    "",
    XPRS_BARORDERTHREADS, 0, INT_MAX);

  AddSolverOption("bar:output",
    "",
    XPRS_BAROUTPUT, INT_MIN, INT_MAX);

  AddSolverOption("bar:presolve barpresolve",
    "",
    XPRS_BARPRESOLVEOPS, INT_MIN, INT_MAX);

  AddSolverOption("bar:primalstop",
    "",
    XPRS_BARPRIMALSTOP, 0.0, DBL_MAX);

  AddSolverOption("bar:refiter",
    "",
    XPRS_BARREFITER, INT_MIN, INT_MAX);

  AddSolverOption("bar:regularize barreg barrregularize",
    "",
    XPRS_BARREGULARIZE, INT_MIN, INT_MAX);

  AddSolverOption("bar:stepstop",
    "",
    XPRS_BARSTEPSTOP, 1e-10, DBL_MAX);

  AddSolverOption("bar:threads",
    "",
    XPRS_BARTHREADS, -1, INT_MAX);
  //endbarrier
  // ****************************
  // SIMPLEX RELATED
  // ****************************
  AddSolverOption("lp:bigmmethod",
    "",
    XPRS_BIGMMETHOD, INT_MIN, INT_MAX);

  AddSolverOption("lp:bigm bigmpenalty",
    "",
    XPRS_BIGM, 1024.0, DBL_MAX);

  AddSolverOption("lp:crash",
    "",
    XPRS_CRASH, INT_MIN, INT_MAX);

  AddSolverOption("lp:dualgradient",
    "",
    XPRS_DUALGRADIENT, INT_MIN, INT_MAX);

  AddSolverOption("lp:dualize",
    "",
    XPRS_DUALIZE, INT_MIN, INT_MAX);

  AddSolverOption("lp:dualstrategy",
    "",
    XPRS_DUALSTRATEGY, INT_MIN, INT_MAX);

  AddSolverOption("lp:dualthreads",
    "",
    XPRS_DUALTHREADS, -1, INT_MAX);

  AddSolverOption("lp:dualizeops",
    "",
    XPRS_DUALIZEOPS, INT_MIN, INT_MAX);

  AddSolverOption("lp:dualperturb",
    "",
    XPRS_DUALPERTURB, -1.0, DBL_MAX);

  AddSolverOption("lp:etatol",
    "",
    XPRS_ETATOL, 0.0, DBL_MAX);

  AddSolverOption("lp:dualforceparallel dualforceparallel",
    "",
    XPRS_FORCEPARALLELDUAL, INT_MIN, INT_MAX);



    AddSolverOption("lp:invertfreq",
      "",
      XPRS_INVERTFREQ, -1, INT_MAX);

    AddSolverOption("lp:invertmin",
      "",
      XPRS_INVERTMIN, 3, INT_MAX);

    AddSolverOption("lp:keepbasis",
      "",
      XPRS_KEEPBASIS, INT_MIN, INT_MAX);

    AddSolverOption("lp:keepnrows",
      "",
      XPRS_KEEPNROWS, INT_MIN, INT_MAX);

    AddSolverOption("lp:log",
      "",
      XPRS_LPLOG, -INT_MAX, INT_MAX);

    AddSolverOption("lp:netstalllimit",
      "",
      XPRS_NETSTALLLIMIT, INT_MIN, INT_MAX);

   
    AddSolverOption("lp:opttol lp:optimalitytol opttol",
      "",
      XPRS_OPTIMALITYTOL, 0.0, DBL_MAX);


    AddSolverOption("lp:optimalitytoltarget",
      "",
      XPRS_OPTIMALITYTOLTARGET, 0.0, DBL_MAX);

    AddSolverOption("lp:penalty",
      "",
      XPRS_PENALTY, 0.0, DBL_MAX);

    AddSolverOption("lp:pricingalg",
      "",
      XPRS_PRICINGALG, INT_MIN, INT_MAX);

    AddSolverOption("lp:primalunshift",
      "",
      XPRS_PRIMALUNSHIFT, INT_MIN, INT_MAX);

    AddSolverOption("lp:relpivottol",
      "",
      XPRS_RELPIVOTTOL, 0.0, Infinity());

    AddSolverOption("lp:sifting",
      "",
      XPRS_SIFTING, INT_MIN, INT_MAX);

    AddSolverOption("lp:siftpasses",
      "",
      XPRS_SIFTPASSES, 1, INT_MAX);

    AddSolverOption("lp:siftpresolveops",
      "",
      XPRS_SIFTPRESOLVEOPS, INT_MIN, INT_MAX);

    AddSolverOption("lp:siftswitch",
      "",
      XPRS_SIFTSWITCH, INT_MIN, INT_MAX);


    // endlp
  // MIP
  // ****************************
  AddSolverOption("mip:branchchoice",
    "",
    XPRS_BRANCHCHOICE, INT_MIN, INT_MAX);

  AddSolverOption("mip:branchdisj",
    "",
    XPRS_BRANCHDISJ, INT_MIN, INT_MAX);

  AddSolverOption("mip:branchstructural branchstruct",
    "",
    XPRS_BRANCHSTRUCTURAL, INT_MIN, INT_MAX);

  AddSolverOption("mip:breadthfirst",
    "",
    XPRS_BREADTHFIRST, 11, INT_MAX);
  
  AddSolverOption("mip:deterministic",
    "",
    XPRS_DETERMINISTIC, INT_MIN, INT_MAX);

  AddSolverOption("mip:feasibilityjump",
    "",
    XPRS_FEASIBILITYJUMP, INT_MIN, INT_MAX);

  AddSolverOption("mip:feasibilitypump",
    "",
    XPRS_FEASIBILITYPUMP, INT_MIN, INT_MAX);

    AddSolverOption("mip:heurbeforelp",
      "",
      XPRS_HEURBEFORELP, INT_MIN, INT_MAX);

    AddSolverOption("lim:heurdiveiterlimit mip:heurdiveiterlimit",
      "",
      XPRS_HEURDIVEITERLIMIT, -DBL_MAX, DBL_MAX);


    AddSolverOption("mip:heurshiftprop",
      "",
      XPRS_HEURSHIFTPROP, INT_MIN, INT_MAX);
      
    AddSolverOption("mip:heurdiverandomize hdive_rand",
      "",
      XPRS_HEURDIVERANDOMIZE, 0.0, 1.0);

    AddSolverOption("mip:heurdivesoftrounding hdive_rounding",
      "",
      XPRS_HEURDIVESOFTROUNDING, INT_MIN, INT_MAX);

    AddSolverOption("mip:heurdivespeedup hdive_speed",
      "",
      XPRS_HEURDIVESPEEDUP, INT_MIN, INT_MAX);

    AddSolverOption("mip:heurdivestrategy hdive_strategy",
      "",
      XPRS_HEURDIVESTRATEGY, INT_MIN, INT_MAX);

    AddSolverOption("mip:heuremphasis",
      "",
      XPRS_HEUREMPHASIS, INT_MIN, INT_MAX);

    AddSolverOption("alg:numericalemphasis alg:numericfocus numericfocus numfocus numericemphasis",
      "",
      XPRS_NUMERICALEMPHASIS, INT_MIN, INT_MAX);

    AddSolverOption("mip:heurforcespecialobj heurforcespecobj",
      "",
      XPRS_HEURFORCESPECIALOBJ, INT_MIN, INT_MAX);

    AddSolverOption("mip:heurfreq",
      "",
      XPRS_HEURFREQ, -1, INT_MAX);

    AddSolverOption("mip:heursearcheffort",
      "",
      XPRS_HEURSEARCHEFFORT, 1.0, DBL_MAX);

    AddSolverOption("mip:heursearchfreq",
      "",
      XPRS_HEURSEARCHFREQ, INT_MIN, INT_MAX);

    AddSolverOption("mip:heursearchrootcutfreq heurrootcutfreq",
      "",
      XPRS_HEURSEARCHROOTCUTFREQ, INT_MIN, INT_MAX);


    AddSolverOption("mip:heursearchrootselect",
      "",
      XPRS_HEURSEARCHROOTSELECT, INT_MIN, INT_MAX);

    AddSolverOption("mip:heursearchtreeselect",
      "",
      XPRS_HEURSEARCHTREESELECT, INT_MIN, INT_MAX);

    AddSolverOption("mip:heurthreads",
      "",
      XPRS_HEURTHREADS, INT_MIN, INT_MAX);

    AddSolverOption("mip:historycosts",
      "",
      XPRS_HISTORYCOSTS, INT_MIN, INT_MAX);

    AddSolverOption("mip:localchoice",
      "",
      XPRS_LOCALCHOICE, INT_MIN, INT_MAX);
    
AddSolverOption("mip:maxlocalbacktrack maxlocalbt",
      "",
      XPRS_MAXLOCALBACKTRACK, -1, INT_MAX);

AddSolverOption("mip:maxtasks",
  "",
  XPRS_MAXMIPTASKS, -1, INT_MAX);

AddSolverOption("mip:gap mipgap",
  "",
  XPRS_MIPRELSTOP, 1e-4, DBL_MAX);

AddSolverOption("mip:gapabs mipgapabs",
  "",
  XPRS_MIPABSSTOP, 0.0, DBL_MAX);

AddSolverOption("mip:components",
  "",
  XPRS_MIPCOMPONENTS, INT_MIN, INT_MAX);

AddSolverOption("mip:concurrentnodes",
  "",
  XPRS_MIPCONCURRENTNODES, INT_MIN, INT_MAX);

AddSolverOption("mip:concurrentsolves",
  "",
  XPRS_MIPCONCURRENTSOLVES, INT_MIN, INT_MAX);

AddSolverOption("mip:dualreductions",
  "The driver sets 2 when a solution pool is requested "
    "(sol:count or sol:stub).",
  XPRS_MIPDUALREDUCTIONS, INT_MIN, INT_MAX);


AddSolverOption("mip:kappafreq",
  "",
  XPRS_MIPKAPPAFREQ, INT_MIN, INT_MAX);

AddSolverOption("mip:log",
  "",
  XPRS_MIPLOG, -INT_MAX, INT_MAX);

AddSolverOption("mip:presolve",
  "",
  XPRS_MIPPRESOLVE, INT_MIN, INT_MAX);

AddSolverOption("mip:rampup",
  "",
  XPRS_MIPRAMPUP, INT_MIN, INT_MAX);

AddSolverOption("mip:miprefineiterlimit miprefiterlim",
  "",
  XPRS_MIPREFINEITERLIMIT, -1, INT_MAX);

AddSolverOption("mip:restart",
  "",
  XPRS_MIPRESTART, INT_MIN, INT_MAX);

AddSolverOption("mip:restartgapthreshold",
  "",
  XPRS_MIPRESTARTGAPTHRESHOLD, 0.0, DBL_MAX);

AddSolverOption("mip:restartfactor",
  "",
  XPRS_MIPRESTARTFACTOR, 0.0, DBL_MAX);

AddSolverOption("mip:threads",
  "",
  XPRS_MIPTHREADS, -1, INT_MAX);

AddSolverOption("mip:intfeastol inttol intfeastol",
  "",
  XPRS_MIPTOL, 0.0, DBL_MAX);

AddSolverOption("mip:toltarget",
  "",
  XPRS_MIPTOLTARGET, 0.0, DBL_MAX);

AddSolverOption("mip:nodeprobingeffort",
  "",
  XPRS_NODEPROBINGEFFORT, 0.0, DBL_MAX);

AddSolverOption("mip:nodeselection",
  "",
  XPRS_NODESELECTION, INT_MIN, INT_MAX);

AddSolverOption("mip:pseudocost",
  "",
  XPRS_PSEUDOCOST, 0.0, DBL_MAX);



AddSolverOption("mip:qcrootalg",
  "",
  XPRS_QCROOTALG, INT_MIN, INT_MAX);

AddSolverOption("mip:relaxtreememorylimit relaxtreemem",
  "",
  XPRS_RELAXTREEMEMORYLIMIT, 0.0, DBL_MAX);

AddSolverOption("mip:sbbest",
  "",
  XPRS_SBBEST, INT_MIN, INT_MAX);

AddSolverOption("mip:sbeffort",
  "",
  XPRS_SBEFFORT, 0.0, DBL_MAX);

AddSolverOption("mip:sbestimate",
  "",
  XPRS_SBESTIMATE, INT_MIN, INT_MAX);

AddSolverOption("mip:sbiterlimit",
  "",
  XPRS_SBITERLIMIT, -1, INT_MAX);

AddSolverOption("mip:sbselect",
  "",
  XPRS_SBSELECT, INT_MIN, INT_MAX);

AddSolverOption("mip:symmetry",
  "",
  XPRS_SYMMETRY, INT_MIN, INT_MAX);

AddSolverOption("mip:symselect",
  "",
  XPRS_SYMSELECT, INT_MIN, INT_MAX);

AddSolverOption("mip:varselection",
  "",
  XPRS_VARSELECTION, INT_MIN, INT_MAX);


// endmip
    // ****************************
  // Cuts
  // ****************************
  AddSolverOption("cut:cover",
    "",
    XPRS_COVERCUTS, -1, INT_MAX);

  AddSolverOption("cut:gomory",
    "",
    XPRS_GOMCUTS, -1, INT_MAX);

  AddSolverOption("cut:qccuts",
    "",
    XPRS_QCCUTS, -1, INT_MAX);

    AddSolverOption("cut:lnpbest",
      "",
      XPRS_LNPBEST, 1, INT_MAX);

    AddSolverOption("cut:lnpiterlimit",
      "",
      XPRS_LNPITERLIMIT, -1, INT_MAX);
    
    AddSolverOption("cut:rltcuts",
      "",
      XPRS_RLTCUTS, INT_MIN, INT_MAX);

  AddSolverOption("cut:treecover",
    "",
    XPRS_TREECOVERCUTS, -1, INT_MAX);

  AddSolverOption("cut:treegomory",
    "",
    XPRS_TREEGOMCUTS, -1, INT_MAX);

  AddSolverOption("cut:treeqccuts",
    "",
    XPRS_TREEQCCUTS, -1, INT_MAX);

  AddSolverOption("cut:depth",
    "",
    XPRS_CUTDEPTH, -1, INT_MAX);


  AddSolverOption("cut:factor",
    "",
    XPRS_CUTFACTOR, -1.0, DBL_MAX);

  AddSolverOption("cut:freq",
    "",
    XPRS_CUTFREQ, -1, INT_MAX);

  AddSolverOption("cut:select",
    "",
    XPRS_CUTSELECT, INT_MIN, INT_MAX);

  AddSolverOption("cut:treeselect",
    "",
    XPRS_TREECUTSELECT, INT_MIN, INT_MAX);

  AddSolverOption("cut:strategy",
    "",
    XPRS_CUTSTRATEGY, INT_MIN, INT_MAX);

  // ****************************
  // Q(C)P
  // ****************************
  AddStoredOption("qp:nonconvex nonconvex",
    "Deprecated, use qp:convexitychk. Values 0 and 1 are passed on to "
    "qp:convexitychk (convexity check off/on); the value 2, with which "
    "other AMPL drivers allow nonconvex quadratics, switches the "
    "convexity check off as well.",
    storedOptions_.nonconvex_);

  AddSolverOption("qp:eigenvaluetol",
    "",
    XPRS_EIGENVALUETOL, 0.0, DBL_MAX);

  AddSolverOption("qp:simplexops",
    "",
    XPRS_QSIMPLEXOPS, INT_MIN, INT_MAX);

  AddSolverOption("qp:miqcpalg",
    "",
    XPRS_MIQCPALG, INT_MIN, INT_MAX);

  AddSolverOption("qp:unshift quadunshift",
    "",
    XPRS_QUADRATICUNSHIFT, INT_MIN, INT_MAX);

  AddSolverOption("qp:repairindefiniteq repairindefq",
    "",
    XPRS_REPAIRINDEFINITEQ, INT_MIN, INT_MAX);

  //endqp



  AddStoredOption("alg:readbasis readbasis startbasis basisin",
    "File with the starting basis (Xpress .bas format, see "
    "tech:writesolution). It is read after any incoming .sstatus basis "
    "(see alg:basis) and overrides it. Xpress matches the basis by row "
    "and column names, so the model must be passed with names "
    "(AMPL option auxfiles rc).",
    storedOptions_.readBasis_);

  AddListOption("tech:optionnative optionnative optnative tech:param param",
    "Xpress control by name: a quoted string \"NAME=value\" or "
    "\"NAME value\", for documented and undocumented controls; "
    "can be repeated. \"NAME=?\" prints the current value.",
    storedOptions_.nativeOptions_);

  AddStoredOption("tech:tunebase tunerdir tunebase",
    "Base name for results of running XPRESS's search for best "
    "parameter settings. The search is run only when tunebase "
    "is specified.  This control only defines the root path for the tuner "
    "output. For each problem, the tuner result will be output to a subfolder "
    "underneath this path. For example, by default, the tuner result for a "
    "problem called prob will be located at tuneroutput/prob/",
    storedOptions_.tunebase_);

  AddStoredOption("tech:tunermethodread tunermethodread",
      "Read existing tuner method from the specified .xtm file, "
      "see \"tunermethodwrite\" to obtaing a template file",
      storedOptions_.tunermethodfileread_
  );
  AddStoredOption("tech:tunermethodwrite tunermethodwrite",
      "Write existing tuner method from the specified .xtm file",
      storedOptions_.tunermethodfilewrite_
  );
  AddStoredOption("tech:tunename tunesessionname",
    "Set problem name within the tuner \"tunebase\" is specified.",
    storedOptions_.tunename_);


  AddSolverOption("tech:tuneoutput tuneoutput",
    "",
    XPRS_TUNEROUTPUT, INT_MIN, INT_MAX);

  AddSolverOption("tech:tunerhistory",
    "",
    XPRS_TUNERHISTORY, INT_MIN, INT_MAX);

  AddSolverOption("tech:tunetimelim tunetimelim lim:tunetime",
    "",
    XPRS_TUNERMAXTIME, 0.0, Infinity());

  AddSolverOption("tech:tunerthreads",
    "",
    XPRS_TUNERTHREADS, -1, INT_MAX);

  AddSolverOption("tech:tunermethod",
    "",
    XPRS_TUNERMETHOD, INT_MIN, INT_MAX);

  AddSolverOption("tech:tunertarget",
    "",
    XPRS_TUNERTARGET, INT_MIN, INT_MAX);

  AddSolverOption("tech:tunerverbose",
    "",
    XPRS_TUNERVERBOSE, INT_MIN, INT_MAX);

  /////////////////////// Custom solve results /////////////////////////
//  AddSolveResults({
//                    { sol::NUMERIC, "failure: numeric issue, no feasible solution" }
//                  });

  AddSolverOption("obj:multi:ops",
    "",
    XPRS_MULTIOBJOPS, INT_MIN, INT_MAX);

  AddSolverOption("obj:multi:log",
    "",
    XPRS_MULTIOBJLOG, INT_MIN, INT_MAX);

  // Adding compiled options, see BUILD_XPRESS_OPTION_PARSER
  this->AddOptimizerOptions();
  this->AddNonlinearOptions();
}


double XpressmpBackend::MIPGap() {
  // implementation taken from ASL driver for consistency
  return MIPGapAbs() / (1e-10 + std::abs(ObjectiveValue()));
}
double XpressmpBackend::BestDualBound() {
  return getDblAttr(XPRS_BESTBOUND);
}

double XpressmpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> XpressmpBackend::VarStatii() {

  std::vector<int> vars(NumVarsOrig());
  int status = XPRSgetbasis(model_fixed_, NULL, vars.data());
  if (status)
    vars.clear();
  else 
    for (auto& s : vars) {
      switch (s) {
      case XPRS_BASISSTATUS_BASIC:
        s = (int)BasicStatus::bas;
        break;
      case XPRS_BASISSTATUS_NONBASIC_LOWER:
        s = (int)BasicStatus::low;
        break;
      case XPRS_BASISSTATUS_NONBASIC_UPPER:
        s = (int)BasicStatus::upp;
        break;
      case XPRS_BASISSTATUS_SUPERBASIC:
        s = (int)BasicStatus::sup;
        break;
      default:
        MP_RAISE(fmt::format("Unknown Xpress basis value: {}", s));
      }
    }
  return vars;
}

ArrayRef<int> XpressmpBackend::ConStatii() {
  std::vector<int> cons(NumLinConsOrig());
  int status = XPRSgetbasis(model_fixed_, cons.data(), NULL);
  if (status)
    cons.clear();
  else
    for (auto& s : cons) {
      switch (s) {
      case XPRS_BASISSTATUS_BASIC:
        s = (int)BasicStatus::bas;
        break;
      case XPRS_BASISSTATUS_NONBASIC_LOWER:
        s = (int)BasicStatus::low;
        break;
      case XPRS_BASISSTATUS_NONBASIC_UPPER:
        s = (int)BasicStatus::upp;
        break;
      case XPRS_BASISSTATUS_SUPERBASIC:
        s = (int)BasicStatus::sup;
        break;
      default:
        MP_RAISE(fmt::format("Unknown Xpress basis value: {}", s));
      }
    }
  return cons;
}

std::vector<int> XpressmpBackend::VarStatii(ArrayRef<int> vst) {
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  std::vector<double> lb, ub;
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = XPRS_BASISSTATUS_BASIC;
      break;
    case BasicStatus::low:
      s = XPRS_BASISSTATUS_NONBASIC_LOWER;
      break;
    case BasicStatus::upp:
      s = XPRS_BASISSTATUS_NONBASIC_UPPER;
      break;
    case BasicStatus::equ:
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = XPRS_BASISSTATUS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      if (lb.size() == 0)
      {
        lb.resize(vst.size());
        ub.resize(vst.size());
        XPRESSMP_CCALL(XPRSgetlb(lp(), lb.data(), 0, vst.size()-1));
        XPRESSMP_CCALL(XPRSgetub(lp(), ub.data(), 0, vst.size()-1));
      }
      if (lb[j] >= -1e-6)
        s = XPRS_BASISSTATUS_NONBASIC_LOWER;
      else if (ub[j] <= 1e-6)
        s = XPRS_BASISSTATUS_NONBASIC_UPPER;
      else
        s = XPRS_BASISSTATUS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }
  return stt;
}

std::vector<int> XpressmpBackend::ConStatii(ArrayRef<int> cst) {
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = XPRS_BASISSTATUS_BASIC;
      break;
    case BasicStatus::low:
      s = XPRS_BASISSTATUS_NONBASIC_LOWER;
      break;
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
      s = XPRS_BASISSTATUS_NONBASIC_UPPER;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = XPRS_BASISSTATUS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  return stt;
}

SolutionBasis XpressmpBackend::GetBasis() {
  if (fixed_model_failed_) return {};
  int old_outlev = outlev_;
  outlev_ = 0;
    std::vector<int> varstt = VarStatii();
    std::vector<int> constt = ConStatii();
  outlev_ = old_outlev;
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { varstt,
        {{{ CG_Algebraic, constt }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
  }
  return { varstt,constt};
}

void XpressmpBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto &varstt = mv.GetVarValues()();
  auto &constt = mv.GetConValues()(CG_Algebraic);
//#define XPRESS__ROW_STATS_GENCONS    // Is this valid?
#ifdef XPRESS__ROW_STATS_GENCONS
  // Append general constraints. TODO
  auto& cconstt = mv.GetConValues()(CG_General);
  constt.insert(constt.end(), cconstt.begin(), cconstt.end());
#endif  // XPRESS__ROW_STATS_GENCONS
  if (varstt.size() && constt.size()) {
    auto convertedVarBasis = VarStatii(varstt);
    auto convertedConBasis = ConStatii(constt);
    // Xpress needs to have all constraints, append them as basic
    convertedConBasis.resize(getIntAttr(XPRS_ROWS), XPRS_BASISSTATUS_BASIC);     
    
    XPRESSMP_CCALL(XPRSloadbasis(lp(), convertedConBasis.data(), convertedVarBasis.data()));
  }
}


void XpressmpBackend::ComputeIIS() {
  int status;
  XPRESSMP_CCALL(XPRSiisfirst(lp(), 1, &status));
//  switch (status) {
//  case XPRS_IIS_COMPLETED:
  SetStatus( GetSolveResult() );
}

IIS XpressmpBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}
IISStatus XpressmpBackend::IIS_VarToAMPL(char c) {
  switch (c)
  {
  case 'L':
    return IISStatus::low;
  case 'U':
    return IISStatus::upp;
  case 'F':
    return IISStatus::fix;
  case 'B':
  case 'I':
  case 'P':
    return IISStatus::intvar;
  case 'S':
  case 'R':
    return IISStatus::semi;
  default:
    AddWarning("IISVarStatus",
               fmt::format("Unknown status '{}' of an IIS variable", c));
  }
  return IISStatus::mem;
}

ArrayRef<int> XpressmpBackend::VarsIIS() {
  int nconsiis, nvarsiis;
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, 0, 0, 0, 0, 0, 0, 0, 0));
  std::vector<int> vars(nvarsiis);
  std::vector<char> bounds(nvarsiis);
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, 0,
    vars.data(), 0, bounds.data(), 0, 0, 0, 0));
  std::vector<int> iis(NumVars(), (int)IISStatus::non);
  for (int i = 0; i < nvarsiis; i++)
    iis[vars[i]] = (int)IIS_VarToAMPL(bounds[i]);
  return iis;
}

pre::ValueMapInt XpressmpBackend::ConsIIS() {
  int nconsiis, nvarsiis;
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, 0, 0, 0, 0, 0, 0, 0, 0));
  std::vector<int> cons(nconsiis);
  std::vector<char> contype(nconsiis);
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, cons.data(),
    0, contype.data(), 0, 0, 0, 0, 0));
  std::vector<int> iis_algcon(NumLinCons(), (int)IISStatus::non);
  std::vector<int> iis_sos(NumSOSCons(), (int)IISStatus::non);
  std::vector<int> iis_pl(NumPWLs(), (int)IISStatus::non);
  std::vector<int> iis_gencon(NumGenCons(), (int)IISStatus::non);
  std::vector<int> iis_indcon(NumIndicatorCons(), (int)IISStatus::non);
  for (int i = 0; i < nconsiis; i++)
  {
    switch (contype[i]) {
    case 'G':
      iis_algcon[cons[i]] = (int)IISStatus::low;
      break;
    case 'L':
      iis_algcon[cons[i]] = (int)IISStatus::upp;
      break;
    case 'E':
      iis_algcon[cons[i]] = (int)IISStatus::fix;
      break;
    case '1':
    case '2':
      iis_sos[cons[i]] = (int)IISStatus::mem;
      break;
    case 'W':
      iis_pl[cons[i]] = (int)IISStatus::mem;
      break;
    case 'X':
      iis_gencon[cons[i]] = (int)IISStatus::mem;
      break;
    case 'I':
      iis_indcon[cons[i]] = (int)IISStatus::mem;
      break;
    default:
      AddWarning("IISVarStatus",
                 fmt::format("Unknown status '{}' of an IIS constraint",
                             contype[i]));
    }
  }
  return {{{ CG_Algebraic, iis_algcon },
           { CG_SOS, iis_sos },
           { CG_Piecewiselinear, iis_pl },
           { CG_General, iis_gencon },
           { CG_Indicator, iis_indcon }}};
}

void XpressmpBackend::AddPrimalDualStart(Solution sol0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution(
    { sol0_unpres.primal, sol0_unpres.dual });
  auto& x0 = mv.GetVarValues()();
  auto& pi0 = mv.GetConValues()(CG_Algebraic);

  int status;
  XPRESSMP_CCALL(XPRSloadlpsol(lp(), x0.data(), NULL,
    pi0.data(), NULL, &status));
  if (status)
    Print("warmstart: solution is not loaded "
          "because the problem is in presolved status.\n");
}

void XpressmpBackend::AddMIPStart(
		ArrayRef<double> x0_unpres, ArrayRef<int> s0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution({ x0_unpres });
  auto ms = GetValuePresolver().PresolveGenericInt({ s0_unpres });
  auto& x0 = mv.GetVarValues()();
  auto& s0 = ms.GetVarValues()();
  std::vector<int> idx;                 // Create sparse vector
  idx.reserve(x0.size());
  std::vector<double> val;
  val.reserve(x0.size());
  for (int i = 0; i < (int)x0.size(); ++i) {
    if (s0[i]) {
      idx.push_back(i);
      val.push_back(x0[i]);
    }
  }
  XPRESSMP_CCALL(
      XPRSaddmipsol(lp(), idx.size(), val.data(), idx.data(), "AMPL_initial_guess"));
  XPRESSMP_CCALL(  // in MIPStart() because NLP start is primal-only
      XPRSnlpsetinitval(lp(), idx.size(), idx.data(), val.data()));
}

void XpressmpBackend::xpdisplay(
    XPRSprob prob, void* data, const char* ch, int n, int msglvl)
{
  /*
   msglvl gives the message level as follows:
   * 1 dialogue
   * 2 information
   * 3 warnings
   * 4 errors
   * a negative value indicates the XPRESS is about to finish and
   * buffers should be flushed.
   */
  if (outlev_ == 0)
    return;
  if (msglvl < 0)
    fflush(NULL);
  else if (msglvl >= outlev_ && (msglvl != 4 || strncmp(ch, "?899 ", 5)))
  {
    fmt::print("{}\n", ch);
    fflush(stdout);
  }
}

int XpressmpBackend::xp_mse_display(
    XPRSobject o, void* context, void* thread,
    const char* ch, int msglvl, int msgnumber) {
  if (outlev_ == 0) return 0;
  if (msglvl < 0)
    fflush(NULL);
  else if (msglvl >= outlev_ && (msglvl != 4)) {
    fmt::print("{}\n", ch);
    fflush(stdout);
  }
  return 0;
}

} // namespace mp

// AMPLs
AMPLS_MP_Solver* Open_xpress(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::XpressmpBackend()},
   cb);
}

void AMPLSClose_xpress(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_xpress(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::XpressmpBackend*>(AMPLSGetBackend(slv))->lp();
}

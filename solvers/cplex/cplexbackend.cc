#include <vector>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cassert>
#include <fstream>
#include <cctype>
#include <iterator>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cplexbackend.h"

extern "C" {
#include "cplex-ampls-c-api.h"    // CPLEX AMPLS C API
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
  /// @return CplexModelMgr
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
      auto msg = CPXgeterrorstring(env(), status, errmsg);
      throw std::runtime_error(
        fmt::format("Could not open CPLEX environment.\n{}", msg ? msg : ""));
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

    if (lp() != nullptr) {
      CPLEX_CALL(CPXfreeprob(env(), &lp_ref()));
    }
    /* Free up the CPLEX environment, if necessary */
    if (env() != nullptr) {
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



#define getAndReturnDblParam(function)\
  double value;\
  CPLEX_CALL(function(env(), lp(), &value));\
  return value;

#define getDblParam(function, var)\
  double var;\
  CPLEX_CALL(function(env(), lp(), &var));


  bool CplexBackend::IsMIP() const {
    auto type = CPXgetprobtype(env(), lp());
    return !((type == CPXPROB_LP) || (type == CPXPROB_QP)
      || (type == CPXPROB_QCP));
  }
  bool CplexBackend::IsQCP() const {
    int probtype = CPXgetprobtype(env(), lp());
    return probtype >= (int)CPXPROB_QP;
  }

  bool CplexBackend::HasSolution() {
    if (hasSolution_ == -1){
      auto status = GetSolveResult().first;
      hasSolution_
          = (status == sol::LIMIT_FEAS)
            || (status == sol::SOLVED)
            || (status == sol::UNCERTAIN);
    }
    return hasSolution_;
  }

  double  CplexBackend::MIPGap() {
    if (have_mip_result_before_fixed_) return mip_result_before_fixed_.gap;
    if (HasSolution()) {
      getAndReturnDblParam(CPXgetmiprelgap);
    }
    else return AMPLInf();
  }

  double  CplexBackend::MIPGapAbs() {
    if (have_mip_result_before_fixed_) return mip_result_before_fixed_.absolute_gap;
    auto type = CPXgetprobtype(env(), lp());
    if ((type == CPXPROB_LP) ||
      (type == CPXPROB_QP) || (type == CPXPROB_QCP))
      return 0;

    if (!HasSolution())
      return AMPLInf(); // no solution found

    double obj;
    int status = CPXgetobjval(env(), lp(), &obj);
    if (status)
      return AMPLInf(); // no solution found
    return std::abs(obj - BestDualBound());
  }
  double  CplexBackend::BestDualBound() {
    if (have_mip_result_before_fixed_) return mip_result_before_fixed_.bound;
    auto type = CPXgetprobtype(env(), lp());
    if ((type == CPXPROB_LP) ||
      (type == CPXPROB_QP) || (type == CPXPROB_QCP))
      return 0;
    getDblParam(CPXgetbestobjval, bobj);
    if (bobj == Infinity())
      return AMPLInf();
    if (bobj == -Infinity())
      return -AMPLInf();
    return bobj;
  }

  ArrayRef<int> CplexBackend::VarStatii() {
    std::vector<int> vars(NumVars());
    int status = CPXgetbase(env(), lp(), vars.data(), nullptr);
    if (status)
      vars.clear();
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
    std::vector<int> cons(NumLinCons());
    int status = CPXgetbase(env(), lp(), nullptr, cons.data());
    if (status)
      cons.clear();
    else {
      std::vector<char> sense(NumLinCons());
      CPLEX_CALL(CPXgetsense(env(), lp(), sense.data(), 0, NumLinCons()-1));

      for (auto i = cons.size(); i--; ) {
        switch (cons[i]) {
        case CPX_BASIC:
          cons[i] = (int)BasicStatus::bas;
          break;
        case CPX_AT_LOWER:
          if (sense[i] == 'L')
            cons[i] = (int)BasicStatus::upp;
          else if (sense[i] == 'E')
            cons[i] = (int)BasicStatus::equ;
          else
            cons[i] = (int)BasicStatus::low;
          break;
        case CPX_AT_UPPER: // just for range constraints
          cons[i] = (int)BasicStatus::upp;
          break;
        default:
          MP_RAISE(fmt::format("Unknown CPLEX rstat value: {}", cons[i]));
        }
      }
    }
    return cons;
  }

  void CplexBackend::VarConStatii(ArrayRef<int> vstt, ArrayRef<int> cstt) {
    std::vector<int> vst
        = std::vector<int>(vstt.data(), vstt.data() + vstt.size());
    std::vector<int> cst
        = std::vector<int>(cstt.data(), cstt.data() + cstt.size());
    for (auto j = vst.size(); j--; ) {
      auto& s = vst[j];
      switch ((BasicStatus)s) {
      case BasicStatus::bas:
        s = (int)CPX_BASIC;
        break;
      case BasicStatus::low:
      case BasicStatus::equ:
        s = CPX_AT_LOWER;
        break;
      case BasicStatus::none:
        double lb, ub;
        if (!CPXgetlb(env(), lp(), &lb, j, j) &&
            !CPXgetub(env(), lp(), &ub, j, j))
        {
          if (lb > MinusInfinity())
            s = CPX_AT_LOWER;
          else if (ub < Infinity())
            s = CPX_AT_UPPER;
          else
            s = CPX_FREE_SUPER;  
        }
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
    if (fixed_model_failed_) return {};
    bool outputStatus = silenceOutput_;
    setSilenceOutput(!debug_mode());
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
    setSilenceOutput(outputStatus);
    return { std::move(varstt), std::move(constt) };
  }

  void CplexBackend::SetBasis(SolutionBasis basis) {
    auto mv = GetValuePresolver().PresolveBasis(
      { basis.varstt, basis.constt });
    auto& varstt = mv.GetVarValues()();
    auto& constt = mv.GetConValues()(CG_Linear);
    assert(varstt.size());
    assert(constt.size());
    VarConStatii(varstt, constt);
  }

  void CplexBackend::AddPrimalDualStart(Solution sol0_unpres) {
    auto mv = GetValuePresolver().PresolveSolution(
      { sol0_unpres.primal, sol0_unpres.dual });
    auto& x0 = mv.GetVarValues()();
    auto& pi0 = mv.GetConValues()(CG_Linear);
    CPXcopystart(env(), lp(), nullptr, nullptr, x0.data(), nullptr,
      nullptr, pi0.data());
  }

  void CplexBackend::AddMIPStart(
    ArrayRef<double> x0_unpres, ArrayRef<int> sparsity_unpres) {
    if (!IsMIP()) return;
    auto mv = GetValuePresolver().PresolveSolution({ x0_unpres });
    auto ms = GetValuePresolver().PresolveGenericInt({ sparsity_unpres });
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
    int beg[2] = { 0, static_cast<int>(val.size()) };
    CPLEX_CALL(CPXaddmipstarts(env(), lp(), 1, val.size(), beg, idx.data(), val.data(), 0, nullptr));
  }

ArrayRef<double> CplexBackend::PrimalSolution() {
  if (fixed_model_failed_) return mip_result_before_fixed_.primal;
  // if (!HasSolution())
  //   return std::vector<double>();
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
  if (fixed_model_failed_) return {};
  // if (HasSolution() && ((!IsMIP()) || need_fixed_MIP()) && (NumLinCons()>0))
  {
    int num_cons = NumLinCons();
    std::vector<double> pi(num_cons);
    /* Avoid most error messages on screen */
    CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ScreenOutput, CPX_OFF));
    int error = CPXgetpi(env(), lp(), pi.data(), 0, num_cons - 1);
    if (error)
      pi.clear();
    return pi;
  }
  // else
  //   return std::vector<double>();
}

double CplexBackend::ObjectiveValue() const {
  if (fixed_model_failed_) return mip_result_before_fixed_.objective;
  double objval = -Infinity();
  if(hasSolution_)
    CPLEX_CALL(CPXgetobjval(env(), lp(), &objval));
  return objval;
}

int CplexBackend::NodeCount() const {
  if (have_mip_result_before_fixed_) return mip_result_before_fixed_.nodes;
  return  IsMIP() ? CPXgetnodecnt(env(), lp()) : 0;
}

int CplexBackend::SimplexIterations() const {
  if (have_mip_result_before_fixed_) return mip_result_before_fixed_.simplex_iterations;
  return IsMIP() ? CPXgetmipitcnt(env(), lp()) : CPXgetitcnt(env(), lp());
}

int CplexBackend::BarrierIterations() const {
  if (have_mip_result_before_fixed_) return mip_result_before_fixed_.barrier_iterations;
  if (!IsMIP() && (storedOptions_.cpxMethod_ == CPX_ALG_BARRIER))
    return CPXgetbaritcnt(env(), lp());
  else
    return 0;
}
std::map<std::string, std::variant<int, double, std::string>>
CplexBackend::SolutionStats() {
    std::map<std::string, std::variant<int, double, std::string>> stats;
    stats["simplex_iterations"] = SimplexIterations();
    stats["barrier_iterations"] = BarrierIterations();
    stats["node_count"] = NodeCount();
    return stats;
}
void CplexBackend::DoWriteProblem(const std::string &file) {
  CPLEX_CALL( CPXwriteprob (env(), lp(), file.c_str(), NULL) );
}
void CplexBackend::DoWriteSolution(const std::string& file) {
  auto CheckEnds = [&](std::string e) {
    return ends_with(file, e)
           || ends_with(file, e + ".gz")
           || ends_with(file, e + ".zip");
  };
  if (CheckEnds(".clp")) {
    CPLEX_CALL(CPXclpwrite(env(), lp(), file.c_str()));
    return;
  } else if (CheckEnds(".bas")) {
    CPLEX_CALL(CPXmbasewrite(env(), lp(), file.c_str()));
    return;
  } else if (CheckEnds(".sol")) {
    CPLEX_CALL(CPXsolwrite(env(), lp(), file.c_str()));
    return;
  }
  AddWarning("UnknownResultExt",
             fmt::format("Unknown or non-supported result extension: '{}'.\n"
                         "Trying to write the file as solution", file));
  CPLEX_CALL(CPXsolwrite(env(), lp(), file.c_str()));
}


void CplexBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCplex, nullptr);
  CPLEX_CALL( CPXsetterminate (env(), &terminate_flag) );
}
void CplexBackend::RedirectOutput() {
  // Redirect output
  CPXCHANNELptr cpxresults, cpxwarnings, cpxerrors, cpxlog;

  CPLEX_CALL(CPXgetchannels(env(), &cpxresults, &cpxwarnings,
    &cpxerrors, &cpxlog));
  
  auto msgCallback = [](void* context, const char* msg) {
    CplexBackend* instance = static_cast<CplexBackend*>(context);
    instance->mymsgfunc(msg, "");
    };
  auto warCallback = [](void* context, const char* msg) {
    CplexBackend* instance = static_cast<CplexBackend*>(context);
    instance->mymsgfunc(msg, "");
    };
  auto errCallback = [](void* context, const char* msg) {
    CplexBackend* instance = static_cast<CplexBackend*>(context);
    instance->mymsgfunc(msg, "");
    };
  auto logCallback = [](void* context, const char* msg) {
    CplexBackend* instance = static_cast<CplexBackend*>(context);
    instance->mymsgfunc(msg, "");
    };
  CPXdisconnectchannel(env(), cpxresults);
  CPXdisconnectchannel(env(), cpxwarnings);
  CPXdisconnectchannel(env(), cpxerrors);
  CPXdisconnectchannel(env(), cpxlog);

  CPXaddfuncdest(env(), cpxresults, this, msgCallback);
  CPXaddfuncdest(env(), cpxerrors, this, errCallback);
  CPXaddfuncdest(env(), cpxwarnings, this, warCallback);
  CPXaddfuncdest(env(), cpxlog, this, logCallback);

}

void CplexBackend::Solve() {
  have_mip_result_before_fixed_ = false;
  fixed_model_failed_ = false;
    if (!storedOptions_.tunefileprm.empty())
        DoCPLEXTune();

  if (storedOptions_.noSolve_)
  {
    fmt::print("Not solved because of \"nosolve\"");
    exit(1);
  }

  RedirectOutput();
  setSilenceOutput(storedOptions_.outlev_==0);
  setSolutionMethod();

  if (storedOptions_.dropTol_ > 0) {
    CPLEX_CALL(CPXcleanup(env(), lp(), storedOptions_.dropTol_));
  }
  if (NumObjs() > 1
      || multiobj()==1)                  // User set multiobj=1
    CPLEX_CALL(CPXmultiobjopt(env(), lp(), NULL));
  else {
    auto type = CPXgetprobtype(env(), lp());
    if ((type == CPXPROB_MIQCP) || (type == CPXPROB_MIQP) || (type == CPXPROB_MILP))
    {
      if (storedOptions_.cpxMethod_ == CPX_ALG_BENDERS)
        CPLEX_CALL(CPXbendersopt(env(), lp()));
      else
        CPLEX_CALL(CPXmipopt(env(), lp()));
    }
    else if ((type == CPXPROB_QP) || (type == CPXPROB_QCP))
      CPLEX_CALL(CPXqpopt(env(), lp()));
    else
      CPLEX_CALL(CPXlpopt(env(), lp()));
  }
  if (feasrelax())
  {
    auto solstatus = CPXgetstat(env(), lp());
    if((solstatus==CPX_STAT_INFEASIBLE)||(solstatus==CPXMIP_INFEASIBLE) ||
        (solstatus==CPX_STAT_INForUNBD) || (solstatus== CPXMIP_INForUNBD) ||
        (solstatus== CPXMIP_FAIL_INFEAS))
      DoCplexFeasRelax();

  }
  setSilenceOutput(!debug_mode());
  WindupCPLEXSolve();
}

ArrayRef<double> CplexBackend::GetObjectiveValues() 
{ 
  if(NumObjs()==1)
    return std::vector<double>{ObjectiveValue()}; 
  else {
    std::vector<double> vals(NumObjs());
    for (int i = 0; i < NumObjs(); i++) {
      CPXmultiobjgetobjval(env(), lp(), i, &vals[i]);
    }
    return vals;
  }
} 

void CplexBackend::ObjPriorities(ArrayRef<int> priority) {
  for (int i = 0; i < (int)priority.size(); ++i) {
    CPXmultiobjchgattribs(env(), lp(), i,
      CPX_NO_OFFSET_CHANGE, CPX_NO_WEIGHT_CHANGE,
      priority[i], CPX_NO_ABSTOL_CHANGE, CPX_NO_RELTOL_CHANGE,
      NULL);
  }
}

void CplexBackend::ObjWeights(ArrayRef<double> val) {
  for (int i = 0; i < (int)val.size(); ++i) {
    CPXmultiobjchgattribs(env(), lp(), i,
      CPX_NO_OFFSET_CHANGE, val[i],
      CPX_NO_PRIORITY_CHANGE, CPX_NO_ABSTOL_CHANGE, CPX_NO_RELTOL_CHANGE,
      NULL);
  }
}

void CplexBackend::ObjAbsTol(ArrayRef<double> val) {
  for (int i = 0; i < (int)val.size(); ++i) {
    CPXmultiobjchgattribs(env(), lp(), i,
      CPX_NO_OFFSET_CHANGE, CPX_NO_WEIGHT_CHANGE,
      CPX_NO_PRIORITY_CHANGE, val[i], CPX_NO_RELTOL_CHANGE,
      NULL);
  }
}

void CplexBackend::ObjRelTol(ArrayRef<double> val) {
  for (int i = 0; i < (int)val.size(); ++i) {
    CPXmultiobjchgattribs(env(), lp(), i,
      CPX_NO_OFFSET_CHANGE, CPX_NO_WEIGHT_CHANGE,
      CPX_NO_PRIORITY_CHANGE, CPX_NO_ABSTOL_CHANGE, val[i],
      NULL);
  }
}
double CplexBackend::Kappa() {
  double kappa;
  int res = CPXgetdblquality(env(), lp(), &kappa, CPX_KAPPA);
  if (res) {
    AddToSolverMessage("Basis condition is unavailable.\n");
    return 0;
  }
  else
    return kappa;
}

SensRanges CplexBackend::GetSensRanges() {
  if (fixed_model_failed_) return {};
  SensRanges sensr;

  std::vector<double> objlow(NumVars()), objhigh(NumVars());
  CPLEX_CALL(CPXobjsa(env(), lp(), 0, NumVars() - 1, objlow.data(), objhigh.data()));
  std::vector<double> rhslow(NumLinCons()), rhshigh(NumLinCons());
  CPLEX_CALL(CPXrhssa(env(), lp(), 0, NumLinCons() - 1, rhslow.data(), rhshigh.data()));
  std::vector<double> ublow(NumVars()), ubhigh(NumVars()), lblow(NumVars()), lbhigh(NumVars());
  CPLEX_CALL(CPXboundsa(env(), lp(), 0, NumVars() - 1, lblow.data(), lbhigh.data(), ublow.data(), ubhigh.data()));

  /// Should postsolve to obtain CON suffixes .sens(lb/ub)(lo/hi)
  /// containing sens ranges for range constraints.
  auto mvlbhi = GetValuePresolver().
    PostsolveGenericDbl({ {lbhigh } });
  auto mvlblo = GetValuePresolver().
    PostsolveGenericDbl({ {lblow} });
  auto mvubhi = GetValuePresolver().
    PostsolveGenericDbl({ {ubhigh} });
  auto mvublo = GetValuePresolver().
    PostsolveGenericDbl({ {ublow} });
  auto mvobjhi = GetValuePresolver().
    PostsolveGenericDbl({ {objhigh} });
  auto mvobjlo = GetValuePresolver().
    PostsolveGenericDbl({ {objlow} });
  auto mvrhshi = GetValuePresolver().
    PostsolveGenericDbl({ {},                   // no var values
                           {{{ CG_Linear,rhshigh }}} });
  auto mvrhslo = GetValuePresolver().
    PostsolveGenericDbl({ {},
                           {{{ CG_Linear, rhslow}}} });

  sensr.varlbhi = mvlbhi.GetVarValues()();
  sensr.varlblo = mvlblo.GetVarValues()();
  sensr.varubhi = mvubhi.GetVarValues()();
  sensr.varublo = mvublo.GetVarValues()();
  sensr.varobjhi = mvobjhi.GetVarValues()();
  sensr.varobjlo = mvobjlo.GetVarValues()();
  sensr.conrhshi = mvrhshi.GetConValues()();
  sensr.conrhslo = mvrhslo.GetConValues()();

  {
    std::vector<double> obj(NumVars());
    CPLEX_CALL(CPXgetobj(env(), lp(), obj.data(), 0, NumVars() - 1));
    auto mv = GetValuePresolver().PostsolveGenericDbl({ { obj } });
    sensr.varobj = mv.GetVarValues()();
  }

  /// We rely on the RangeCon2Slack Converter and its specific
  /// way: adding +slack to the body of the equality constraint.
  /// This reliance is a bit of a hack.
  /// Properly, this should be done in RangeCon2SlkLink

  /// For range constraints, their slacks' sensitivities were
  /// propagated into the following ValueNodes:
  const auto& conlbhi = mvlbhi.GetConValues()();
  const auto& conlblo = mvlblo.GetConValues()();
  const auto& conubhi = mvubhi.GetConValues()();
  const auto& conublo = mvublo.GetConValues()();

  /// We need to convert them for the original constraints' domains:
  std::vector<double> rhs_raw(NumLinCons());
  CPLEX_CALL(CPXgetrhs(env(), lp(), rhs_raw.data(), 0, NumLinCons()-1));
  auto mv_rhs = GetValuePresolver().PostsolveGenericDbl({ {},
                                                      {{{ CG_Linear, rhs_raw }}} });
  const auto& rhs = mv_rhs.GetConValues()();

  std::vector<char> sense_raw(rhs_raw.size());
  auto status = CPXgetsense(env(), lp(), sense_raw.data(), 0, rhs.size()-1);
  assert(!status);
  std::vector<int> sense_raw_int(sense_raw.begin(), sense_raw.end());
  auto mv_sense = GetValuePresolver().PostsolveGenericInt({ {},
                                                            {{{ CG_Linear, sense_raw_int }}} });
  const auto& sense = mv_sense.GetConValues()();
  if (rhs.size() == sensr.conrhshi.size() &&            // check for Release
    !status &&
    rhs.size() == conlbhi.size() &&
    rhs.size() == conlblo.size() &&
    rhs.size() == conubhi.size() &&
    rhs.size() == conublo.size()) {
    sensr.conlbhi = sensr.conlblo = sensr.conubhi = sensr.conublo = rhs;
    for (auto i = rhs.size(); i--; ) {
      if (conubhi[i] != conlblo[i]) {                 // anythin' propagated?
        sensr.conlblo[i] -= conubhi[i];               // then it's a range constraint
        sensr.conlbhi[i] -= conublo[i];
        sensr.conublo[i] -= conlbhi[i];
        sensr.conubhi[i] -= conlblo[i];
      }
      else {                                        // non-range constraints
        sensr.conlblo[i] = -1e100;
        sensr.conlbhi[i] = 1e100;
        sensr.conublo[i] = -1e100;
        sensr.conubhi[i] = 1e100;
        if (sense[i] != (int)'G') {
          sensr.conublo[i] = sensr.conrhslo[i];
          sensr.conubhi[i] = sensr.conrhshi[i];
        }
        if (sense[i] != (int)'L') {
          sensr.conlblo[i] = sensr.conrhslo[i];
          sensr.conlbhi[i] = sensr.conrhshi[i];
        }
      }
    }
  }

  return sensr;
}


std::pair<ArrayRef<double>, ArrayRef<double>> CplexBackend::Sensobj() const
{
  std::vector<double> low(NumVars()), high(NumVars());
  CPLEX_CALL(CPXobjsa(env(), lp(), 0, NumVars() - 1, low.data(), high.data()));
  return { low, high };
}
std::pair<ArrayRef<double>, ArrayRef<double>> CplexBackend::Sensrhs() const {
  std::vector<double> low(NumLinCons()), high(NumLinCons());
  CPLEX_CALL(CPXrhssa(env(), lp(), 0, NumLinCons() - 1, low.data(), high.data()));
  return { low, high };

}
std::pair<std::pair<ArrayRef<double>, ArrayRef<double>>, std::pair<ArrayRef<double>, ArrayRef<double>> > CplexBackend::Sensbounds() const {
  std::vector<double> ublow(NumVars()), ubhigh(NumVars()), lblow(NumVars()), lbhigh(NumVars());
  CPLEX_CALL(CPXboundsa(env(), lp(), 0, NumVars() - 1, lblow.data(), lbhigh.data(), ublow.data(), ubhigh.data()));
  return { {lblow, lbhigh}, {ublow, ubhigh} };
}






ArrayRef<double> CplexBackend::Ray() {
  std::vector<double> ray(NumVars());
  double proof_p;
  CPXgetray(env(), lp(), ray.data());
  auto mv = GetValuePresolver().PostsolveSolution({ ray });
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> CplexBackend::DRay() {
  std::vector<double> dd(NumLinCons());
  if (dd.size() == 0)
    return dd;
  double proof_p;
  CPXdualfarkas(env(), lp(), dd.data(), &proof_p);
  auto vm = GetValuePresolver().PostsolveSolution({
                                               {},
                                               {{{CG_Linear, std::move(dd)}}}
    });
  return vm.GetConValues().MoveOut();        // need the vector itself
}


const CplexBackend::CutInfo CplexBackend::Cut_Info[] = {
    { "Boolean-Quadratic-Polytope",			CPX_CUT_BQP },
    { "Benders",			CPX_CUT_BENDERS },
    { "cover",			CPX_CUT_COVER },
    { "GUB-cover",			CPX_CUT_GUBCOVER },
    { "flow-cover",			CPX_CUT_FLOWCOVER },
    { "clique",			CPX_CUT_CLIQUE },
    { "Gomory",			CPX_CUT_FRAC },
    { "mixed-integer rounding",	CPX_CUT_MIR },
    { "flow-path",			CPX_CUT_FLOWPATH },
    { "disjunctive",		CPX_CUT_DISJ },
    { "implied-bound",		CPX_CUT_IMPLBD },
    { "zero-half",			CPX_CUT_ZEROHALF },
    { "multi-commodity flow",	CPX_CUT_MCF },
    {"lift and project", CPX_CUT_LANDP},
    {"local implied bound", CPX_CUT_LOCALIMPLBD},
    {"RLT",  CPX_CUT_RLT},
    { 0, 0 } };

void CplexBackend::WindupCPLEXSolve() { 
  if (storedOptions_.cutstats_ && IsMIP())
  {
      int j;
      for (auto CI = Cut_Info; CI->cutname; ++CI)
        if (!CPXgetnumcuts(env(), lp(), CI->cuttype, &j) && j > 0)
          AddToSolverMessage(fmt::format("{} {} cut{}\n", j, CI->cutname, j == 1 ? "" : "s"));
  }
}

void CplexBackend::ReportResults() {
  ReportCPLEXResults();
  BaseBackend::ReportResults();
}

void CplexBackend::ReportCPLEXResults() {
  SetStatus( GetSolveResult() );
  AddCPLEXMessages();
  if (need_multiple_solutions())
    ReportCPLEXPool();
  if (need_fixed_MIP())
    ConsiderCplexFixedModel();
  if (!storedOptions_.endBasis_.empty())
    CPLEX_CALL(CPXmbasewrite(env(), lp(), storedOptions_.endBasis_.c_str()));
 
  if (storedOptions_.bestnode_) {
    double obj;
    CPLEX_CALL(CPXgetbestobjval(env(), lp(), &obj));
    std::vector<double> dbl(1, obj);
    ReportSuffix(sufBestNodeObj, dbl);
    ReportSuffix(sufBestNodeProb, dbl);
  }
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

std::pair<int, std::string> CplexBackend::GetSolveResult() {
  namespace sol = mp::sol;
  int optimstatus = CPXgetstat(env(), lp());
  switch (optimstatus) {
  default:
    int solcount;
    solcount = CPXgetsolnpoolnumsolns (env(), lp());  // Can we use it without CPXpopulate?
    if (solcount>0) {
      return { sol::LIMIT_FEAS, "feasible solution" };
    }
    if (interrupter()->Stop()) {
      return { sol::LIMIT_NO_FEAS, "interrupted" };
    }
    return { sol::UNKNOWN, "unknown solution status" };
  case CPX_STAT_OPTIMAL:
  case CPXMIP_OPTIMAL:
  case CPXMIP_OPTIMAL_POPULATED:
  case CPX_STAT_MULTIOBJ_OPTIMAL:
    return { sol::SOLVED, "optimal solution" };
  case CPX_STAT_CONFLICT_FEASIBLE:
    return { sol::SOLVED, "conflict refiner reports the problem as feasible" };
  case CPXMIP_OPTIMAL_TOL:
  case CPXMIP_OPTIMAL_POPULATED_TOL:
    return { sol::SOLVED, "optimal solution within tolerance" };
  case CPXMIP_OPTIMAL_RELAXED_SUM:
  case CPXMIP_OPTIMAL_RELAXED_QUAD:
  case CPX_STAT_OPTIMAL_RELAXED_INF:
  case CPX_STAT_OPTIMAL_RELAXED_QUAD:
  case CPX_STAT_OPTIMAL_RELAXED_SUM:
  case CPXMIP_OPTIMAL_RELAXED_INF:
    return { sol::SOLVED, "optimal solution of relaxed problem" };
  case CPX_STAT_INFEASIBLE:
  case CPXMIP_INFEASIBLE:
  case CPX_STAT_MULTIOBJ_INFEASIBLE:
  case CPX_STAT_CONFLICT_MINIMAL:
    return { sol::INFEASIBLE, "infeasible problem" };
  case CPX_STAT_INForUNBD:
  case CPXMIP_INForUNBD:
  case CPX_STAT_MULTIOBJ_INForUNBD:
    return { sol::LIMIT_INF_UNB, "infeasible or unbounded problem. "
                                "Disable dual reductions "
                                "or run IIS finder for definitive answer." };
  case CPX_STAT_UNBOUNDED:
  case CPXMIP_UNBOUNDED:
  case CPX_STAT_MULTIOBJ_UNBOUNDED:
    solcount = CPXgetsolnpoolnumsolns(env(), lp());  // Can we use it without CPXpopulate?
    if (solcount>0)
      return { sol::UNBOUNDED_FEAS,
            "unbounded problem, feasible solution returned" };
    return { sol::UNBOUNDED_NO_FEAS,
          "unbounded problem, no feasible solution returned" };
  case CPX_STAT_OPTIMAL_FACE_UNBOUNDED:
    return { sol::UNBOUNDED_NO_FEAS, "model has an unbounded optimal face" };
  case CPX_STAT_FEASIBLE:
  case CPXMIP_FEASIBLE:
  case CPX_STAT_FEASIBLE_RELAXED_INF:
  case CPX_STAT_FEASIBLE_RELAXED_QUAD:
  case CPX_STAT_FEASIBLE_RELAXED_SUM:
  case CPXMIP_FEASIBLE_RELAXED_INF:
  case CPXMIP_FEASIBLE_RELAXED_QUAD:
  case CPXMIP_FEASIBLE_RELAXED_SUM:
    return { sol::LIMIT_FEAS, "limit, feasible solution" };
  case CPX_STAT_CONFLICT_ABORT_DETTIME_LIM:
  case CPX_STAT_CONFLICT_ABORT_TIME_LIM:
    return { sol::INFEASIBLE_IIS, "non-minimal IIS due to time limit" };
  case CPX_STAT_CONFLICT_ABORT_MEM_LIM :
    return { sol::INFEASIBLE_IIS, "non-minimal IIS due to memory limit" };
  case  CPX_STAT_CONFLICT_ABORT_NODE_LIM:
    return { sol::INFEASIBLE_IIS, "non-minimal IIS due to node limit" };
  case CPX_STAT_CONFLICT_ABORT_OBJ_LIM:
    return { sol::INFEASIBLE_IIS, "non-minimal IIS due to objective limit" };
  case CPX_STAT_CONFLICT_ABORT_CONTRADICTION:
    return { sol::INFEASIBLE_IIS, "non-minimal IIS due to contradiction" };
  case CPX_STAT_NUM_BEST:
  case CPX_STAT_OPTIMAL_INFEAS:
  case CPXMIP_OPTIMAL_INFEAS:
    return { sol::UNCERTAIN, "reported feasible or optimal but numeric issue" };
  }
}



void CplexBackend::ComputeIIS() {
  int status;
  int cs;
  int nr, nc, nr2, nc2;
  static const char* abort_reason[7] = {
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
      int val = CPX_CONFLICT_EXCLUDED;
      if (lbs[i] > MinusInfinity()) {
        switch (grpStat[j++]) {
        case CPX_CONFLICT_EXCLUDED:
          break;
        case CPX_CONFLICT_MEMBER:
        case CPX_CONFLICT_LB:
          val = CPX_CONFLICT_LB;
          break;
        case CPX_CONFLICT_POSSIBLE_MEMBER:
        case CPX_CONFLICT_POSSIBLE_LB:
          val = CPX_CONFLICT_POSSIBLE_LB;
          break;
        default:
          break;
        }
      }
      if (ubs[i] < Infinity()) {
        switch (grpStat[j++]) {
        case CPX_CONFLICT_EXCLUDED:
          break;
        case CPX_CONFLICT_MEMBER:
        case CPX_CONFLICT_UB:
          val =
              (CPX_CONFLICT_EXCLUDED==val) ?
              CPX_CONFLICT_UB : CPX_CONFLICT_MEMBER;
          break;
        case CPX_CONFLICT_POSSIBLE_MEMBER:
        case CPX_CONFLICT_POSSIBLE_UB:
          val =
              (CPX_CONFLICT_EXCLUDED==val) ?
                  CPX_CONFLICT_POSSIBLE_UB :
                    CPX_CONFLICT_POSSIBLE_MEMBER;
          break;
        default:
          break;
        }
      }
      if (CPX_CONFLICT_EXCLUDED != val) {
        iisColIndices.push_back(i);
        iisColValues.push_back(val);
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
  switch (i) {
  case CPX_CONFLICT_EXCLUDED: return (int)IISStatus::non;
  case CPX_CONFLICT_MEMBER: return (int)IISStatus::mem;
  case CPX_CONFLICT_LB: return (int)IISStatus::low;
  case CPX_CONFLICT_UB: return (int)IISStatus::upp;
  case CPX_CONFLICT_POSSIBLE_MEMBER: return (int)IISStatus::pmem;
  case CPX_CONFLICT_POSSIBLE_LB: return (int)IISStatus::plow;
  case CPX_CONFLICT_POSSIBLE_UB: return (int)IISStatus::pupp;
  }
  return (int)IISStatus::bug;
}
IIS CplexBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}
ArrayRef<int> CplexBackend::VarsIIS() {
  std::vector<int> iis(NumVars(), (int)IISStatus::non);
  for (int i = 0; i < (int)iisColIndices.size(); i++)
    iis[iisColIndices[i]] = (int)IISCplexToAMPL(iisColValues[i]);
  return iis;
}
pre::ValueMapInt CplexBackend::ConsIIS() {
  std::vector<int>
      iis_lincon(NumLinCons(), (int)IISStatus::non),
      iis_qc(NumQPCons(), (int)IISStatus::non),
      iis_indcon(NumIndicatorCons(), (int)IISStatus::non),
      iis_soscon(NumSOSCons(), (int)IISStatus::non);

  for (int i = 0; i < (int)iisRowIndices.size(); i++)
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

void CplexBackend::InputExtras() {
  BaseBackend::InputExtras();
  InputCPLEXExtras();
}

void CplexBackend::InputCPLEXExtras() {
  // Set behaviour for solultion pool related options
  if (!need_multiple_solutions()) {
    storedOptions_.populate_ = -1;
    storedOptions_.poolIntensity_ = -1;
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
  CplexPlayObjNParams();
  SetSolverOption(CPX_PARAM_SOLNPOOLINTENSITY, storedOptions_.poolIntensity_ < 0 ? 0 :
                                                   storedOptions_.poolIntensity_);
}

void CplexBackend::ConsiderCplexFixedModel() {
  if (!IsMIP())
    return;
  MIPResultBeforeFixed result;
  result.primal = PrimalSolution();
  if (result.primal.empty()) {
    AddToSolverMessage("Fixed MIP for mip:basis skipped: no integer-feasible incumbent.");
    mip_result_before_fixed_ = {};
    mip_result_before_fixed_.objective = ObjectiveValue();
    fixed_model_failed_ = true;
    return;
  }
  result.objective = ObjectiveValue();
  result.gap = MIPGap();
  result.absolute_gap = MIPGapAbs();
  result.bound = BestDualBound();
  result.nodes = NodeCount();
  result.simplex_iterations = SimplexIterations();
  result.barrier_iterations = BarrierIterations();
  mip_result_before_fixed_ = std::move(result);
  have_mip_result_before_fixed_ = true;
  // CPXcloneprob copies problem data, not solution data. Save both to
  // retain the MIP result when the supplementary solve reaches a limit.
  int error = 0;
  CPXLPptr original = CPXcloneprob(env(), lp(), &error);
  if (!original) {
    fixed_model_failed_ = true;
    AddToSolverMessage(fmt::format("Cannot save MIP before fixed solve (CPLEX error {}).", error));
    return;
  }
  std::string msg;
  try {
    msg = DoCplexFixedModel();
  } catch (const std::exception& e) {
    msg = e.what();
  }
  if (!msg.empty()) {
    fixed_model_failed_ = true;
    AddToSolverMessage(msg + " failed in DoCplexFixedModel().");
    CPXfreeprob(env(), &lp_ref());
    lp_ref() = original;
    copy_common_info_to_other();
  } else {
    CPXfreeprob(env(), &original);
  }
}


std::string CplexBackend::DoCplexFixedModel() {
  typedef int (CPXPUBLIC* Optalg)(CPXCENVptr, cpxlp*);
  int status;
  Optalg Contopt;

  original_model_type_ = CPXgetprobtype(env(), lp());
  int new_model_type = CPXPROB_FIXEDMILP;
  Contopt = CPXlpopt;
  switch (original_model_type_) {
  case CPXPROB_MIQP:
    new_model_type = CPXPROB_FIXEDMIQP;
    Contopt = CPXqpopt;
    break;
  case CPXPROB_MIQCP: // CPLEX does support changing MIQCPs
    new_model_type = CPXPROB_MIQCP;
    break;
  default:
    new_model_type = CPXPROB_FIXEDMILP;
  }
  if (new_model_type == original_model_type_) // notably MIQCP
  {
    // fix variables manually:
    auto sol = PrimalSolution();
    std::vector<char> types(NumVars());
    CPLEX_CALL(CPXgetctype(env(), lp(), types.data(), 0, NumVars()));
    int nint = CPXgetnumbin(env(), lp()) + CPXgetnumint(env(), lp());
    std::vector<int> indices(nint);
    std::vector<double> bound(nint);
    std::vector<char> type(nint, 'B');
    std::size_t j = 0;
    for (auto i = 0; i < types.size(); i++) {
      if ((types[i] == CPX_BINARY) || (types[i] == CPX_INTEGER))
      {
        bound[j] = sol[i];
        indices[j++] = i;
      }
    }
    CPLEX_CALL(CPXchgbds(env(), lp(), nint, indices.data(), type.data(), bound.data()));
  }
  else
    CPLEX_CALL(CPXchgprobtype(env(), lp(), new_model_type));
  double old_time_limit;
  CPXLONG old_iteration_limit;
  int old_barrier_limit;
  CPLEX_CALL(CPXgetdblparam(env(), CPXPARAM_TimeLimit, &old_time_limit));
  CPLEX_CALL(CPXgetlongparam(env(), CPXPARAM_Simplex_Limits_Iterations, &old_iteration_limit));
  CPLEX_CALL(CPXgetintparam(env(), CPXPARAM_Barrier_Limits_Iteration, &old_barrier_limit));
  auto restore_limits = [&](void*) {
    CPXsetdblparam(env(), CPXPARAM_TimeLimit, old_time_limit);
    CPXsetlongparam(env(), CPXPARAM_Simplex_Limits_Iterations, old_iteration_limit);
    CPXsetintparam(env(), CPXPARAM_Barrier_Limits_Iteration, old_barrier_limit);
  };
  std::unique_ptr<void, decltype(restore_limits)> limits_guard(env(), restore_limits);
  CPLEX_CALL(CPXsetdblparam(env(), CPXPARAM_TimeLimit,
      fixed_model_time_limit() < 0 ? 1e75 : fixed_model_time_limit()));
  CPLEX_CALL(CPXsetlongparam(env(), CPXPARAM_Simplex_Limits_Iterations,
      fixed_model_iter_limit() < 0 ? CPX_BIGLONG : fixed_model_iter_limit()));
  CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_Barrier_Limits_Iteration,
      fixed_model_iter_limit() < 0 ? INT_MAX : fixed_model_iter_limit()));
  CPLEX_CALL(Contopt(env(), lp()));

  int optimstatus = CPXgetstat(env(), lp());
  if (optimstatus != CPX_STAT_OPTIMAL)
    return fmt::format("Fixed model status {}", optimstatus);
  int cnt = CPXgetitcnt(env(), lp());
  AddToSolverMessage(fmt::format("Fixed MIP for mip:basis: {} simplex iteration{}",
    cnt, "s"[cnt == 1.]));
  return {};
}

void CplexBackend::DoCplexFeasRelax() {
  int reltype;
  switch (feasrelax()) {
  case 1:
    reltype = CPX_FEASOPT_MIN_SUM;
    break;
  case 2:
    reltype = CPX_FEASOPT_MIN_QUAD;
    break;
  case 3:
    reltype = CPX_FEASOPT_MIN_INF;
    break;
  case 4:
    reltype = CPX_FEASOPT_OPT_SUM;
    break;
  case 5:
    reltype = CPX_FEASOPT_OPT_QUAD;
    break;
  case 6:
    reltype = CPX_FEASOPT_OPT_INF;
    break;
  default:
    throw std::runtime_error("Unexpected feasrelax value");
  }

  SetCPLEXParam(CPXPARAM_Feasopt_Mode, reltype);
  auto mv = GetValuePresolver().PresolveSolution({
                                              {},
                                              feasrelax().rhspen()
    });
  const auto& rhspen = mv.GetConValues()(CG_Linear);
  std::vector<double> lbpen = feasrelax().lbpen();
  if (lbpen.size() && lbpen.size() < (size_t)NumVars())
    lbpen.resize(NumVars());
  std::vector<double> ubpen = feasrelax().ubpen();
  if (ubpen.size() && ubpen.size() < (size_t)NumVars())
    ubpen.resize(NumVars());
  CPLEX_CALL(CPXfeasopt(env(), lp(),
    (double*)data_or_null(rhspen), (double*)data_or_null(rhspen),
    (double*)data_or_null(lbpen), (double*)data_or_null(ubpen)));
}


const mp::OptionValueInfo values_noautoyes_[3] = {
  {     "-1", "No", -1 },
  {     "0", "Automatic choice (default)", 0 },
  {     "1", "Yes.", 1}
};


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
static const mp::OptionValueInfo optimalitytarget_values_[] = {
    { "0", "automatic (default)", 0},
    { "1", "assume convex and search for global optimum" , 1},
    { "2", "search for first order optimality (not valid for QMIP)" , 2},
    {"3", "solve non-convex to global optimality" ,3 }
};
static const mp::OptionValueInfo outlev_values_[] = {
  { "0", "no output (default)", 0},
  { "1", "equivalent to \"bardisplay\"=1, \"display\"=1, \"mipdisplay\"=2", 1},
  { "2", "equivalent to \"bardisplay\"=2, \"display\"=2, \"mipdisplay\"=5", 2}
};
static const mp::OptionValueInfo auxrootthreads_values_[] = {
  { "-1", "Do not use additional threads for auxiliary tasks", 0},
  { "0", "Automatic (default)", 1},
  { "n < N", "Use n threads for auxiliary root tasks", 2}
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

static const mp::OptionValueInfo values_method[] = {
    { "0", "Automatic (default)",0},
    { "1", "Primal simplex", 1},
    { "2", "Dual simplex", 2},
    { "3", "Network simplex", 3},
    { "4", "Barrier", 4},
    { "5", "Sifting", 5},
    { "6", "Concurrent (dual, barrier and primal in opportunistic mode; "
          "dual and barrier in deterministic mode; 4 is used for MIPQs).", 6}
};

static const mp::OptionValueInfo values_advance_[] = {
    { "0", "Do not use advanced start information", 0},
    { "1", "Use an advanced basis supplied by the user (default)" , 1},
    { "2", "Crush an advanced basis or starting vector supplied by the user.", 2}
};

static const mp::OptionValueInfo values_netopt[] = {
    { "0", "(Default) never invoke the network optimizer", 0},
    { "1", "Compatibility value; same as 3", 1},
    { "2", "Compatibility value; same as 3", 2},
    { "3", "Invoke the network optimizer by setting "
          "CPLEX' LP/QP/MIPNodeMethod to CPX_ALG_NET telling CPLEX "
          "to search for network (sub)structures in the model. "
          "CPLEX presolve might influence automatic recognition of "
          "network structures", 3}
};

static const mp::OptionValueInfo values_netdisplay[] = {
    { "0", "No display", 0},
    { "1", "Display true objective values (can be non-monotonic)", 1},
    { "2", "(Default) display penalized objective values", 2}
};

static const mp::OptionValueInfo values_netfind[] = {
    { "1", "Extract pure network only", 1},
    { "2", "Try reflection scaling (default)", 2},
    { "3", "Try general scaling", 3}
};

static const mp::OptionValueInfo values_netpricing[] = {
    { "0", "(Default) automatic", 0},
    { "1", "Partial pricing", 1},
    { "2", "Multiple partial pricing", 2},
    { "3", "Multiple partial pricing with sorting", 3}
};

static const mp::OptionValueInfo values_nodemethod[] = {
  { "0", "Automatic (default)",0},
  { "1", "Primal simplex", 1},
  { "2", "Dual simplex", 2},
  { "3", "Network simplex", 3},
  { "4", "Barrier", 4},
  { "5", "Sifting", 5}
};
static const mp::OptionValueInfo values_fpheur[] = {
  { "-1", "No",0},
  { "0", "Automatic (default)",0},
  { "1", "Yes, focus on finding a feasible solution", 1},
  { "2", "Yes, focus on finding a good objective value at a feasible solution", 2}
};


static const mp::OptionValueInfo values_dgradient[] = {
  { "0", "Automatic (default)",0},
  { "1", "Standard dual pricing", 1},
  { "2", "Steepest-edge pricing", 2},
  { "3", "Steepest-edge pricing in slack space", 3},
  { "4", "Steepest-edge with unit initial norms", 4},
  { "5", "Devex pricing", 5}
};

static const mp::OptionValueInfo values_benders_strategy[] = {
  { "0", "Automatic (default): if suffix benders is present on variables, "
         "variables that have .benders = 0 go into the master and CPLEX "
				 "assigns other variables to workers; otherwise integer variables "
         "go into the master and variables into workers",0},
  { "1", "Use suffix benders to determine which variables are for the master "
         "(.benders = 0) and which for workers (.benders = n > 0 ==> worker n", 1},
  { "2", "Similar to 0, but suffix benders is required", 2},
  { "3", "Similar to 0, but ignore suffix benders", 3}
};

static const mp::OptionValueInfo values_barcrossover[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "None: return an interior solution", 0},
  { "1", "Primal crossover", 1},
  { "2", "Dual crossover", 2}
};
static const mp::OptionValueInfo values_solutiontype[] = {
  { "0", "Automatic - seeks a solution with basis (default)", 0},
  { "1", "Yes (equivalent to 0)", 1},
  { "2", "No", 2}
};
static const mp::OptionValueInfo values_baralg[] = {
  { "0",  "Default (= 1 for MIP subproblems, else 3)", 0},
  { "1", "Infeasibility-estimate start", 1},
  { "2", "Infeasibility-constant start", 2},
  { "3", "Standard start", 3}
};
static const mp::OptionValueInfo values_barmaxcor[] = {
  { "-1",  "Automatic (default)", -1},
  { "0", "None", 0},
  { "n>0", "Maximum number of centering corrections per iteration", 2}
};
static const mp::OptionValueInfo values_modisplay[] = {
  { "0", "No output", 0},
  { "1", "(Default) summary display after each subproblem", 1},
  { "2", "Summary display after each subproblem, as well as subproblem logs", 2}
};
static const mp::OptionValueInfo values_barstart[] = {
  { "1", "Assume dual is 0 (default)", 1},
  { "2", "Estimate dual", 2},
  { "3", "Average of primal estimate, 0 dual", 3},
  { "4", "Average of primal and dual estimates", 4}
};


static const mp::OptionValueInfo values_pgradient[] = {
  { "-1", "Reduced-cost pricing", -1},
  { "0", "Hybrid reduced-cost and devex pricing (default)", 0},
  { "1", "Devex pricing", 1},
  { "2", "Steepest-edge pricing", 2},
  { "3", "Steepest-edge pricing with slack initial norms", 3},
  { "4", "Full pricing", 4}
};
static const mp::OptionValueInfo values_probe[] = {
  { "-1", "No probing", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Moderate probing", 1},
  { "2", "Aggressive probing", 2},
  { "3", "Very aggressive probing", 3}
};
static const mp::OptionValueInfo values_symmetry[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Moderate symmetry breaking", 1},
  { "2", "Aggressive", 2},
  { "3", "Very aggressive", 3},
  { "4", "Highly aggressive", 4},
  { "5", "Extremely aggressive", 5}
};
static const mp::OptionValueInfo values_repeatpresolve[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Repeat presolve without cuts", 1},
  { "2", "Repeat presolve with cuts", 2},
  { "3", "Repeat presolve with cuts and allow new root cuts", 3}
};
static const mp::OptionValueInfo values_kappastats[] = {
  { "-1", "No MIP kappa statistics", -1},
  { "0", "Automatic (default): compute for a sample of subproblems "
         "when alg:numericalemphasis is on", 0},
  { "1", "Compute for a sample of subproblems", 1},
  { "2", "Compute for all subproblems (possibly expensive)", 2}
};
static const mp::OptionValueInfo values_miqcpstrat[] = {
  { "0", "Automatic (default)", 0},
  { "1", "Solve a QCP relaxation at each node", 1},
  { "2", "Solve an LP relaxation at each node", 2}
};
static const mp::OptionValueInfo values_ordertype[] = {
  { "0", "Do not generate priorities (default)", 0},
  { "1", "Use decreasing costs", 1},
  { "2", "Use increasing bound range", 2},
  { "3", "Use coefficient count", 3}
};
static const mp::OptionValueInfo values_barorder[] = {
  { "0", "Automatic (default)", 0},
  { "1", "Approximate minimum degree (AMD)", 1},
  { "2", "Approximate minimum fill (AMF)", 2},
  { "3", "Nested dissection (ND)", 3}
};
static const mp::OptionValueInfo values_clocktype[] = {
  { "0", "Automatic", 0},
  { "1", "CPU time", 1},
  { "2", "Wall clock time (default)", 2}
};
static const mp::OptionValueInfo values_parallelmode[] = {
  { "-1", "Opportunistic", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Deterministic", 1}
};
static const mp::OptionValueInfo values_datacheck[] = {
  { "0", "Off (default)", 0},
  { "1", "Check data consistency and issue warnings", 1},
  { "2", "Also offer modeling assistance (warnings about bad scaling etc.)", 2}
};
static const mp::OptionValueInfo values_heurfreq[] = {
  { "-1", "Never", -1},
  { "0", "Automatic (default)", 0},
  { "n > 0", "Every n nodes", 1}
};
static const mp::OptionValueInfo values_rinsheur[] = {
  { "-1", "Never", -1},
  { "0", "Automatic (default)", 0},
  { "n > 0", "Every n nodes", 1}
};
static const mp::OptionValueInfo values_repairtries[] = {
  { "-1", "None", -1},
  { "0", "Automatic (default)", 0},
  { "n > 0", "That many times", 1}
};

static const mp::OptionValueInfo values_folding[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Moderate folding", 1},
  { "2", "Aggressive", 2},
  { "3", "Very aggressive", 3},
  { "4", "Highly aggressive", 4},
  { "5", "Extremely aggressive", 5}
};
static const mp::OptionValueInfo values_dive[] = {
  { "0", "Automatic (default)", 0},
  { "1", "Traditional dive (never probing dives)", 1},
  { "2", "Probing dive (look ahead at both children)", 2},
  { "3", "Guided dive (explore solutions similar to the incumbent)", 3}
};
static const mp::OptionValueInfo values_cardls[] = {
  { "-1", "Off", -1},
  { "0", "Automatic (default)", 0},
  { "1, 2", "Apply with increasing effort", 1}
};
static const mp::OptionValueInfo values_qcpduals[] = {
  { "0", "Do not calculate dual values for QCPs", 0},
  { "1", "Calculate them unless this interferes with presolve reductions (default)", 1},
  { "2", "Calculate them, disabling interfering presolve reductions", 2}
};
static const mp::OptionValueInfo values_writelevel[] = {
  { "0", "Automatic (default): all variables in SOL files, "
         "discrete variables in MST files", 0},
  { "1", "All variables", 1},
  { "2", "Discrete variables only", 2},
  { "3", "Variables with nonzero values only", 3},
  { "4", "Discrete variables with nonzero values only", 4}
};
static const mp::OptionValueInfo values_tunemeasure[] = {
  { "1", "Mean time over the models (default)", 1},
  { "2", "Minmax approach", 2}
};
static const mp::OptionValueInfo values_siftalg[] = {
  { "0", "Automatic (default)", 0},
  { "1", "Primal simplex", 1},
  { "2", "Dual simplex", 2},
  { "3", "Network simplex", 3},
  { "4", "Barrier", 4}
};
static const mp::OptionValueInfo values_siftdisplay[] = {
  { "0", "No display", 0},
  { "1", "Display major iterations (default)", 1},
  { "2", "Display LP subproblem information", 2}
};
static const mp::OptionValueInfo values_clonelog[] = {
  { "-1", "No clone logs", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Write clone logs for parallel optimization", 1}
};

static const mp::OptionValueInfo values_presolvenode[] = {
  { "-1", "No node presolve", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Always apply presolve at nodes", 1},
  { "2", "Probing on integer-infeasible variables", 2},
  { "3", "Aggressive node probing", 3}
};

static const mp::OptionValueInfo values_prereduce[] = {
  { "0", "No primal or dual reductions", 0},
  { "1", "Only primal reductions", 1},
  { "2", "Only dual reductions", 2},
  { "3", "Both primal and dual reductions (default)", 3}
};

static const mp::OptionValueInfo values_prereformulations[] = {
  { "0", "None", 0},
  { "1", "Allow reformulations that interfere with crushing forms", 1},
  { "2", "Allow reformulations that interfere with uncrushing forms", 2},
  { "3", "All reformulations (default)", 3}
};

static const mp::OptionValueInfo values_prescale[] = {
  {"-1", "No scaling", -1},
  { "0", "Equilibration (default)", 0},
  { "1", "More aggressive", 1}
};

static const mp::OptionValueInfo values_presosenc[] = {
  {"-1", "No reformulation", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Reformulate as linear constraints, with a reformulation "
         "which is logarithmic in the size of the SOSs.", 1}
};

static const mp::OptionValueInfo values_preaggregator[] = {
  { "-1", "Automatic (default) (1 for LP, infinite for MIP)", 0},
  { "0", "Do not use any aggregator", 1},
  { "n > 0", "Apply aggregator n times", 2}
};
static const mp::OptionValueInfo values_precoeffreduce[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "No", 0},
  { "1", "Reduce only integral coefficients", 1},
  { "2", "Reduce all potential coefficients", 2},
  { "3", "Reduce aggresively with tiling", 3}
}; 
static const mp::OptionValueInfo values_predependency[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "No", 0},
  { "1", "Turn on only at start of preprocessing", 1},
  { "2", "Turn on only at end of preprocessing", 2},
  { "3", "Turn on at both start and end of preprocessing", 3}
}; 
static const mp::OptionValueInfo values_cutpasses[] = {
  { "-1", "None", -1},
  { "0", "Automatic", 0},
  { "n > 0", "Number of passes to perform", 2}
};
static const mp::OptionValueInfo values_cutsfactor[] = {
  { "<0", "No limit", -2},
  { "-1", "No limit (default)", -1},
  { "0 <= n <= 1", "No MIP cuts", 0},
  { "n > 1", "(n-1)*m, where m is the original number of rows(after presolve)", 2}
};
static const mp::OptionValueInfo values_bqpcuts[] = {
  {"-1", "Disallow", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Enable moderate cuts generation", 1},
  { "2", "Enable aggressive cuts generation.", 2},
  { "3", "Enable very aggressive cuts generation.", 3 }
};

static const mp::OptionValueInfo values_cuts2[] = {
  {"-1", "Disallow", -1},
  { "0", "Automatic (default)", 0},
  { "1", "Enable moderate cuts generation", 1},
  { "2", "Enable aggressive cuts generation", 2}
};
static const mp::OptionValueInfo values_branchdir[] = {
  {"-1", "Explore \"down\" branch first", -1},
  { "0", "Explore \"most promising\" branch first (default)", 0},
  { "1", "Explore \"up\" branch first.", 1}
};

static const mp::OptionValueInfo conflictalg_values[] = {
  {"0", "Automatic choice (default)", 0},
  {"1", "Fast", 1},
  {"2", "Propagate", 2},
  {"3", "Presolve", 3},
  {"4", "IIS", 4},
  {"5", "Limited solve", 5},
 {"6", "Full solve", 6}
};

static const mp::OptionValueInfo conflictdisplay_values[] = {
  {"0", "Nothing", 0},
  {"1", "Summary (default)", 1},
  {"2", "Detailed", 2}
};

static const mp::OptionValueInfo values_nodeselect[] = {
  { "0", "Depth-first search", 0},
  { "1", "Breadth-first search (default)", 1},
  { "2", "Best-estimate search", 2},
  { "3", "Alternative best-estimate search", 3 }
};

static const mp::OptionValueInfo values_varselect[] = {
  { "-1", "Branch on variable with smallest integer infeasibility", -1},
  { "0", "Algorithm decides (default)", 0},
  { "1", "Branch on variable with largest integer infeasibility", 1},
  { "2", "Branch based on pseudo costs", 2},
  { "3", "Strong branching", 3 },
  { "4", "Branch based on pseudo reduced costs", 4}
};

static const mp::OptionValueInfo values_mipsearch[] = {
  { "0", "Automatic (default)",0},
  { "1", "Traditional branch and cut", 1},
  { "2", "Dynamic search", 2}
};
  
static const mp::OptionValueInfo values_nodefile[] = {
  { "0", "No", 0},
  { "1", "Compressed node file in memory (default)", 1},
  { "2", "Node file on disk", 2},
  { "3", "Compressed node file on disk", 3}
};

static const mp::OptionValueInfo values_mipemphasis[] = {
  { "0", "Balance finding good feasible solutions and "
    "proving optimality (default)", 0},
  { "1", "Favor finding feasible solutions", 1},
  { "2", "Favor providing optimality", 2},
  { "3", "Focus on improving the best objective bound", 3},
  { "4", "Focus on finding hidden feasible solutions", 4},
  { "5", "Focus on finding high quality solutions earlier (heuristic).", 5}


};
static const mp::OptionValueInfo values_tunedisplay_[] = {
  { "0", "Nothing", 0},
  { "1", "Minimal printing (default)", 1},
  { "2", "Show parameters being tried", 2},
  { "3", "Exhaustive printing", 3}
};
  
void CplexBackend::setSolutionMethod() {
  if (storedOptions_.netopt_>0)         // netopt: lowest priority
    storedOptions_.cpxMethod_ = CPX_ALG_NET;
  int nFlags = bool(storedOptions_.fBarrier_)
    + bool(storedOptions_.fPrimal_)
    + bool(storedOptions_.fDual_)
    + bool(storedOptions_.fBenders_)
    + bool(storedOptions_.fNetwork_)
    + bool(storedOptions_.fSifting_);
  if (nFlags>= 2) 
    AddWarning("Ambiguous solving method",
      "Only one of barrier/primalopt/dualopt/network/sifting/benders should be specified.");
  int flagMethod = CPX_ALG_AUTOMATIC;   // the flags apply to root and nodes
  if (nFlags >= 1)                  // 2nd-lowest priority
  {
    if (storedOptions_.fPrimal_)
      flagMethod = CPX_ALG_PRIMAL;
    if (storedOptions_.fDual_)
      flagMethod = CPX_ALG_DUAL;
    if (storedOptions_.fBarrier_)
      flagMethod = CPX_ALG_BARRIER;
    if (storedOptions_.fNetwork_)
      flagMethod = CPX_ALG_NET;
    if (storedOptions_.fSifting_)
      flagMethod = CPX_ALG_SIFTING;
    storedOptions_.cpxMethod_ = flagMethod;
  }
  if (storedOptions_.algMethod_)   // Highest priority
    storedOptions_.cpxMethod_ = storedOptions_.algMethod_;
  auto nodeMethod = flagMethod;
  if (storedOptions_.nodeMethod_)  // Highest priority
    nodeMethod = storedOptions_.nodeMethod_;
  if (IsMIP()) {
    SetSolverOption(CPX_PARAM_STARTALG, storedOptions_.cpxMethod_);
    SetSolverOption(CPX_PARAM_SUBALG, nodeMethod);
  }
  else if (IsQP())
    SetSolverOption(CPX_PARAM_QPMETHOD, storedOptions_.cpxMethod_);
  else
    SetSolverOption(CPX_PARAM_LPMETHOD, storedOptions_.cpxMethod_);

	if (storedOptions_.dualprob_)
		storedOptions_.predual_ *= -1;
	if (storedOptions_.predual_)
		SetSolverOption(CPXPARAM_Preprocessing_Dual, storedOptions_.predual_);

  if (storedOptions_.cpxMethod_ == CPX_ALG_BARRIER) {
    SetSolverOption(CPX_PARAM_SOLUTIONTYPE, storedOptions_.solutionType_);
    if (storedOptions_.crossover_ == -1) // translate from MP to cplex
      SetSolverOption(CPXPARAM_Barrier_Crossover, 0);
    else if ((storedOptions_.crossover_ >= 1) && (storedOptions_.crossover_ <= 2))
      SetSolverOption(CPXPARAM_Barrier_Crossover, storedOptions_.crossover_);
    else if (storedOptions_.crossover_ == 0)
      SetSolverOption(CPX_PARAM_SOLUTIONTYPE, CPX_NONBASIC_SOLN);
    else AddToSolverMessage("Invalid crossover value specified, see -= output");
  }

  if (storedOptions_.fBenders_)             // at last
  {
    storedOptions_.cpxMethod_ = CPX_ALG_BENDERS;
    ReadBendersSuffix();
  }
}

////////////////////////////// OPTIONS /////////////////////////////////

void CplexBackend::FinishOptionParsing() {
  {
    // Set output on screen
    int lp, mip, bar;
    GetSolverOption(CPX_PARAM_SIMDISPLAY, lp);
    GetSolverOption(CPX_PARAM_MIPDISPLAY, mip);
    GetSolverOption(CPX_PARAM_BARDISPLAY, bar);
    // GetSolverOption(CPXPARAM_MultiObjective_Display, mo);
    // GetSolverOption(CPX_PARAM_NETDISPLAY, netw);
    if (storedOptions_.outlev_ < 0)
      storedOptions_.outlev_ = 0;
    if (storedOptions_.outlev_ > 2)
      storedOptions_.outlev_ = 2;
    int olp[] = { 0, 1, 2 };
    int omip[] = { 0, 2, 5 };
    lp = lp ? lp : olp[storedOptions_.outlev_];
    mip = mip ? mip : omip[storedOptions_.outlev_];
    bar = bar ? bar : olp[storedOptions_.outlev_];
    // mo = mo ? mo : olp[storedOptions_.outlev_];  -- no they are >0
    // netw = netw ? netw : olp[storedOptions_.outlev_];
    // Set options before turning on output:
    if (lp)
      SetSolverOption(CPX_PARAM_SIMDISPLAY, lp);
    if (mip)
      SetSolverOption(CPX_PARAM_MIPDISPLAY, mip);
    if (bar)
      SetSolverOption(CPX_PARAM_BARDISPLAY, bar);
    // SetSolverOption(CPXPARAM_MultiObjective_Display, mo);
    // SetSolverOption(CPX_PARAM_NETDISPLAY, netw);
    if (lp || mip || bar
        // || mo || netw  -- no they are >0
        ) {
      if (!storedOptions_.outlev_)
        storedOptions_.outlev_ = 1;      // for setSilenceOutput()
      /* Log messages on screen */
      CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ScreenOutput, CPX_ON));
      /* Echo changed params before solve */
      CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ParamDisplay, 1));
      set_verbose_mode(true);
    }
    if (!storedOptions_.logFile_.empty())
    {
      if (lp < 1) SetSolverOption(CPX_PARAM_SIMDISPLAY, 1);
      if (mip < 1) SetSolverOption(CPX_PARAM_MIPDISPLAY, 1);
      if (bar < 1) SetSolverOption(CPX_PARAM_BARDISPLAY, 1);
      CPLEX_CALL(CPXsetlogfilename(env(), storedOptions_.logFile_.data(), "w"));
    }
    set_verbose_mode(storedOptions_.outlev_ > 0);
  }

  // Global cut control (cut:cuts): apply to every cut family
  // which the user has not set individually (i.e., still at default 0).
  // Families whose maximum is 2 receive min(cuts, 2).
  if (storedOptions_.cuts_ != -2) {
    static const int op[] = {
      CPXPARAM_MIP_Cuts_BQP,
      CPXPARAM_MIP_Cuts_Cliques,
      CPXPARAM_MIP_Cuts_Covers,
      CPXPARAM_MIP_Cuts_Disjunctive,
      CPXPARAM_MIP_Cuts_FlowCovers,
      CPXPARAM_MIP_Cuts_Gomory,
      CPXPARAM_MIP_Cuts_GUBCovers,
      CPXPARAM_MIP_Cuts_Implied,
      CPXPARAM_MIP_Cuts_LiftProj,
      CPXPARAM_MIP_Cuts_LocalImplied,
      CPXPARAM_MIP_Cuts_MCFCut,
      CPXPARAM_MIP_Cuts_MIRCut,
      CPXPARAM_MIP_Cuts_Nodecuts,
      CPXPARAM_MIP_Cuts_PathCut,
      CPXPARAM_MIP_Cuts_RLT,
      CPXPARAM_MIP_Cuts_ZeroHalfCut
    };
    for (size_t f = 0; f < sizeof(op) / sizeof(int); f++) {
      int local = 0;
      CPXgetintparam(env(), op[f], &local);
      if (0 == local) {
        CPXINT dflt = 0, vmin = -1, vmax = 3;
        CPXinfointparam(env(), op[f], &dflt, &vmin, &vmax);
        SetSolverOption(op[f], std::min(storedOptions_.cuts_, (int)vmax));
      }
    }
  }

  if (!storedOptions_.mipStart_.empty())
  {
    // Write MIP starts
    auto nms = CPXgetnummipstarts(env(), lp());
    if (nms == 0)
      AddToSolverMessage("No MIP start is available.\n");
    CPLEX_CALL(CPXwritemipstarts(env(), lp(),
                                 storedOptions_.mipStart_.c_str(), 0, nms - 1));
  }

  if (!storedOptions_.cpuMask_.empty())
    CPLEX_CALL(CPXsetstrparam(env(),
                              CPXPARAM_CPUmask, storedOptions_.cpuMask_.c_str()));

  if (storedOptions_.numcores_) {
    CPLEX_CALL(CPXgetnumcores(env(), &storedOptions_.numcores_));
    AddToSolverMessage(
        fmt::format("{} logical cores are available.\n", storedOptions_.numcores_));
  }

  if (!storedOptions_.workDir_.empty())
    CPLEX_CALL(CPXsetstrparam(env(),
                              CPX_PARAM_WORKDIR, storedOptions_.workDir_.c_str()));

  // Native params
  if (storedOptions_.paramRead_.size())
    CPLEX_CALL(
        CPXreadcopyparam(env(), storedOptions_.paramRead_.c_str() ));
  // Write native params
  if (storedOptions_.paramWrite_.size())
    CPLEX_CALL(
        CPXwriteparam(env(), storedOptions_.paramWrite_.c_str() ));

}
void CplexBackend::InitCustomOptions() {

  set_option_header(
      "IBM ILOG CPLEX Optimizer Options for AMPL\n"
      "-----------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cplex_options``. For example::\n"
      "\n"
      "  ampl: option cplex_options 'mipgap=1e-6';\n");
  

  AddToOptionDescription(
      "tech:writesolution",
      "File name extensions can be "
      "``.sol[.gz|.zip]``, ``.clp``, ``.bas``.");



  // Multi objective controls
  AddIntOption("obj:*:priority obj_*_priority", "Priority for objective with index *",
    &CplexBackend::CplexGetObjIntParam, &CplexBackend::CplexSetObjIntParam);

  AddDblOption("obj:*:weight obj_*_weight", "Weight for objective with index *",
    &CplexBackend::CplexGetObjDblParam, &CplexBackend::CplexSetObjDblParam);

  AddDblOption("obj:*:reltol obj_*_reltol", "Relative tolerance for objective with index *",
    &CplexBackend::CplexGetObjDblParam, &CplexBackend::CplexSetObjDblParam);

  AddDblOption("obj:*:abstol obj_*_abstol", "Absolute tolerance for objective with index *. "
    "Can only be applied on a multi-objective problem with obj:multi=1",
    &CplexBackend::CplexGetObjDblParam, &CplexBackend::CplexSetObjDblParam);




  // Solution method
  AddStoredOption("alg:method method lpmethod qpmethod simplex mip:method mipstartalg startalg startalgorithm",
                  "Which algorithm to use for non-MIP problems or for the root node "
                  "of MIP problems:\n"
                  "\n.. value-table::\n"
                  "For MIQP problems (quadratic objective, linear constraints), setting 5 "
                  "is treated as 0 and 6 as 4. For MIQCP problems (quadratic objective & "
                  "constraints), all settings are treated as 4. See also mip:nodemethod."
                  "\n"
                  "Overrides netopt option and "
                  "primalopt/dualopt/barrier/network/sifting flags.",
                  storedOptions_.algMethod_, values_method);

  AddToOptionDescription("alg:start",
                         "For IBM ILOG CPLEX, choices can be refined via alg:start:advance.");

  AddCplexParamOption("alg:start:advance alg:start:advanced advance advind",
                  "If set to 1 or 2, this parameter specifies that CPLEX should use advanced "
                  "starting information when it initiates optimization. "
                  "Exact meaning differs by problem type. Here for LP:"
                  "\n\n.. value-table::",
                  CPXPARAM_Advance, values_advance_, 1);

  AddStoredOption("alg:barrier barrier baropt",
    "Solve (MIP node) LP/QPs by barrier method.",
    storedOptions_.fBarrier_);

	AddStoredOption("alg:primalopt alg:primal primalopt optimize",
    "Solve (MIP node) LPs by primal simplex method.",
    storedOptions_.fPrimal_);

	AddStoredOption("alg:dualopt alg:dual dualopt tranopt",
    "Solve (MIP node) LPs by dual simplex method.",
    storedOptions_.fDual_);

  AddStoredOption("alg:sifting sifting siftopt siftingopt",
    "Solve (MIP node) LPs by sifting method.",
    storedOptions_.fSifting_);

  AddStoredOption("alg:network network",
    "Solve (substructure of) (MIP node's) LP/QP "
    "by network simplex method. Synonym for alg:netopt=3.",
    storedOptions_.fNetwork_);

  AddStoredOption("alg:netopt netopt",
                  "Whether to use network simplex method for non-MIP problems "
                  "or for the continuous relaxations of MIP nodes, unless "
                  "alg:(node)method or "
                  "primalopt/dualopt/barrier/network/sifting flags are specified. "
                  "Options alg:(node)method override (for MIP, in the root or "
                  "subnodes):\n"
                  "\n.. value-table::\n",
                  storedOptions_.netopt_, values_netopt);

  AddCplexParamOption("tech:netdisplay netdisplay",
                  "Decides what CPLEX reports to the screen during network optimization:\n"
                  "\n.. value-table::\n",
                  CPXPARAM_Network_Display,
                  values_netdisplay, 2);

  AddCplexParamOption("alg:netfeasibility netfeasibility netfeastol",
                  "Feasibility tolerance for network primal optimization. "
                  "Can be any number from 1e-11 to 1e-1; default: 1e-6.",
                  CPXPARAM_Network_Tolerances_Feasibility,
                  1e-11, 1e-1);

  AddCplexParamOption("alg:netfind netfind netfinder",
                  "Level of network extraction for network simplex optimization:\n"
                  "\n.. value-table::\n",
                  CPXPARAM_Network_NetFind,
                  values_netfind, 2);

  AddCplexParamOption("alg:netoptimality netoptimality",
                  "Specifies the optimality tolerance for network optimization; "
                  "that is, the amount a reduced cost may violate the criterion "
                  "for an optimal solution. "
                  "Can be any number from 1e-11 to 1e-1; default: 1e-6.",
                  CPXPARAM_Network_Tolerances_Optimality,
                  1e-11, 1e-1);

  AddCplexParamOption("alg:netpricing netpricing",
                  "Pricing algorithm for network simplex optimization:\n"
                  "\n.. value-table::\n",
                  CPXPARAM_Network_Pricing,
                  values_netpricing, 0);


  AddStoredOption("alg:benders benders bendersopt",
    "Solve MIP using Benders algorithm. Both integer and continuous "
    "variables must be present.",
    storedOptions_.fBenders_);
  

  AddCplexParamOption("alg:benders_worker benders_worker",
    "Designate the algorithm that CPLEX applies to solve the "
    "subproblems when using Benders decomposition:\n"
    "\n.. value-table::\n",
    CPXPARAM_Benders_WorkerAlgorithm,
    values_nodemethod, 0);

  AddCplexParamOption("alg:benders_feascuttol benders_feascut_tol",
    "Tolerance for violations of feasibility cuts in Benders "
    "algorithm (default 1e-6; must be in [1e-9, 0.1]).",
    CPXPARAM_Benders_Tolerances_FeasibilityCut,
    1e-9, 0.1);

  AddCplexParamOption("alg:benders_optcut_tol benders_optcut_tol",
    "Tolerance for violations of optimality cuts in Benders "
    "algorithm (default 1e-6; must be in [1e-9, 0.1]).",
    CPXPARAM_Benders_Tolerances_OptimalityCut,
    1e-9, 0.1);

  AddCplexParamOption("alg:benders_strategy benders_strategy",
    "How to decompose the problem for Benders algorithm:\n"
    "\n.. value-table::\n",
    CPXPARAM_Benders_Strategy,
    values_benders_strategy, 0);

  AddCplexParamOption("alg:feastol feastol feasibility",
    "Primal feasibility tolerance (default 1e-6, possible "
    "values are between 1e-9 and 0.1).",
    CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9, 1e-1);




  AddCplexParamOption("lp:crash crash",
    "Crash strategy (used to obtain starting basis in simplex); "
    "possible values = -1, 0, 1; default = 1. "
    "The best setting is problem-dependent and can only be found "
    "by experimentation. 0 completely ignores the objective.",
    CPXPARAM_Simplex_Crash, -1, 1);

  AddCplexParamOption("lp:dgradient dgradient",
    "Pricing algorithm for dual simplex (default 0):\n"
    "\n.. value-table::\n",
    CPXPARAM_Simplex_DGradient, 
    values_dgradient, 0);

  AddCplexParamOption("lp:doperturb doperturb perturb",
    "Whether to perturb the problem initially when using simplex "
    "(by the amount given by lp:perturb); setting it "
    "to 1 is occasionally helpful for highly degenerate problems:\n"
    "\n.. value-table::\n",
    CPXPARAM_Simplex_Perturbation_Indicator,
    values_01_noyes_0default_, 0);

  AddCplexParamOption("lp:perturblim perturblim perturblimit limitperturb",
    "Number of stalled simplex iterations before the "
    "problem is perturbed; default=0 => automatic.",
    CPXPARAM_Simplex_Limits_Perturbation, 0, CPXINT_MAX);

  AddCplexParamOption("lp:perturb perturbation perturbconst",
    "Amount by which to perturb variable bounds when "
    "perturbing problems (see \"lp:doperturb\"); default=1e-6, "
    "must be in [1e-8, 1e75].",
    CPXPARAM_Simplex_Perturbation_Constant, 1e-8, DBL_MAX);

  AddCplexParamOption("lp:opttol optimality opttol optimalitytolerance",
    "Dual optimality tolerance: simplex will stop when all "
    "reduced costs are smaller than this value in the improving "
    "direction; default=1e-6.",
    CPXPARAM_Simplex_Tolerances_Optimality, 1e-9, 1e-1);


  // Cut generation
  
  // Hidden param
  #define CPX_PARAM_EPSAGG 1055
  AddCplexParamOption("cut:aggtol aggtol",
    "Pivot tolerance for aggregating (undocumented CPLEX parameter 1055).  "
    "It seldom needs fiddling.  Default = .05; must be in [1e-10, .99].",
    CPX_PARAM_EPSAGG, 1e-10, .99);

  AddCplexParamOption("cut:aggforcut aggforcut aggcutlim",
    "Bound on the number of constraints aggregated to generate flow-cover "
    "and mixed-integer-rounding cuts (default 3).",
    CPXPARAM_MIP_Limits_AggForCut, 0, CPXINT_MAX);


  AddCplexParamOption("cut:bqp bqpcuts",
    "Whether to enable Boolean Quadric Polytope cut generation for nonconvex QP and "
    "MIQP problems, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_BQP, values_bqpcuts, 0);

  AddCplexParamOption("cut:clique cliquecuts cliques",
    "Whether to use clique cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Cliques, values_bqpcuts, 0);

  AddCplexParamOption("cut:cover covercuts covers",
    "Whether to use cover cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Covers, values_bqpcuts, 0);

  AddCplexParamOption("cut:disj disjcuts",
    "Whether to use disjunctive cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Disjunctive, values_bqpcuts, 0);

  AddCplexParamOption("cut:eachcutlim eachcutlim",
    "Limit on the number of cuts of each type; default=2100000000.", 
    CPXPARAM_MIP_Limits_EachCutLimit, 0, CPXINT_MAX);

  AddCplexParamOption("cut:flowcover flowcovercuts flowcuts",
    "Whether to use flowcover cuts in solving MIPs:\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_Cuts_FlowCovers, values_cuts2, 0);


  AddCplexParamOption("cut:flowpath flowpathcuts",
    "Whether to use flow path cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_PathCut, values_cuts2, 0);

  AddCplexParamOption("cut:frac gomory fraccuts",
    "Whether to use fractional Gomory cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_Gomory, values_cuts2, 0);

  AddCplexParamOption("cut:fracpass fracpass gomorypass fractionalcuts",
    "Limit on number of passes to generate MIP fractional Gomory cuts, default=0 => automatic choice, "
    "positive = at most that many passes.",
    CPXPARAM_MIP_Limits_GomoryPass, 0, CPXINT_MAX);

  AddCplexParamOption("cut:fraccand fraccand gomorycand",
    "Limit on number of candidate variables when generating Gomory cuts "
    "for MIP problems; default=200.",
    CPXPARAM_MIP_Limits_GomoryCand, 1, CPXINT_MAX);


  AddCplexParamOption("cut:gubcover gubcover gubcuts",
    "Whether to use fractional GUB cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_GUBCovers, values_cuts2, 0);

  AddCplexParamOption("cut:implied implied impliedcuts",
    "Whether to use implied cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_Implied, values_cuts2, 0);

  AddCplexParamOption("cut:liftandproject liftandproject liftandprojectcuts splitcuts",
    "Whether to use lift and project cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_LiftProj, values_bqpcuts, 0);

  AddCplexParamOption("cut:localimplied localimplied localimpliedcuts",
    "Whether to use locally valid implied cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_LocalImplied, values_bqpcuts, 0);

  AddCplexParamOption("cut:mcf mcfcuts cut:mfc mfccuts",
    "Whether to use multi-commodity flow (MCF) cuts in solving MIPs, choices as for \"cut:flowcover\". "
    "(cut:mfc/mfccuts are legacy spellings.)",
    CPXPARAM_MIP_Cuts_MCFCut, values_cuts2, 0);

  AddCplexParamOption("cut:mir mircuts",
    "Whether to use mixed-integer rounding (MIR) cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_MIRCut, values_cuts2, 0);

  AddCplexParamOption("cut:node nodecuts",
    "Decides whether or not cutting planes are separated at the nodes of the branch-and-bound tree, "
    "choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Nodecuts, values_bqpcuts, 0);

  AddCplexParamOption("cut:rlt rltcuts",
    "Whether to use the Relaxation Linearization Technique (RLT) to generate cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_RLT, values_bqpcuts, 0);

  AddCplexParamOption("cut:zerohalf zerohalfcuts",
    "Whether to use zero-half cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_ZeroHalfCut, values_cuts2, 0);

  AddStoredOption("cut:cuts cuts mipcuts",
    "Global cut generation control, applied to all cut families "
    "(cut:bqp, cut:clique, cut:cover, cut:disj, cut:flowcover, cut:frac, "
    "cut:gubcover, cut:implied, cut:liftandproject, cut:localimplied, "
    "cut:mcf, cut:mir, cut:node, cut:flowpath, cut:rlt, cut:zerohalf) "
    "unless overridden by the individual cut-type controls; "
    "value 3 is reduced to 2 for families accepting at most 2:\n"
    "\n.. value-table::\n",
    storedOptions_.cuts_, values_bqpcuts);

  AddCplexParamOption("cut:factor cutsfactor cutfactor",
    " Limit on MIP cuts added:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Limits_CutsFactor, values_cutsfactor, 0);

  AddCplexParamOption("cut:passes cutpasses cutpass",
    "Maximum number of cutting-plane passes during root-cut generation:\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_Limits_CutPasses, values_cutpasses, 0);

  AddStoredOption("cut:stats cutstats",
    " Whether the solve_message report the numbers and kinds of cuts used.\n"
    "\n.. value-table::\n",
    storedOptions_.cutstats_, values_01_noyes_0default_);

  // MIP related options

  AddCplexParamOption("mip:backtrack backtrack",
    "Tolerance (>0, default 0.9999) for when to backtrack during "
    "branch & bound.  Low values tend to pure best-bound search. "
    "High values(~1) tend to pure depth-first search. Values less "
    "than the default are often good when subproblems are expensive.",
    CPXPARAM_MIP_Strategy_Backtrack, 0.0, 1.0);

  AddCplexParamOption("mip:bbinterval bbinterval",
    "For nodeselect = 2, select the best-bound node, rather than the "
    "best - estimate node, every bbinterval iterations (default 7); "
    "0 means always use the best-estimate node.",
    CPXPARAM_MIP_Strategy_BBInterval, 0, CPXINT_MAX);

  AddStoredOption("mip:bestnode bestnode",
    "Request return of suffix .bestnode on the objective and problem for the "
    "objective value at the best feasible MIP node. "
    "If a feasible node has not yet been found, this value "
    "is +Infinity for minimization problems and -Infinity "
    "for maximization problems.",
    storedOptions_.bestnode_, 0, 1);

  AddCplexParamOption("mip:boundstr boundstr bndstrenind",
    "Whether to apply bound strengthening in mixed integer programs (MIPs):\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_BoundStrength, 
    values_autonoyes_, -1);

  AddCplexParamOption("mip:branchdir branchdir branch",
    "Which child node to explore first when branching:\n"
    "\n.. value-table::",
    CPXPARAM_MIP_Strategy_Branch, values_branchdir, 0);

  AddCplexParamOption("mip:emphasis mipemphasis mip:focus mipfocus",
    "MIP solution strategy:\n" "\n.. value-table::\n",
    CPXPARAM_Emphasis_MIP, values_mipemphasis, 0);


  // TECH:

  AddStoredOption("tech:optionnativeread optionnativeread tech:param:read param:read",
                  "Name of CPLEX parameter file (surrounded by 'single' or "
                  "\"double\" quotes if the name contains blanks). "
                  "The easiest approach to use PRM files is to export an initial file "
                  "with tech:optionnativewrite.",
                  storedOptions_.paramRead_);
  AddStoredOption("tech:optionnativewrite optionnativewrite tech:param:write param:write",
                  "Name of CPLEX parameter file (surrounded by 'single' or \"double\" quotes if the "
                  "name contains blanks) to be written. Only parameters with non-default "
                  "values are written.",
                  storedOptions_.paramWrite_);

  AddCplexParamOption("tech:memoryemphasis memoryemphasis",
                  "0*/1: Whether to compress data to reduce the memory used, "
                  "which may make some information (e.g., basis condition) "
                  "unavailable.",
                  CPXPARAM_Emphasis_Memory, 0, 1);

  AddCplexParamOption(
      "alg:numericalemphasis numfocus numericemphasis numericalemphasis "
      "alg:numericfocus numericfocus",
      "0*/1: Whether to try to improve numerical accuracy (at a "
      "possible cost of time or memory).",
      CPXPARAM_Emphasis_Numerical, 0, 1);


  AddCplexParamOption("mip:gapabs mipgapabs absmipgap",
    "Max. absolute MIP optimality gap (default 1e-6; must be in [0, 1e20]).",
    CPXPARAM_MIP_Tolerances_AbsMIPGap, 0.0, 1e20);

  AddCplexParamOption("mip:gap mipgap",
    "Max. relative MIP optimality gap (default 1e-4; must be in [0, 1]).",
    CPXPARAM_MIP_Tolerances_MIPGap, 0.0, 1.0);


  AddCplexParamOption("mip:fpheur fpheur",
    "Whether to use the feasibility pump heuristic on MIP "
    "problems:\n\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_FPHeur, values_fpheur, 0);

  AddCplexParamOption("mip:inttol inttol intfeastol integrality",
    "Feasibility tolerance for integer variables (default 1e-05, "
    "must be in [0.0, 0.5])",
    CPXPARAM_MIP_Tolerances_Integrality, 0.0, 0.5);


  AddCplexParamOption("mip:nodefile nodefile",
    "Whether to save node information in a temporary file:\n"
    "\n.. value-table::\n", CPXPARAM_MIP_Strategy_File, 
    values_nodefile, 0);



  AddCplexParamOption("mip:search mipsearch",
    "Search strategy for mixed-integer problems:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_Search, values_mipsearch, 0);
  

  AddCplexParamOption("mip:submipalg submipalg",
    "Choice of algorithm used to solve the subproblems of a subMIP: "
    "not a subproblem, but an auxiliary MIP that CPLEX sometimes forms "
    "and solves, e.g., when dealing with a partial MIP start, "
    "repairing an infeasible MIP start, using the RINS heuristic, "
    "branching locally or polishing a solution. "
    "Possible values (only 0 is allowed for MIQCPs):\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_SubMIP_SubAlg, values_nodemethod, 0);

  AddCplexParamOption("mip:submipscale submipscale",
    "Rarely used choice of scaling for auxiliary subMIPs "
    "(described with \"submipalg\"):\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_SubMIP_Scale, values_prescale, 0);


  AddCplexParamOption("mip:submipstartalg submipstartalg",
    "Rarely used choice of algorithm for the initial relaxation of a "
    "subMIP (described with \"submipalg\"); for MIQPs, 5 is treated as 0; "
    "only 0 is allowed for MIQCPs:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_SubMIP_StartAlg, values_nodemethod, 0);


  AddStoredOption("mip:nodemethod nodemethod mipalg mipalgorithm mipsubalg subalg subalgorithm",
    "Algorithm used to solve relaxed MIP node problems "
                  "(after the initial relaxation); for MIQP problems "
    "(quadratic objective, linear constraints), settings other than 3 and 5 " 
    "are treated as 0. For MIQCP problems (quadratic objective and "
		"constraints), only 0 is permitted.\n"
    "\n.. value-table::\n"
                  "\n"
                  "See also alg:method. Overrides netopt option and "
                  "primalopt/dualopt/barrier/network/sifting flags.",
                  storedOptions_.nodeMethod_, values_nodemethod);

  AddCplexParamOption("mip:nodesel nodesel nodeselect",
    "Strategy for choosing next node while optimizing\n\
		integer variables:\n\n.. value-table::",
    CPXPARAM_MIP_Strategy_NodeSelect, values_nodeselect, 0);

  AddCplexParamOption("mip:varbranch varbranch varsel varselect" ,
    "MIP branch variable selection strategy:\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_Strategy_VariableSelect, values_varselect, 0);

  AddCplexParamOption("bar:baralg baralg",
    "How to start the barrier algorithm:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_Algorithm, values_baralg, 0);

  AddCplexParamOption("bar:comptol comptol",
    "Convergence tolerance for barrier algorithm: "
		"the algorithm stops when the relative "
		"complementarity is < bartol (default 1e-8).",
    CPXPARAM_Barrier_ConvergeTol, 1e-12, DBL_MAX);

  AddStoredOption("bar:crossover crossover mipcrossover",
    "How to transform a barrier solution to a basic one:\n"
   "\n.. value-table::\n", storedOptions_.crossover_, values_barcrossover);

  AddCplexParamOption("bar:corrections barcorr barmaxcor",
    "Limit on centering corrections in each iteration of the barrier algorithm:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_Limits_Corrections, values_barmaxcor, -1);

  AddCplexParamOption("bar:dense bar:densecol dense densecol",
    "If positive, minimum nonzeros in a column for the barrier algorithm "
    "to consider the column dense. If 0 (default), this tolerance is "
    "selected automatically.", 
    CPXPARAM_Barrier_ColNonzeros, 0, INT_MAX);


  AddCplexParamOption("tech:modisplay modisplay multiobjdisplay",
    "Level of display during multiobjective optimization:"
    "\n\n.. value-table::\n",
    CPXPARAM_MultiObjective_Display, values_modisplay, 1);

  AddCplexParamOption("bar:growth bargrowth growth",
    "Tolerance for detecting unbounded faces in the barrier algorithm: "
    "higher values make the test for unbounded faces harder to satisfy "
    "(default 1e12; must be >= 1).",
    CPXPARAM_Barrier_Limits_Growth, 1.0, 1e75);

  AddCplexParamOption("bar:iterlim bariterlim",
    "Maximum barrier iterations allowed (default: no limit).",
    CPXPARAM_Barrier_Limits_Iteration, 0, CPXINT_MAX);

  AddCplexParamOption("bar:objrange barobjrange",
    "Limit on the absolute objective value before the barrier algorithm "
    "considers the problem unbounded (default 1e20).",
    CPXPARAM_Barrier_Limits_ObjRange, 0.0, 1e75);

  AddCplexParamOption("bar:start barstart barstartalg",
    "Barrier starting-point algorithm:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_StartAlg, values_barstart, 1);

  AddStoredOption("lp:solutiontype solutiontype",
    "Whether to seek a basic solution when solving an LP:\n"
    "\n.. value-table::\n", storedOptions_.solutionType_, values_solutiontype);

  // Solution pool controls
  AddCplexParamOption("sol:poolgap ams_eps poolgap",
    "Relative tolerance for reporting alternate MIP solutions "
    "(default: 1e75).", CPX_PARAM_SOLNPOOLGAP, 0.0, DBL_MAX);
  AddCplexParamOption("sol:poolgapabs ams_epsabs poolagap poolgapabs",
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
  AddOptionSynonyms_Inline_Front("poolstub", "sol:stub");

  AddCplexParamOption("sol:poollimit ams_limit poolcapacity poollimit solnlimit",
    "Limit on the number of alternate MIP solutions written. Default: 2100000000.",
    CPX_PARAM_SOLNPOOLCAPACITY, 1, 2100000000);

  AddCplexParamOption("sol:poolpopulatelim populatelim",
    "Limit on number of solutions added to the solution pool by the populate algorithm. "
    "Default: 20.",
    CPX_PARAM_POPULATELIM, 1, CPXINT_MAX);

  AddCplexParamOption("sol:poolreplace poolreplace",
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

  AddCplexParamOption("tech:auxrootthreads auxrootthreads",
    "Controls the number of threads used for auxiliary chores when solving "
    "the root node of a MIP problem. When N threads are available (possibly "
    "limited by \"threads\"), auxrootthreads must be less than N.\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Limits_AuxRootThreads, auxrootthreads_values_, 0);

  
  AddCplexParamOption("alg:conflictalg conflictalg",
    "Choice of algorithm used by the CPLEX's conflict "
		"refiner:\n"
    "\n.. value-table::\n"
    "\nSettings 1, 2, and 3 are fast but may not discard "
    "many constraints; 5 and 6 work harder at this. "
    "Setting 4 searches for an Irreducible Infeasible "
    "Set of linear constraints(e.g., ignoring quadratic "
    "constraints).",
    CPXPARAM_Conflict_Algorithm, conflictalg_values, 0);

  AddCplexParamOption("tech:conflictdisplay conflictdisplay",
    "What to report when the conflict finder is working:\n"
    "\n.. value-table::\n",
    CPXPARAM_Conflict_Display, conflictdisplay_values, 0);

  AddStoredOption("tech:cpumask cpumask",
    "Whether and how to bind threads to cores on systems where this is "
    "possible: off=no CPU binding, auto=automatic binding(default). "
    "Values other than \"off\" and \"auto\" must be a hexadecimal string."
    "The lowest order bit is for the first logical CPU. For example, \"a5\" " 
    "and \"A5\" indicate that CPUs 0, 2, 5, and 7 are available for binding "
    "to threads, since hex value a5 =2^7+2^5+2^2+2^0.",
    storedOptions_.cpuMask_);

  AddStoredOption("alg:droptol droptol",
    "If droptol > 0 is specified, linear constraint and objective coefficients "
    "less than droptol in magnitude are treated as zero.",
    storedOptions_.dropTol_);

  AddStoredOption("tech:numcores numcores",
    "Write number of logical cores to stdout:\n"
    "\n.. value-table::\n",
    storedOptions_.numcores_, values_01_noyes_0default_);

  AddStoredOption("tech:outlev outlev",
    "Whether to write CPLEX log lines (chatter) to stdout, "
    "for granular control see \"tech:lpdisplay\", \"tech:mipdisplay\", \"tech:bardisplay\". "
    "Values:\n"
    "\n.. value-table::\n",
    storedOptions_.outlev_, outlev_values_);

  AddCplexParamOption("tech:bardisplay bardisplay",
    "Specifies how much the barrier algorithm chatters:\n"
    "\n.. value-table::\n",
    CPX_PARAM_BARDISPLAY, bardisplay_values_, 0);

  AddCplexParamOption("tech:lpdisplay display lpdisplay",
    "Frequency of displaying LP progress information:\n"
    "\n.. value-table::\n",
    CPX_PARAM_SIMDISPLAY, display_values_, 0);

  AddCplexParamOption("tech:mipdisplay mipdisplay",
    "Frequency of displaying branch-and-bound information:\n"
    "\n.. value-table::\n",
    CPX_PARAM_MIPDISPLAY, mipdisplay_values_, 0);

    AddCplexParamOption("tech:mipinterval mipinterval",
      "Frequency of node logging for \"tech::mipdisplay\" >=2:\n"
      "\n.. value-table::\n",
      CPX_PARAM_MIPINTERVAL, mipinterval_values_, 0);

  AddStoredOption("tech:logfile logfile",
    "Log file name.", storedOptions_.logFile_);

  AddCplexParamOption("tech:seed seed",
    "Seed for random number generator used internally "
    "by CPLEX. Use \"seed=?\" to see the default, which "
    "depends on the CPLEX release.",
    CPX_PARAM_RANDOMSEED, 0, CPXINT_MAX);

  AddCplexParamOption("tech:threads threads",
      "Default maximum number of threads for any of the parallel CPLEX "
      "optimizers (simplex, barrier, MIP, concurrent); default 0 ==> "
      "automatic choice (maximum available, subject to the license).",
      CPXPARAM_Threads, 0, CPXINT_MAX);

  AddStoredOption("tech:workdir nodefiledir workfiledir workdir",
    "Directory where CPLEX creates a temporary subdirectory "
    "for temporary files, e.g., for node information and Cholesky factors.",
    storedOptions_.workDir_
  );

  AddCplexParamOption("tech:workfilelim workfilelim",
    "Maximum size in megabytes for in-core work \"files\" "
    "(central memory used before node files are written); default 2048.",
    CPX_PARAM_WORKMEM, 0.0, 1e75);



  AddCplexParamOption("lim:iter iterlim iterlimit iterations lpiterlim",
    "LP iteration limit (default:  large).",
    CPXPARAM_Simplex_Limits_Iterations, 0, CPXINT_MAX);

  AddCplexParamOption("lim:netiterations netiterations netiter",
    "Limit on network simplex iterations (default: large).",
    CPXPARAM_Network_Iterations, 0, CPXINT_MAX);

  AddCplexParamOption("lim:nodes nodes node nodelim nodelimit",
    "Maximum MIP nodes to explore (default: no limit).",
    CPX_PARAM_NODELIM, 0, CPXINT_MAX);

  AddCplexParamOption("lim:dettime dettimelim",
    "Time limit in platform-dependent \"ticks\".",
    CPXPARAM_DetTimeLimit, 0.0, DBL_MAX);

  AddCplexParamOption("lim:sol sollimit solutionlimit mipsolutions solutionlim",
    "Limit the number of feasible MIP solutions found, causing early "
    "termination if exceeded; must be >= 1; default = no limit.",
    CPXPARAM_MIP_Limits_Solutions , 1, CPXINT_MAX);

  AddCplexParamOption("lim:time timelim timelimit time lptimelim",
                  "Limit on solve time (in seconds; default: no limit).",
                  CPXPARAM_TimeLimit, 0.0, DBL_MAX);

  AddCplexParamOption("lim:upperobj upperobjlim upperobj",
                  "Sets an upper limit on the value of the objective "
                  "function in the simplex algorithms. Setting an upper "
                  "objective function limit causes CPLEX to halt the "
                  "optimization process when the maximum objective "
                  "function value limit has been reached. This limit "
                  "applies only during Phase II of the simplex algorithm "
                  "in maximization problems. Default: 1e+75.",
                  CPXPARAM_Simplex_Limits_UpperObj, -DBL_MAX, DBL_MAX);

  AddCplexParamOption("lim:lowerobj lowerobjlim lowerobj",
                  "See lim:upperobjlim.",
                  CPXPARAM_Simplex_Limits_LowerObj, -DBL_MAX, DBL_MAX);

  AddCplexParamOption("lim:upperobjstop upperobjstop mip:upperobj",
                  "In a maximization MILP or MIQP, the solver will abort "
                  "the optimization process as soon it finds a solution "
                  "of value greater than or equal to the specified value. "
                  "It is ignored if the problem is multiobjective. "
                  "Default: 1e+75.",
                  CPXPARAM_MIP_Limits_UpperObjStop, -DBL_MAX, DBL_MAX);

  AddCplexParamOption("lim:lowerobjstop lowerobjstop mip:lowerobj",
                  "See lim:upperobjstop.",
                  CPXPARAM_MIP_Limits_LowerObjStop, -DBL_MAX, DBL_MAX);

  AddCplexParamOption("lim:uppercutoff uppercutoff mip:uppercutoff",
                  "Sets the upper cutoff tolerance. When the problem is "
                  "a minimization problem, CPLEX cuts off or discards "
                  "any solutions that are greater than the specified "
                  "upper cutoff value. If the model has no solution with "
                  "an objective value less than or equal to the cutoff "
                  "value, CPLEX declares the model infeasible. In other "
                  "words, setting an upper cutoff value c for a "
                  "minimization problem is similar to adding this "
                  "constraint to the objective function of the model: "
                  "obj <= c. "
                  "Default: 1e+75.\n"
                  "\n"
                  "Only effective in the branch and bound algorithm, "
                  "for example, in a mixed integer program (MIP).",
                  CPXPARAM_MIP_Tolerances_UpperCutoff, -DBL_MAX, DBL_MAX);

  AddCplexParamOption("lim:lowercutoff lowercutoff mip:lowercutoff",
                  "See lim:uppercutoff.",
                  CPXPARAM_MIP_Tolerances_LowerCutoff, -DBL_MAX, DBL_MAX);

  AddCplexParamOption("mip:objdiff objdiff objdif objdifference",
      "Used to update the cutoff each time a solution is found. "
      "This absolute value is subtracted from (added to) the newly "
      "found integer objective value when minimizing (maximizing). "
      "This forces the mixed integer optimization to ignore integer "
      "solutions that are not at least this amount better than the "
      "best one found so far. Positive values may speed "
      "the search -- and may cause the optimal solution "
      "to be missed.", CPXPARAM_MIP_Tolerances_ObjDifference,
      0.0, DBL_MAX);

  AddCplexParamOption("mip:relobjdiff relobjdiff relobjdif",
      "If objdiff parameter is 0, relobjdiff times the absolute "
      "value is used as a cutoff (see objdiff)", 
      CPXPARAM_MIP_Tolerances_RelObjDifference, 0.0, DBL_MAX);

  AddCplexParamOption("qp:target optimalitytarget",
    "Type of solution to compute for a (MI)QP (not QCP) problem:"
                  "\n\n.. value-table::\n",
    CPXPARAM_OptimalityTarget, optimalitytarget_values_, 0);

  
  AddCplexParamOption("pre:aggfill aggfill agglim",
    "Variables that appear in more than agglim rows (default 10) "
    "will not be substituted away by the aggregator.",
    CPXPARAM_Preprocessing_Fill, 0, CPXINT_MAX);

  AddCplexParamOption("pre:aggregate aggregate",
    "Whether to make substitutions to reduce the number of "
		"rows:\n\n.. value-table::\n",
    CPXPARAM_Preprocessing_Aggregator, values_preaggregator, -1);

  AddCplexParamOption("pre:coeffreduce coeffreduce",
    "Whether to use coefficient reduction when "
    "preprocessing MIPS\n\n.. value-table::\n",
    CPXPARAM_Preprocessing_CoeffReduce, values_precoeffreduce, -1);

  AddCplexParamOption("pre:dependency dependency",
    "Whether to use CPLEX's presolve dependency checker:\n\n.. value-table::\n",
    CPXPARAM_Preprocessing_Dependency, values_predependency, -1);

  AddCplexParamOption("pre:solve presolve", "Whether to use CPLEX's presolve:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Presolve, values_01_noyes_1default_, 1);

	AddStoredOption("pre:dual predual presolvedual",
		"Whether CPLEX' presolve phase should present the "
    "CPLEX solution algorithm with the primal(-1) or "
    "dual(1) problem or (default = 0) should decide "
		"automatically.\n"
									"\n"
									"For compatibility, pre:dual interacts with alg:dualproblem.",
		storedOptions_.predual_);

  AddStoredOption("alg:dualproblem dual",
									"Compatibility option with the legacy cplexasl driver. "
									"Interacts with pre:dual as follows:\n"
									"\n"
									"| alg:dualproblem - alone, equivalent to pre:dual=0\n"
									"| alg:dualproblem pre:dual=-1 - together, equivalent to pre:dual=1\n"
									"| alg:dualproblem pre:dual=0 - together, equivalent to pre:dual=0\n"
									"| alg:dualproblem pre:dual=1 - together, equivalent to pre:dual=-1.",
									storedOptions_.dualprob_);

  AddStoredOption("alg:primalproblem primal",
                  "Compatibility option with the legacy cplexasl driver. No effect.",
									storedOptions_.primalprob_dummy_);

  AddCplexParamOption("pre:node presolvenode",
    "Whether to run presolve at each node of the MIP branch-and-bound:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_PresolveNode, values_presolvenode, 0);

  AddCplexParamOption("pre:passes prepasses prepass",
    "Limit on the number of CPLEX presolve passes:\n"
    "\n"
    "| -1 - Automatic choice (default)\n"
    "| n>=0 - At most n passes.",
    CPXPARAM_Preprocessing_NumPass, -1, CPXINT_MAX);

  AddCplexParamOption("pre:reduce prereduce",
    "Kinds of reductions permitted during CPLEX presolve:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Reduce, values_prereduce, 3);

  AddCplexParamOption("pre:reformulations prereformulations",
    "Reformulations applied during CPLEX presolve:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Reformulations, values_prereformulations, 3);

  AddCplexParamOption("pre:relax prerelax relaxpresolve relpresolve",
    "Whether to use CPLEX's presolve on the initial LP relaxation of a MIP: \n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Relax, values_autonoyes_, 0);


  AddCplexParamOption("pre:scale scale",
    "How to scale the problem:\n"
    "\n.. value-table::\n",
    CPXPARAM_Read_Scale, values_prescale, 0);

  AddCplexParamOption("pre:sos1enc presos1enc presos1reform",
    "Encoding used for SOS1 reformulation:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_SOS1Reform, values_presosenc, 0);


  AddCplexParamOption("pre:sos2enc presos2enc presos2reform",
    "Encoding used for SOS2 reformulation, see pre:sos1enc.",
    CPXPARAM_Preprocessing_SOS2Reform, -1, 1);

  AddStoredOption("tech:nosolve nosolve",
    "Stop after loading the problem and honoring any "
    "\"writeprob\" or \"writemipstart\" directives.",
    storedOptions_.noSolve_);


    AddStoredOption("tech:pretunefileprm pretunefileprm pretunefile",
      "File to which nondefault keyword settings are written "
    "in CPLEX PRM format before tuning.",
      storedOptions_.pretunefileprm);


    AddCplexParamOption("tech:tunedisplay tunedisplay",
        "How much to print during tuning:\n"
        "\n.. value-table::\n",
        CPXPARAM_Tune_Display, values_tunedisplay_, 1);


  AddCplexParamOption("tech:tunerepeat tunerepeat",
      "How many times to perturb the problem during tuning (default = 1).",
      CPXPARAM_Tune_Repeat, 1, CPXINT_MAX);



  AddCplexParamOption("tech:tunetimelim tunetimelim tunetime",
      "Limit (in seconds) on tuning time; meaningful if < time (default = 1e75)",
      CPXPARAM_Tune_TimeLimit, 0.0, DBL_MAX);

  AddCplexParamOption("tech:tunetimedet tunetimedet",
      "Limit (in \"ticks\") on tuning time; meaningful if < time (default = 1e75)",
      CPXPARAM_Tune_DetTimeLimit, 0.0, DBL_MAX);


  AddStoredOption("tech:tunefix tunefix",
      "List of option names not to tune (their current values are kept), "
      "separated by commas or white space; enclose in quotes if the list "
      "contains white space.",
      storedOptions_.tunefix);

  AddStoredOption("tech:tunefixfile tunefixfile",
      "Name of a file containing option names not to tune, merged with "
      "tech:tunefix.",
      storedOptions_.tunefixfile);

  AddStoredOption("tech:tunefile tunefile tunefileprm tech:tunebase tunebase",
      "Name of file for tuning results in CPLEX PRM format. "
      "If specified, CPLEX will experiment with parameter "
      "settings  that would make solving faster.  This can "
      "significantly increase execution time of the current "
      "invocation, but the settings it finds might save time "
                  "in future runs.",
      storedOptions_.tunefileprm);


  AddToOptionDescription("tech:writemodel",
    "Cplex-specific file name extensions are \".sav\", \".mps\", "
    "\".lp\", \".rmp\",  \".rew\", \".rlp\"");

  AddStoredOption("tech:endbasis writebas endbasis",
    "Write the final basis to the specified file (in BAS format).",
    storedOptions_.endBasis_);

  AddStoredOption("tech:writemipstart writemipstart",
    "The name of a file to which the MIP starting guess(if any) "
    "is written in \".mst\" format.",
    storedOptions_.mipStart_);


 

  // --------- Parameters known from the cplexasl driver ---------
  AddCplexParamOption("lp:markowitz markowitz",
    "Markowitz pivot tolerance; default = 0.01; must be between "
    "1e-4 and 0.99999.  Bigger values may improve numerical properties "
    "of the solution (and may take more time).",
    CPXPARAM_Simplex_Tolerances_Markowitz, 1e-4, 0.99999);

  AddCplexParamOption("lp:pgradient pgradient",
    "Pricing algorithm for primal simplex (default 0):\n"
    "\n.. value-table::\n",
    CPXPARAM_Simplex_PGradient, values_pgradient, 0);

  AddCplexParamOption("lp:pricing pricing",
    "Size of the pricing candidate list for partial pricing; "
    "default 0 ==> automatic choice.",
    CPXPARAM_Simplex_Pricing, 0, CPXINT_MAX);

  AddCplexParamOption("lp:refactor refactor",
    "LP iterations between refactorizing the basis; "
    "default 0 ==> automatic choice; maximum 10000.",
    CPXPARAM_Simplex_Refactor, 0, 10000);

  AddCplexParamOption("lp:singular singular singularlim",
    "Maximum number of times CPLEX tries to repair the basis "
    "when it encounters singularities; default 10.",
    CPXPARAM_Simplex_Limits_Singularity, 0, CPXINT_MAX);

  AddCplexParamOption("alg:siftsimplex simplexsifting",
    "0/1*: Whether the simplex algorithm may use sifting when appropriate.",
    CPXPARAM_Sifting_Simplex, 0, 1);

  AddCplexParamOption("mip:heurfreq heurfreq heuristicfreq",
    "How often to apply the periodic (node) heuristic for MIPs:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_HeuristicFreq, values_heurfreq, 0);

  AddCplexParamOption("mip:heureffort heureffort",
    "Effort spent on MIP heuristics: 0 disables them, values < 1 "
    "reduce and values > 1 increase the effort; default 1. "
    "Do not combine with a non-default mip:heurfreq.",
    CPXPARAM_MIP_Strategy_HeuristicEffort, 0.0, 1e75);

  AddCplexParamOption("mip:rinsheur rinsheur",
    "Frequency of the relaxation induced neighborhood search (RINS) "
    "heuristic:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_RINSHeur, values_rinsheur, 0);

  AddCplexParamOption("mip:lbheur lbheur",
    "0*/1: Whether to use a local branching heuristic to try to improve "
    "new incumbents found during a MIP search.",
    CPXPARAM_MIP_Strategy_LBHeur, 0, 1);

  AddCplexParamOption("mip:probe probe",
    "How much variable probing to do before MIP branching (sometimes "
    "dramatically affects performance, for better or worse):\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_Probe, values_probe, 0);

  AddCplexParamOption("lim:probetime probetime",
    "Limit in seconds on time spent probing; default 1e75.",
    CPXPARAM_MIP_Limits_ProbeTime, 0.0, 1e75);

  AddCplexParamOption("lim:probetimedet probetimedet",
    "Limit in \"ticks\" on time spent probing; default 1e75.",
    CPXPARAM_MIP_Limits_ProbeDetTime, 0.0, 1e75);

  AddCplexParamOption("mip:symmetry symmetry",
    "Whether to break symmetry during preprocessing of MIP problems:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Symmetry, values_symmetry, -1);

  AddCplexParamOption("pre:repeat repeatpresolve",
    "Whether to repeat CPLEX's presolve after root processing:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_RepeatPresolve, values_repeatpresolve, -1);

  AddCplexParamOption("mip:repairtries repairtries",
    "How many times to try to repair an infeasible MIP start:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Limits_RepairTries, values_repairtries, 0);

  AddCplexParamOption("mip:submipnodelim submipnodelim",
    "Limit on nodes explored when CPLEX solves a subMIP (RINS, "
    "solution polishing, MIP start processing); default 500.",
    CPXPARAM_MIP_SubMIP_NodeLimit, 1, CPXINT_MAX);

  AddCplexParamOption("mip:strongcand strongcand",
    "Length of the candidate list for strong branching; default 10.",
    CPXPARAM_MIP_Limits_StrongCand, 1, CPXINT_MAX);

  AddCplexParamOption("mip:strongit strongit",
    "Number of simplex iterations on each variable in the candidate list "
    "during strong branching; default 0 ==> automatic choice.",
    CPXPARAM_MIP_Limits_StrongIt, 0, CPXINT_MAX);

  AddCplexParamOption("lim:treememory treememory treelimit treememlim",
    "Maximum megabytes of memory for the branch-and-bound tree "
    "(default 1e75 ==> no limit).",
    CPXPARAM_MIP_Limits_TreeMemory, 0.0, 1e75);

  AddCplexParamOption("mip:polishafter:absgap polishafter_absmipgap",
    "Start polishing integer solutions after the absolute MIP gap "
    "is at most this value; default 0.",
    CPXPARAM_MIP_PolishAfter_AbsMIPGap, 0.0, 1e20);

  AddCplexParamOption("mip:polishafter:gap polishafter_mipgap",
    "Start polishing integer solutions after the relative MIP gap "
    "is at most this value; default 0.",
    CPXPARAM_MIP_PolishAfter_MIPGap, 0.0, 1.0);

  AddCplexParamOption("mip:polishafter:nodes polishafter_nodes",
    "Start polishing integer solutions after processing this many nodes; "
    "default: no limit.",
    CPXPARAM_MIP_PolishAfter_Nodes, 0, CPXINT_MAX);

  AddCplexParamOption("mip:polishafter:sols polishafter_intsol",
    "Start polishing integer solutions after finding this many "
    "integer-feasible solutions; default: no limit.",
    CPXPARAM_MIP_PolishAfter_Solutions, 1, CPXINT_MAX);

  AddCplexParamOption("mip:polishafter:time polishafter_time",
    "Start polishing integer solutions after finding at least one "
    "integer-feasible solution and spending this many seconds; "
    "default 1e75.",
    CPXPARAM_MIP_PolishAfter_Time, 0.0, 1e75);

  AddCplexParamOption("mip:polishafter:dettime polishafter_timedet",
    "Start polishing integer solutions after finding at least one "
    "integer-feasible solution and spending this many \"ticks\"; "
    "default 1e75.",
    CPXPARAM_MIP_PolishAfter_DetTime, 0.0, 1e75);

  AddCplexParamOption("mip:kappastats mipkappa",
    "Whether to compute MIP kappa statistics, which summarize the "
    "condition numbers of the optimal bases seen while solving:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_KappaStats, values_kappastats, 0);

  AddCplexParamOption("mip:miqcpstrat miqcpstrat",
    "Strategy for solving quadratically-constrained MIPs (MIQCP):\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_MIQCPStrat, values_miqcpstrat, 0);

  AddCplexParamOption("mip:ordertype ordertype mipordertype",
    "How to generate default priorities for integer variables "
    "when no priority order is given:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_OrderType, values_ordertype, 0);

  AddCplexParamOption("bar:ordering ordering",
    "Ordering algorithm used by the barrier algorithm:\n"
    "\n.. value-table::\n",
    CPXPARAM_Barrier_Ordering, values_barorder, 0);

  AddCplexParamOption("bar:qcpcomptol qcpconvergetol qcpcomptol",
    "Convergence tolerance on relative complementarity for problems "
    "with quadratic constraints; default 1e-7.",
    CPXPARAM_Barrier_QCPConvergeTol, 1e-12, 1e75);

  AddCplexParamOption("pre:qtolin qtolin",
    "Whether to linearize products of bounded variables in quadratic "
    "objectives during preprocessing:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_QToLin, values_autonoyes_, -1);

  AddCplexParamOption("tech:clocktype clocktype",
    "Kind of times CPLEX reports and uses for time limits:\n"
    "\n.. value-table::\n",
    CPXPARAM_ClockType, values_clocktype, 2);

  AddCplexParamOption("tech:parallelmode parallelmode",
    "Parallel optimization mode:\n"
    "\n.. value-table::\n",
    CPXPARAM_Parallel, values_parallelmode, 0);

  AddCplexParamOption("tech:datacheck datacheck",
    "Data consistency checking and modeling assistance:\n"
    "\n.. value-table::\n",
    CPXPARAM_Read_DataCheck, values_datacheck, 0);

  AddCplexParamOption("tech:warninglimit warninglimit",
    "Limit on the number of warnings per issue given when "
    "tech:datacheck=2; default 10.",
    CPXPARAM_Read_WarningLimit, 0, CPXINT_MAX);

  AddCplexParamOption("tech:record record",
    "0*/1: Whether to record CPLEX library calls for debugging use by IBM "
    "in a file with an automatically chosen name of the form cplexXXXXXXX.db.",
    CPXPARAM_Record, 0, 1);

  // --------- Parameters in neither driver before ---------
  AddCplexParamOption("pre:folding folding",
    "Whether to try LP folding during preprocessing:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Folding, values_folding, -1);

  AddCplexParamOption("mip:dive dive",
    "MIP dive strategy (how to choose between child nodes):\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_Dive, values_dive, 0);

  AddCplexParamOption("mip:cardls cardls",
    "Cardinality local search heuristic:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_CardLs, values_cardls, 0);

  AddCplexParamOption("lp:dynamicrows dynamicrows",
    "Whether CPLEX manages the rows of the basis matrix dynamically "
    "during dual simplex:\n"
    "\n.. value-table::\n",
    CPXPARAM_Simplex_DynamicRows, values_autonoyes_, -1);

  AddCplexParamOption("pre:qcpduals qcpduals",
    "Whether to compute dual values for quadratically constrained "
    "problems:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_QCPDuals, values_qcpduals, 1);

  AddCplexParamOption("pre:qpmakepsd qpmakepsd",
    "0/1*: Whether to reformulate a binary-only MIQP or MIQCP so that "
    "its quadratic matrix becomes positive semi-definite.",
    CPXPARAM_Preprocessing_QPMakePSD, 0, 1);

  AddCplexParamOption("mip:eplin eplin",
    "Epsilon used in CPLEX's linearization of logical constraints, "
    "e.g., x < a becomes x <= a - eplin; default 1e-3.",
    CPXPARAM_MIP_Tolerances_Linearization, 0.0, 1e75);

  AddCplexParamOption("alg:feasrelax:tol feasopttol",
    "Tolerance for constraint violations in the relaxed problem of "
    "alg:feasrelax; default 1e-6.",
    CPXPARAM_Feasopt_Tolerance, 0.0, 1e20);

  AddCplexParamOption("tech:writelevel writelevel",
    "Level of detail written to SOL and MST files "
    "(tech:writesolution, tech:writemipstart):\n"
    "\n.. value-table::\n",
    CPXPARAM_Output_WriteLevel, values_writelevel, 0);

  AddCplexParamOption("tech:intsolfileprefix intsolfileprefix",
    "If set, CPLEX writes every MIP integer solution it finds to a file "
    "whose name starts with this prefix.",
    CPXPARAM_Output_IntSolFilePrefix);

  AddCplexParamOption("tech:mpslong mpslong",
    "0/1*: Whether to write numbers with full precision to MPS and REW files.",
    CPXPARAM_Output_MPSLong, 0, 1);

  AddCplexParamOption("tech:clonelog clonelog",
    "Whether clone log files are created during parallel optimization:\n"
    "\n.. value-table::\n",
    CPXPARAM_Output_CloneLog, values_clonelog, 0);

  AddCplexParamOption("tech:tunemeasure tunemeasure",
    "How tuning compares parameter sets across models:\n"
    "\n.. value-table::\n",
    CPXPARAM_Tune_Measure, values_tunemeasure, 1);

  AddCplexParamOption("alg:siftalg siftalg",
    "Algorithm for sifting subproblems:\n"
    "\n.. value-table::\n",
    CPXPARAM_Sifting_Algorithm, values_siftalg, 0);

  AddCplexParamOption("lim:siftiter siftiterations siftiter",
    "Limit on sifting iterations (default: no limit).",
    CPXPARAM_Sifting_Iterations, 0, CPXINT_MAX);

  AddCplexParamOption("tech:siftdisplay siftdisplay",
    "How much the sifting algorithm reports:\n"
    "\n.. value-table::\n",
    CPXPARAM_Sifting_Display, values_siftdisplay, 1);

  /// Custom solve results here'...
  //  AddSolveResults({
  //                    { sol::NUMERIC, "failure: numeric issue, no feasible solution" }
  //                  });
}

void CplexBackend::CplexSetObjIntParam(const SolverOption& opt, int val) {
  objnparam_int_.push_back({ {opt.wc_tail(), opt.wc_keybody_last()}, val });
}
void CplexBackend::CplexSetObjDblParam(const SolverOption& opt, double val) {
  objnparam_dbl_.push_back({ {opt.wc_tail(), opt.wc_keybody_last()}, val });
}
int CplexBackend::CplexGetObjIntParam(const SolverOption& opt) const {
  auto it = std::find_if(objnparam_int_.rbegin(), objnparam_int_.rend(),
    [&](const ObjNParam<int>& prm) {
      return prm.first == std::make_pair(opt.wc_tail(), opt.wc_keybody_last());
    });
  if (objnparam_int_.rend() == it) {
    if (debug_mode())
      return INT_MIN;
    throw std::runtime_error("Failed to find recorded option " +
      opt.wc_key_last__std_form());
  }
  return it->second;
}
double CplexBackend::CplexGetObjDblParam(const SolverOption& opt) const {
  auto it = std::find_if(objnparam_dbl_.rbegin(), objnparam_dbl_.rend(),
    [&](const ObjNParam<int>& prm) {
      return prm.first == std::make_pair(opt.wc_tail(), opt.wc_keybody_last());
    });
  if (objnparam_dbl_.rend() == it) {
    if (debug_mode())
      return NAN;
    throw std::runtime_error("Failed to find recorded option " +
      opt.wc_key_last__std_form());
  }
  return it->second;
}


void CplexBackend::ReadBendersSuffix()
{
  int strategy = 0;
  CPXgetintparam(env(), CPXPARAM_Benders_Strategy, &strategy);
  auto suf_mask = suf::Kind::VAR_BIT;
  if (auto mv0 = ReadModelSuffixInt({ "benders", suf_mask }))
  {
    auto mv = GetValuePresolver().PresolveGenericInt(mv0);
    auto values = mv.GetVarValues()();
    int nmaster = 0, nsub = 0;
    // Count the vars in master and in non-master
    for (auto v : values)
      if (v) ++nsub; else ++nmaster;
    if (!nmaster || !nsub) {
      AddToSolverMessage(fmt::format("Ignoring .benders because {} variables have .benders = 0.\n",
        nmaster ? "all" : "no"));
      storedOptions_.cpxMethod_ = CPX_ALG_MIP;
      return;
    }

    // Annotate
    CPLEX_CALL(CPXnewlongannotation(env(), lp(), CPX_BENDERS_ANNOTATION, 0));
    int a;
    CPLEX_CALL(CPXgetlongannotationindex(env(), lp(), CPX_BENDERS_ANNOTATION, &a));
    
    // Get subproblem indices
    std::vector<int> indices(nsub);
    std::vector<CPXLONG> subproblem(nsub);
    for (int i = 0, j = 0; i < values.size(); ++i)
    {
      if (values[i] > 0)
      {
        indices[j] = i;
        subproblem[j++] = values[i];
      }
    }
    CPLEX_CALL(CPXsetlongannotations(env(), lp(), a, CPX_ANNOTATIONOBJ_COL,
      nsub, indices.data(), subproblem.data()));
  }
  else
  {
    if (strategy)
    {
      AddToSolverMessage(fmt::format("Ignoring bendersopt because suffix benders not present\n"
        "and benders_strategy = {}, which requires .benders.\n", strategy));
      storedOptions_.cpxMethod_ = CPX_ALG_MIP;
    }
  }
    
   

}
/// What to do on certain "obj:*:..." option
static std::tuple<int, CplexBackend::CplexObjParams>
CplexGetObjParamAction(const CplexBackend::ObjNParamKey& key) {
  int n;
  try {
    n = std::stoi(key.second) - 1;  // subtract 1 for 0-based indexing
  }
  catch (...) {
    throw std::runtime_error("Could not parse index '" + key.second +
      "' of option 'obj:" + key.second + key.first + "'");
  }
  if (":priority" == key.first)
    return { n, CplexBackend::OBJ_PRIORITY};
  if (":weight" == key.first)
    return { n, CplexBackend::OBJ_WEIGHT };
  if (":abstol" == key.first)
    return { n, CplexBackend::OBJ_ABSTOL };
  if (":reltol" == key.first)
    return { n, CplexBackend::OBJ_RELTOL};
  throw std::runtime_error(
    "Unknown wildcard option 'obj:" + key.second + key.first + "'");
  return { -1,  CplexBackend::OBJ_NOTVALID };
}

/// env() is only for error reporting
static void CplexDoSetObjParam(
  const CplexBackend::ObjNParam<int>& prm,
  CPXLPptr model, CPXENVptr env) {
  auto action = CplexGetObjParamAction(prm.first);
  auto iobj = std::get<0>(action);
  auto prm_attr = std::get<1>(action);
  if (prm_attr != CplexBackend::OBJ_PRIORITY)
    return;

  int status = CPXmultiobjchgattribs(env, model, iobj,
    CPX_NO_OFFSET_CHANGE, CPX_NO_WEIGHT_CHANGE, prm.second,
    CPX_NO_ABSTOL_CHANGE, CPX_NO_RELTOL_CHANGE, NULL);
  if (status)
    throw CplexCommon::GetException("CPXmultiobjchgattribs", status, env);
}

static void CplexDoSetObjParam(
  const CplexBackend::ObjNParam<double>& prm,
  CPXLPptr model, CPXENVptr env) {
  auto action = CplexGetObjParamAction(prm.first);
  auto iobj = std::get<0>(action);
  auto prm_attr = std::get<1>(action);
  double weight = CPX_NO_WEIGHT_CHANGE;
  double abstol = CPX_NO_ABSTOL_CHANGE;
  double reltol = CPX_NO_RELTOL_CHANGE;
  
  switch (prm_attr) {
  case CplexBackend::OBJ_WEIGHT:
    weight = prm.second;
    break;
  case CplexBackend::OBJ_ABSTOL:
    abstol = prm.second;
    break;
  case CplexBackend::OBJ_RELTOL:
    reltol = prm.second;
    break;
  }
    int status = CPXmultiobjchgattribs(env, model, iobj,
      CPX_NO_OFFSET_CHANGE, weight,CPX_NO_PRIORITY_CHANGE,
      abstol, reltol, NULL);
    if (status)
      throw CplexCommon::GetException("CPXmultiobjchgattribs", status, env);
}

template <class T>
static void DoPlayCplexObjNParams(
  const std::vector< CplexBackend::ObjNParam<T> >& objnp,
  CPXLPptr model, CPXENVptr env) {
  for (const auto& p : objnp)
    CplexDoSetObjParam(p, model, env);
}

void CplexBackend::CplexPlayObjNParams() {
  DoPlayCplexObjNParams(objnparam_int_, lp(), env());
  DoPlayCplexObjNParams(objnparam_dbl_, lp(), env());
}
void CplexBackend::RecordCplexParamNames(const char* name_list, int key) {
  std::string names(name_list);
  size_t pos = 0;
  while (pos < names.size()) {
    auto end = names.find(' ', pos);
    if (std::string::npos == end)
      end = names.size();
    if (end > pos)
      cplexParamOfOption_[names.substr(pos, end - pos)] = key;
    pos = end + 1;
  }
}

/// Resolve the option names in tech:tunefix / tech:tunefixfile
/// to CPLEX parameters with their current values.
void CplexBackend::CollectTuneFixed(const std::string& list, const char* what,
                                    std::vector<int>& intnum, std::vector<int>& intval,
                                    std::vector<int>& dblnum, std::vector<double>& dblval) {
  size_t pos = 0;
  while (pos < list.size()) {
    while (pos < list.size() && (std::isspace((unsigned char)list[pos]) || list[pos] == ','))
      ++pos;
    auto end = pos;
    while (end < list.size() && !std::isspace((unsigned char)list[end]) && list[end] != ',')
      ++end;
    if (end == pos)
      break;
    auto name = list.substr(pos, end - pos);
    pos = end;
    auto it = cplexParamOfOption_.find(name);
    if (cplexParamOfOption_.end() == it)
      MP_RAISE(fmt::format("Unrecognized option \"{}\" in {}", name, what));
    int key = it->second, type = CPX_PARAMTYPE_NONE;
    CPXgetparamtype(env(), key, &type);
    switch (type) {
    case CPX_PARAMTYPE_INT:
    case CPX_PARAMTYPE_LONG: {
      if (std::find(intnum.begin(), intnum.end(), key) == intnum.end()) {
        int v = 0;
        CPXgetintparam(env(), key, &v);
        intnum.push_back(key);
        intval.push_back(v);
      }
      break;
    }
    case CPX_PARAMTYPE_DOUBLE: {
      if (std::find(dblnum.begin(), dblnum.end(), key) == dblnum.end()) {
        double v = 0.0;
        CPXgetdblparam(env(), key, &v);
        dblnum.push_back(key);
        dblval.push_back(v);
      }
      break;
    }
    default:
      MP_RAISE(fmt::format("Option \"{}\" in {} is not a tunable "
                           "numeric parameter", name, what));
    }
  }
}

void CplexBackend::DoCPLEXTune() {
  int status = 0;

  if (!storedOptions_.pretunefileprm.empty())
    CPXwriteparam(env(), storedOptions_.pretunefileprm.c_str());

  std::vector<int> intnum, intval, dblnum;
  std::vector<double> dblval;
  if (!storedOptions_.tunefix.empty())
    CollectTuneFixed(storedOptions_.tunefix, "tech:tunefix",
                     intnum, intval, dblnum, dblval);
  if (!storedOptions_.tunefixfile.empty()) {
    std::ifstream f(storedOptions_.tunefixfile);
    if (!f)
      MP_RAISE(fmt::format("Could not open tech:tunefixfile \"{}\"",
                           storedOptions_.tunefixfile));
    std::string contents((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());
    CollectTuneFixed(contents, "tech:tunefixfile",
                     intnum, intval, dblnum, dblval);
  }

  int res = CPXtuneparam(env(), lp(),
                         (int)intnum.size(), intnum.data(), intval.data(),
                         (int)dblnum.size(), dblnum.data(), dblval.data(),
                         0, nullptr, nullptr, &status);
  if (res)
    AddWarning("TuneFailed",
               fmt::format("CPXtuneparam returned error {}", res));
  else if (status)
    AddToSolverMessage(fmt::format(
        "Tuning stopped early (tuning status {}).\n", status));

  if (!storedOptions_.tunefileprm.empty())
    CPXwriteparam(env(), storedOptions_.tunefileprm.c_str());
}

} // namespace mp


AMPLS_MP_Solver* Open_cplex(CCallbacks cb = {}) {
  AMPLS_MP_Solver* slv = 
    AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::CplexBackend()},
    cb);
  return slv;
}

void AMPLSClose_cplex(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_cplex(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CplexBackend*>(AMPLSGetBackend(slv))->lp();
}

void* AMPLSGetEnv_cplex(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CplexBackend*>(AMPLSGetBackend(slv))->env();
}

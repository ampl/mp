#include <vector>
#include <climits>
#include <cfloat>
#include <cassert>
#include <fstream> 

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
      hasSolution_ = (status == sol::LIMIT_FEAS) || (status == sol::SOLVED) || (status == sol::UNCERTAIN);
    }
    return hasSolution_;
  }

  double  CplexBackend::MIPGap() {
    if (HasSolution()) {
      getAndReturnDblParam(CPXgetmiprelgap);
    }
    else return AMPLInf();
  }

  double  CplexBackend::MIPGapAbs() {
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
    bool outputStatus = silenceOutput_;
    setSilenceOutput(true);
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
    auto varstt = mv.GetVarValues()();
    auto constt = mv.GetConValues()(CG_Linear);
    assert(varstt.size());
    assert(constt.size());
    VarConStatii(varstt, constt);
  }

  void CplexBackend::AddPrimalDualStart(Solution sol0_unpres) {
    auto mv = GetValuePresolver().PresolveSolution(
      { sol0_unpres.primal, sol0_unpres.dual });
    auto x0 = mv.GetVarValues()();
    auto pi0 = mv.GetConValues()(CG_Linear);
    CPXcopystart(env(), lp(), nullptr, nullptr, x0.data(), nullptr,
      nullptr, pi0.data());
  }

  void CplexBackend::AddMIPStart(
    ArrayRef<double> x0_unpres, ArrayRef<int> sparsity_unpres) {
    if (!IsMIP()) return;
    auto mv = GetValuePresolver().PresolveSolution({ x0_unpres });
    auto ms = GetValuePresolver().PresolveGenericInt({ sparsity_unpres });
    auto x0 = mv.GetVarValues()();
    auto s0 = ms.GetVarValues()();
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
  // if (HasSolution() && ((!IsMIP()) || need_fixed_MIP()) && (NumLinCons()>0))
  {
    int num_cons = NumLinCons();
    std::vector<double> pi(num_cons);
    int error = CPXgetpi(env(), lp(), pi.data(), 0, num_cons - 1);
    if (error)
      pi.clear();
    return pi;
  }
  // else
  //   return std::vector<double>();
}

double CplexBackend::ObjectiveValue() const {
  double objval = -Infinity();
  if(hasSolution_)
    CPLEX_CALL(CPXgetobjval(env(), lp(), &objval));
  return objval;
}

double CplexBackend::NodeCount() const {
  return  IsMIP() ? CPXgetnodecnt(env(), lp()) : 0;
}

double CplexBackend::SimplexIterations() const {
  return IsMIP() ? CPXgetmipitcnt(env(), lp()) : CPXgetitcnt(env(), lp());
}

int CplexBackend::BarrierIterations() const {
  if (!IsMIP() && (storedOptions_.cpxMethod_ == CPX_ALG_BARRIER))
    return CPXgetbaritcnt(env(), lp());
  else
    return 0;
}

void CplexBackend::DoWriteProblem(const std::string &file) {
  CPLEX_CALL( CPXwriteprob (env(), lp(), file.c_str(), NULL) );
}
void CplexBackend::DoWriteSolution(const std::string& file) {
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
    if (NumObjs() > 1)
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
    setSilenceOutput(!verbose_mode());
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

void CplexBackend::InputExtras() {
  BaseBackend::InputExtras();
  InputCPLEXExtras();
}

void CplexBackend::InputCPLEXExtras() {
  // Set output on screen
  int lp, mip, bar, mo, netw;
  GetSolverOption(CPX_PARAM_SIMDISPLAY, lp);
  GetSolverOption(CPX_PARAM_MIPDISPLAY, mip);
  GetSolverOption(CPX_PARAM_BARDISPLAY, bar);
  GetSolverOption(CPXPARAM_MultiObjective_Display, mo);
  GetSolverOption(CPX_PARAM_NETDISPLAY, netw);
  if (storedOptions_.outlev_ > 2)
    storedOptions_.outlev_ = 2;
  int olp[] = { 0, 1, 2 };
  int omip[] = { 0, 3, 5 };
  lp = lp ? lp : olp[storedOptions_.outlev_];
  mip = mip ? mip : omip[storedOptions_.outlev_];
  bar = bar ? bar : olp[storedOptions_.outlev_];
  mo = mo ? mo : olp[storedOptions_.outlev_];
  netw = netw ? netw : olp[storedOptions_.outlev_];
  if (lp || mip || bar || mo || netw) {
    /* Log messages on screen */
    CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ScreenOutput, CPX_ON));
    /* Echo changed params before solve */
    CPLEX_CALL(CPXsetintparam(env(), CPXPARAM_ParamDisplay, 1));
  }
  SetSolverOption(CPX_PARAM_SIMDISPLAY, lp);
  SetSolverOption(CPX_PARAM_MIPDISPLAY, mip);
  SetSolverOption(CPX_PARAM_BARDISPLAY, bar);
  SetSolverOption(CPXPARAM_MultiObjective_Display, mo);
  SetSolverOption(CPX_PARAM_NETDISPLAY, netw);
  if (!storedOptions_.logFile_.empty())
  {
    if (lp < 1) SetSolverOption(CPX_PARAM_SIMDISPLAY, 1);
    if (mip < 1) SetSolverOption(CPX_PARAM_MIPDISPLAY, 1);
    CPLEX_CALL(CPXsetlogfilename(env(), storedOptions_.logFile_.data(), "w"));
  }
  set_verbose_mode(storedOptions_.outlev_ > 0);

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
  auto msg = DoCplexFixedModel();
  if (!msg.empty()) {
    AddToSolverMessage(msg +
      " failed in DoCplexFixedModel().");
    CPXfreeprob(env(),  &lp_ref());
    CPXchgprobtype(env(), lp(), original_model_type_);
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
  CPLEX_CALL(Contopt(env(), lp()));

  int optimstatus = CPXgetstat(env(), lp());
  if (optimstatus != CPX_STAT_OPTIMAL){
  }
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
  { "1", "equivalent to \"bardisplay\"=1, \"display\"=1, \"mipdisplay\"=3", 1},
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
  { "-1", "Automatic (default)", -1},
  { "0", "Primal simplex", 0},
  { "1", "Dual simplex", 1},
  { "2", "Barrier", 2},
  { "3", "Nondeterministic concurrent (several solves in parallel)", 3},
  { "4", "Network simplex", 4},
  { "5", "Sifting", 5}
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

static const mp::OptionValueInfo values_submipstartalg[] = {
  { "0", "Automatic (default)",0},
  { "1", "Primal simplex", 1},
  { "2", "Dual simplex", 2},
  { "3", "Network simplex", 3},
  { "4", "Barrier", 4},
  { "5", "Sifting", 5},
  { "6", "Concurrent (dual, barrier and primal in opportunistic mode; "
         "dual and barrier in deterministic mode; 4 is used for MIPQs).", 6}
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
static const mp::OptionValueInfo values_bardisplay[] = {
  { "0", "No output", 0},
  { "1", "Normal setup and iteration information (default)", 1},
  { "2", "Diagnostic information", 2}
};
static const mp::OptionValueInfo values_barstart[] = {
  { "1", "Assume dual is 0 (default)", 1},
  { "2", "Estimate dual", 2},
  { "3", "Average of primal estimate, 0 dual", 3},
  { "4", "Average of primal and dual estimates", 4}
};


static const mp::OptionValueInfo values_prereformulations[] = {
  { "0", "None", 0},
  { "1", "Allow reformulations that interphere with crushing forms", 1},
  { "2", "Allow reformulations that interphere with uncrushing forms", 2},
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
  { "0", "no", 0},
  { "1", "compressed node file in memory (default)", 1},
  { "2", "node file on disk", 2},
  { "3", "compressed node file on disk", 3}
};

  
void CplexBackend::setSolutionMethod() {
  int nFlags = bool(storedOptions_.fBarrier_)
    + bool(storedOptions_.fPrimal_)
    + bool(storedOptions_.fDual_)
    + bool(storedOptions_.fBenders_)
    + bool(storedOptions_.fNetwork_)
    + bool(storedOptions_.fSifting_);
  if (nFlags>= 2) 
    AddWarning("Ambiguous LP method",
      "Only one of barrier/primal/dual/network/sifting/benders should be specified.");
  if (nFlags >= 1)
  {
    if (storedOptions_.fPrimal_)
      storedOptions_.cpxMethod_ = CPX_ALG_PRIMAL;
    if (storedOptions_.fDual_)
      storedOptions_.cpxMethod_ = CPX_ALG_DUAL;
    if (storedOptions_.fBarrier_)
      storedOptions_.cpxMethod_ = CPX_ALG_BARRIER;
    if (storedOptions_.fNetwork_)
      storedOptions_.cpxMethod_ = CPX_ALG_NET;
    if (storedOptions_.fSifting_)
      storedOptions_.cpxMethod_ = CPX_ALG_SIFTING;
  }
  else {
    int mapMethods[] = {
      CPX_ALG_AUTOMATIC,
      CPX_ALG_PRIMAL,
      CPX_ALG_DUAL,
      CPX_ALG_BARRIER,
      CPX_ALG_CONCURRENT,
      CPX_ALG_NET,
      CPX_ALG_SIFTING
    };
    storedOptions_.cpxMethod_ = mapMethods[storedOptions_.algMethod_ + 1];
  }
  if (IsMIP())
    SetSolverOption(CPX_PARAM_STARTALG, storedOptions_.cpxMethod_);
  else if (IsQP())
    SetSolverOption(CPX_PARAM_QPMETHOD, storedOptions_.cpxMethod_);
  else
    SetSolverOption(CPX_PARAM_LPMETHOD, storedOptions_.cpxMethod_);

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

  if (storedOptions_.fBenders_)
  {
    storedOptions_.cpxMethod_ = CPX_ALG_BENDERS;
    ReadBendersSuffix();
  }
}

////////////////////////////// OPTIONS /////////////////////////////////

void CplexBackend::FinishOptionParsing() {
  // Apply to all cuts when not overriden
  if (storedOptions_.cuts_ != 2) {
    static int op[] = {
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
    int f;
    int local;
    for (f = 0; f < sizeof(op) / sizeof(int); f++)
    {
      CPXgetintparam(env(), op[f], &local);
      if(local!=0)
        CPXsetintparam(env(), op[f], storedOptions_.cuts_);
    }
  }

  if (!storedOptions_.mipStart_.empty())
  {
    // Write MIP starts
    auto nms = CPXgetnummipstarts(env(), lp());
    if (nms == 0)
      AddToSolverMessage("No MIP start is available.\n");
    CPLEX_CALL(CPXwritemipstarts(env(), lp(), storedOptions_.mipStart_.c_str(), 0, nms - 1));
  }

  if (!storedOptions_.cpuMask_.empty())
    CPLEX_CALL(CPXsetstrparam(env(), CPXPARAM_CPUmask, storedOptions_.cpuMask_.c_str()));

  int numcores;
  if (storedOptions_.numcores_) {
    CPLEX_CALL(CPXgetnumcores(env(), &storedOptions_.numcores_));
    AddToSolverMessage(fmt::format("{} logical cores are available.\n", storedOptions_.numcores_));
  }

  if (!storedOptions_.workDir_.empty())
    CPLEX_CALL(CPXsetstrparam(env(), CPX_PARAM_WORKDIR, storedOptions_.workDir_.c_str()));


}
void CplexBackend::InitCustomOptions() {

  set_option_header(
      "IBM ILOG CPLEX Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cplex_options``. For example::\n"
      "\n"
      "  ampl: option cplex_options 'mipgap=1e-6';\n");
  

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
  AddStoredOption("alg:method method lpmethod simplex mipstartalg",
    "Which algorithm to use for non-MIP problems or for the root node of MIP problems, unless "
    "primal/dual/barrier/network/sifting flags are specified:\n"
    "\n.. value-table::\n"
    "For MIQP problems (quadratic objective, linear constraints), setting 5 "
    "is treated as 0 and 6 as 4. For MIQCP problems(quadratic objective & "
    "constraints), all settings are treated as 4.", 
    storedOptions_.algMethod_, values_method);

  AddStoredOption("alg:barrier barrier baropt",
    "Solve (MIP root) LPs by barrier method.",
    storedOptions_.fBarrier_);

  AddStoredOption("alg:primal primalopt",
    "Solve (MIP root) LPs by primal simplex method.",
    storedOptions_.fPrimal_);

  AddStoredOption("alg:dual dual dualopt",
    "Solve (MIP root) LPs by dual simplex method.",
    storedOptions_.fDual_);

  AddStoredOption("alg:sifting sifting siftopt siftingopt",
    "Solve (MIP root) LPs by sifting method.",
    storedOptions_.fSifting_);

  AddStoredOption("alg:network network netopt",
    "Solve (substructure of) (MIP node) LPs "
    "by network simplex method. For best recognition of "
                  "network (sub)structures, switch off CPLEX presolve.",
    storedOptions_.fNetwork_);


  AddStoredOption("alg:benders benders bendersopt",
    "Solve MIP using Benders algorithm. Both integer and continuous "
    "variables must be present.",
    storedOptions_.fBenders_);
  

  AddSolverOption("alg:benders_worker benders_worker",
    "Designate the algorithm that CPLEX applies to solve the "
    "subproblems when using Benders decomposition:\n"
    "\n.. value-table::\n",
    CPXPARAM_Benders_WorkerAlgorithm,
    values_nodemethod, 0);

  AddSolverOption("alg:benders_feascuttol benders_feascut_tol",
    "Tolerance for violations of feasibility cuts in Benders "
    "algorithm(default 1e-6).",
    CPXPARAM_Benders_Tolerances_feasibilitycut,
    1e-9, 1.0);

  AddSolverOption("alg:benders_optcut_tol benders_optcut_tol",
    "Tolerance for violations of optimality cuts in Benders "
    "algorithm(default 1e-6).",
    CPXPARAM_Benders_Tolerances_optimalitycut,
    1e-9, 1.0);

  AddSolverOption("alg:benders_strategy benders_strategy",
    "How to decompose the problem for Benders algorithm:\n"
    "\n.. value-table::\n",
    CPXPARAM_Benders_Strategy,
    values_benders_strategy, 0);

  AddSolverOption("alg:feastol feastol",
    "Primal feasibility tolerance (default 1e-6, possible "
    "values are between 1e-9 and 0.1).",
    CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9, 1e-1);




  AddSolverOption("lp:crash crash",
    "Crash strategy (used to obtain starting basis in simplex); "
    "possible values = -1, 0, 1; default = 1. "
    "The best setting is problem-dependent and can only be found "
    "by experimentation. 0 completely ignores the objective.",
    CPXPARAM_Simplex_Crash, -1, 1);

  AddSolverOption("lp:dgradient dgradient",
    "Pricing algorithm for dual simplex (default 0):\n"
    "\n.. value-table::\n",
    CPXPARAM_Simplex_DGradient, 
    values_dgradient, 0);

  AddSolverOption("lp:doperturb doperturb",
    "Decides whether to perturb problems when using simplex, setting it "
    "to 1 is occasionally helpful for highly degenerate problems.:\n"
    "\n.. value-table::\n",
    CPXPARAM_Simplex_Perturbation_Indicator,
    values_01_noyes_0default_, 0);

  AddSolverOption("lp:perturblim perturblim perturblimit",
    "Number of stalled simplex iterations before the "
    "problem is perturbed; default=0 => automatic.",
    CPXPARAM_Simplex_Limits_Perturbation, 0, CPXINT_MAX);

  AddSolverOption("lp:perturb perturb",
    "Amount by which to perturb variable bounds when "
    "perturbing problems (see \"perturb\"); default=1e-6, "
    "must be positive.",
    CPXPARAM_Simplex_Perturbation_Constant, 1e-8, DBL_MAX);

  // Cut generation
  
  // Hidden param
  #define CPX_PARAM_EPSAGG 1055
  AddSolverOption("cut:aggtol aggtol",
    "Pivot tolerance for aggregating.  It seldom needs "
		"fiddling.  Default = .05; must be in [1e-10, .99].",
    CPX_PARAM_EPSAGG, 1e-10, .99);

  AddSolverOption("cut:aggforcut aggforcut",
    "Bound on the number of constraints aggregated to generate flow-cover "
    "and mixed-integer-rounding cuts (default 3).",
    CPXPARAM_MIP_Limits_AggForCut, 0, CPXINT_MAX);


  AddSolverOption("cut:bqp bqpcuts",
    "Whether to enable Boolean Quadric Polytope cut generation for nonconvex QP and "
    "MIQP problems, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_BQP, values_bqpcuts, 0);

  AddSolverOption("cut:clique cliquecuts",
    "Whether to use clique cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Cliques, values_bqpcuts, 0);

  AddSolverOption("cut:cover covercuts covers",
    "Whether to use cover cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Covers, values_bqpcuts, 0);

  AddSolverOption("cut:disj disjcuts",
    "Whether to use disjunctive cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Disjunctive, values_bqpcuts, 0);

  AddSolverOption("cut:eachcutlim eachcutlim",
    "Limit on the number of cuts of each time; default=2100000000.", 
    CPXPARAM_MIP_Limits_EachCutLimit, 0, CPXINT_MAX);

  AddSolverOption("cut:flowcover flowcovercuts flowcuts",
    "Whether to use flowcover cuts in solving MIPs:\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_Cuts_FlowCovers, values_cuts2, 0);


  AddSolverOption("cut:flowpath flowpathcuts",
    "Whether to use flow path cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_PathCut, values_cuts2, 0);

  AddSolverOption("cut:frac gomory fraccuts",
    "Whether to use fractional Gomory cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_Gomory, values_cuts2, 0);

  AddSolverOption("cut:fracpass fracpass gomorypass fractionalcuts",
    "Limit on number of passes to generate MIP fractional Gomory cuts, default=0 => automatic choice, "
    "positive = at most that many passes.",
    CPXPARAM_MIP_Limits_GomoryPass, 0, CPXINT_MAX);

  AddSolverOption("cut:fraccand fraccand gomorycand",
    "Limit on number of candidate variables when generating Gomory cuts "
    "for MIP problems; default=200.",
    CPXPARAM_MIP_Limits_GomoryCand, 1, CPXINT_MAX);


  AddSolverOption("cut:gubcover gubcover gubcuts",
    "Whether to use fractional GUB cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_GUBCovers, values_cuts2, 0);

  AddSolverOption("cut:implied implied impliedcuts",
    "Whether to use implied cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_Implied, values_cuts2, 0);

  AddSolverOption("cut:liftandproject liftandproject liftandprojectcuts",
    "Whether to use lift and project cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_LiftProj, values_bqpcuts, 0);

  AddSolverOption("cut:localimplied localimplied localimpliedcuts",
    "Whether to use locally valid implied cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_LocalImplied, values_bqpcuts, 0);

  AddSolverOption("cut:mfc mfccuts",
    "Whether to use multi-commodity flow (MCF) cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_MCFCut, values_cuts2, 0);

  AddSolverOption("cut:mir mircuts",
    "Whether to use MIP roundind cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_MIRCut, values_cuts2, 0);

  AddSolverOption("cut:node nodecuts",
    "Decides whether or not cutting planes are separated at the nodes of the branch-and-bound tree, "
    "choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_Nodecuts, values_bqpcuts, 0);

  AddSolverOption("cut:rlt rltcuts",
    "Whether to use the Relaxation Linearization Technique (RLT) to generate cuts in solving MIPs, choices as for \"cut:cuts\".",
    CPXPARAM_MIP_Cuts_RLT, values_bqpcuts, 0);

  AddSolverOption("cut:zerohalf zerohalfcuts",
    "Whether to use zero-half cuts in solving MIPs, choices as for \"cut:flowcover\".",
    CPXPARAM_MIP_Cuts_ZeroHalfCut, values_cuts2, 0);

  AddStoredOption("cut:cuts cuts mipcuts",
    "Global cut generation control, valid unless overridden "
    "by individual cut-type controls:\n"
    "\n.. value-table::\n",
    storedOptions_.cuts_, values_bqpcuts);

  AddSolverOption("cut:factor cutsfactor cutfactor",
    " Limit on MIP cuts added:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Limits_CutsFactor, values_cutsfactor, 0);

  AddSolverOption("cut:passes cutpasses cutpass",
    "Maximum number of cutting-plane passes during root-cut generation:\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_Limits_CutPasses, values_cutpasses, 0);

  AddStoredOption("cut:stats cutstats",
    " Whether the solve_message report the numbers and kinds of cuts used.\n"
    "\n.. value-table::\n",
    storedOptions_.cutstats_, values_01_noyes_0default_);

  // MIP related options

  AddSolverOption("mip:backtrack backtrack",
    "Tolerance (>0, default 0.9999) for when to backtrack during "
    "branch & bound.  Low values tend to pure best-bound search. "
    "High values(~1) tend to pure depth-first search. Values less "
    "than the default are often good when subproblems are expensive.",
    CPXPARAM_MIP_Strategy_Backtrack, 0.0, 1.0);

  AddSolverOption("mip:bbinterval bbinterval",
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

  AddSolverOption("mip:boundstr boundstr bndstrenind",
    "Whether to apply bound strengthening in mixed integer programs (MIPs):\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_BoundStrength, 
    values_autonoyes_, -1);

  AddSolverOption("mip:branchdir branchdir branch",
    "Which child node to explore first when branching:\n"
    "\n.. value-table::",
    CPXPARAM_MIP_Strategy_Branch, values_branchdir, 0);



  AddSolverOption("mip:gapabs mipgapabs absmipgap",
    "Max. absolute MIP optimality gap (default 1e-6).",
    CPXPARAM_MIP_Tolerances_AbsMIPGap, 0.0, DBL_MAX);

  AddSolverOption("mip:gap mipgap",
    "Max. relative MIP optimality gap (default 1e-4).",
    CPXPARAM_MIP_Tolerances_MIPGap, 1e-9, 1.0);


  AddSolverOption("mip:fpheur fpheur",
    "Whether to use the feasibility pump heuristic on MIP "
    "problems:\n\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_FPHeur, values_fpheur, 0);

  AddSolverOption("mip:inttol inttol intfeastol integrality",
    "Feasibility tolerance for integer variables (default 1e-05, "
    "must be in [0.0, 0.5])",
    CPXPARAM_MIP_Tolerances_Integrality, 0.0, 0.5);


  AddSolverOption("mip:nodefile nodefile",
    "Whether to save node information in a temporary file:\n"
    "\n.. value-table::\n", CPXPARAM_MIP_Strategy_File, 
    values_nodefile, 0);



  AddSolverOption("mip:search mipsearch",
    "Search strategy for mixed-integer problems:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_Search, values_mipsearch, 0);
  

  AddSolverOption("mip:submipalg submipalg",
    "Choice of algorithm used to solve the subproblems of a subMIP: "
    "not a subproblem, but an auxiliary MIP that CPLEX sometimes forms "
    "and solves, e.g., when dealing with a partial MIP start, "
    "repairing an infeasible MIP start, using the RINS heuristic, "
    "branching locally or polishing a solution. "
    "Possible values (only 0 is allowed for MIQCPs):\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_SubMIP_SubAlg, values_nodemethod, 0);

  AddSolverOption("mip:submipscale submipscale",
    "Rarely used choice of scaling for auxiliary subMIPs "
    "(described with \"submipalg\"):\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_SubMIP_Scale, values_prescale, 0);


  AddSolverOption("mip:submipstartalg submipstartalg",
    "Rarely used choice of algorithm for the initial relaxation of a "
    "subMIP (described with \"submipalg\"):\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_SubMIP_StartAlg, values_submipstartalg, 0);


  AddSolverOption("mip:nodemethod nodemethod mipalg mipalgorithm",
    "Algorithm used to solve relaxed MIP node problems; for MIQP problems "
    "(quadratic objective, linear constraints), settings other than 3 and 5 " 
    "are treated as 0. For MIQCP problems (quadratic objective and "
		"constraints), only 0 is permitted.\n"
    "\n.. value-table::\n", CPXPARAM_MIP_Strategy_SubAlgorithm, values_nodemethod, 0);

  AddSolverOption("mip:nodesel nodesel nodeselect",
    "Strategy for choosing next node while optimizing\n\
		integer variables:\n\n.. value-table::",
    CPXPARAM_MIP_Strategy_NodeSelect, values_nodeselect, 0);

  AddSolverOption("mip:varbranch varbranch varsel varselect" ,
    "MIP branch variable selection strategy:\n"
    "\n.. value-table::\n", 
    CPXPARAM_MIP_Strategy_VariableSelect, values_varselect, 0);

  AddSolverOption("bar:baralg baralg",
    "How to start the barrier algorithm:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_Algorithm, values_baralg, 0);

  AddSolverOption("bar:comptol comptol",
    "Convergence tolerance for barrier algorithm: "
		"the algorithm stops when the relative "
		"complementarity is < bartol (default 1e-8).",
    CPXPARAM_Barrier_ConvergeTol, 1e-12, DBL_MAX);

  AddStoredOption("bar:crossover crossover mipcrossover",
    "How to transform a barrier solution to a basic one:\n"
   "\n.. value-table::\n", storedOptions_.crossover_, values_barcrossover);

  AddSolverOption("bar:corrections barcorr barmaxcor",
    "Limit on centering corrections in each iteration of the barrier algorithm:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_Limits_Corrections, values_barmaxcor, -1);

  AddSolverOption("bar:dense bar:densecol dense densecol",
    "If positive, minimum nonzeros in a column for the barrier algorithm "
    "to consider the column dense. If 0 (default), this tolerance is "
    "selected automatically.", 
    CPXPARAM_Barrier_ColNonzeros, 0, INT_MAX);


  AddSolverOption("bar:display bardisplay",
    "Specifies how much the barrier algorithm chatters:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_Display, values_bardisplay, 1);

  AddSolverOption("bar:growth bargrowth growth",
    "Tolerance for detecting unbounded faces in the barrier algorithm: "
    "higher values make the test for unbounded faces harder to satisfy "
		"(default 1e12).",
    CPXPARAM_Barrier_Limits_Growth, 1, CPXINT_MAX);

  AddSolverOption("bar:iterlim bariterlim",
    "Maximum barrier iterations allowed (default 9223372036800000000).",
    CPXPARAM_Barrier_Limits_Iteration, 0, CPXINT_MAX);

  AddSolverOption("bar:objrange barobjrange",
    "Limit on the absolute objective value before the barrier algorithm "
    "considers the problem unbounded (default 1e20).",
    CPXPARAM_Barrier_Limits_ObjRange, 0, CPXINT_MAX);

  AddSolverOption("bar:start barstart barstartalg",
    "Barrier starting-point algorithm:"
    "\n\n.. value-table::\n",
    CPXPARAM_Barrier_StartAlg, values_barstart, 1);

  AddStoredOption("lp:solutiontype solutiontype",
    "Whether to seek a basic solution when solving an LP:\n"
    "\n.. value-table::\n", storedOptions_.solutionType_, values_solutiontype);

  // Solution pool controls
  AddSolverOption("sol:poolgap ams_eps poolgap",
    "Relative tolerance for reporting alternate MIP solutions "
    "(default: 1e75).", CPX_PARAM_SOLNPOOLGAP, 0.0, DBL_MAX);
  AddSolverOption("sol:poolgapabs ams_epsabs poolagap poolgapabs",
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

  AddSolverOption("tech:auxrootthreads auxrootthreads",
    "Controls the number of threads used for auxiliary chores when solving "
    "the root node of a MIP problem. When N threads are available (possibly "
    "limited by \"threads\"), auxrootthreads must be less than N.\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Limits_AuxRootThreads, auxrootthreads_values_, 0);

  
  AddSolverOption("alg:conflictalg conflictalg",
    "Choice of algorithm used by the CPLEX's conflict "
		"refiner:\n"
    "\n.. value-table::\n"
    "\nSettings 1, 2, and 3 are fast but may not discard "
    "many constraints; 5 and 6 work harder at this. "
    "Setting 4 searches for an Irreducible Infeasible "
    "Set of linear constraints(e.g., ignoring quadratic "
    "constraints).",
    CPXPARAM_Conflict_Algorithm, conflictalg_values, 0);

  AddSolverOption("alg:conflictdisplay conflictdisplay",
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

  AddSolverOption("tech:seed seed",
    "Seed for random number generator used internally "
    "by CPLEX.Use \"seed=?\" to see the default, which "
    "depends on the CPLEX release.",
    CPX_PARAM_RANDOMSEED, 0, CPXINT_MAX);

  AddSolverOption("tech:threads threads",
      "How many threads to use when using the barrier algorithm\n"
      "or solving MIP problems; default 0 ==> automatic choice.",
      CPXPARAM_Threads, 0, CPXINT_MAX);

  AddStoredOption("tech:workdir nodefiledir workfiledir workdir",
    "Directory where CPLEX creates a temporary subdirectory "
    "for temporary files, e.g., for node information and Cholesky factors.",
    storedOptions_.workDir_
  );

  AddSolverOption("tech:workfilelim workfilelim",
    "Maximum size in megabytes for in-core work \"files\". Default 2048.",
    CPX_PARAM_WORKMEM, 2048, 2048);



  AddSolverOption("lim:iter iterlim iterlimit iterations",
    "LP iteration limit (default:  9223372036800000000).",
    CPXPARAM_Simplex_Limits_Iterations, 0, CPXINT_MAX);

  AddSolverOption("lim:netiterations netiterations",
    "Limit on network simplex iterations (default:  9223372036800000000).",
    CPXPARAM_Network_Iterations, 0, CPXINT_MAX);

  AddSolverOption("lim:nodes node nodelim nodelimit",
    "Maximum MIP nodes to explore (default: 2^31 - 1).",
    CPX_PARAM_NODELIM, 0, CPXINT_MAX);

  AddSolverOption("lim:dettime dettimelim",
    "Time limit in platform-dependent \"ticks\".",
    CPXPARAM_DetTimeLimit, 0.0, DBL_MAX);

  AddSolverOption("lim:sol sollimit solutionlimit mipsolutions",
    "Limit the number of feasible MIP solutions found, causing early "
    "termination if exceeded; default = 2e31-1",
    CPXPARAM_MIP_Limits_Solutions , 0, 2000000000);

  AddSolverOption("lim:time timelim timelimit time",
      "limit on solve time (in seconds; default: no limit).",
      CPXPARAM_TimeLimit, 0.0, DBL_MAX);

  AddSolverOption("qp:target optimalitytarget",
    "Type of solution to compute for a (MI)QP (not QCP) problem:"
                  "\n\n.. value-table::\n",
    CPXPARAM_OptimalityTarget, optimalitytarget_values_, 0);

  
  AddSolverOption("pre:aggfill aggfill agglim",
    "Variables that appear in more than agglim rows (default 10) "
    "will not be substituted away by the aggregator.",
    CPXPARAM_Preprocessing_Fill, 0, CPXINT_MAX);

  AddSolverOption("pre:aggregate aggregate",
    "Whether to make substitutions to reduce the number of "
		"rows:\n\n.. value-table::\n",
    CPXPARAM_Preprocessing_Aggregator, values_preaggregator, -1);

  AddSolverOption("pre:coeffreduce coeffreduce",
    "Whether to use coefficient reduction when "
    "preprocessing MIPS\n\n.. value-table::\n",
    CPXPARAM_Preprocessing_CoeffReduce, values_precoeffreduce, -1);

  AddSolverOption("pre:dependency dependency",
    "Whether to use CPLEX's presolve dependency checker:\n\n.. value-table::\n",
    CPXPARAM_Preprocessing_Dependency, values_predependency, -1);

  AddSolverOption("pre:solve presolve", "Whether to use CPLEX's presolve:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Presolve, values_01_noyes_1default_, 1);

  AddSolverOption("pre:dual predual",
    "Whether CPLEX's presolve phase should present the "
    "CPLEX solution algorithm with the primal(-1) or "
    "dual(1) problem or (default = 0) should decide"
    "which automatically.Specifying \"predual=1\" often "
    "gives better performance than specifying just \"dual\", "
    "but sometimes \"dual predual=1\" is still better.",
    CPXPARAM_Preprocessing_Dual, -1, 1);

  AddSolverOption("pre:node presolvenode",
    "Whether to run presolve at each node of the MIP branch-and-bound:\n"
    "\n.. value-table::\n",
    CPXPARAM_MIP_Strategy_PresolveNode, values_noautoyes_, 0);

  AddSolverOption("pre:passes prepasses prepass",
    "Limit on the number of CPLEX presolve passes:\n"
    "\n"
    "| -1 - Automatic choice (default)\n"
    "| n>=0 - At most n passes.",
    CPXPARAM_Preprocessing_NumPass, -1, CPXINT_MAX);

  AddSolverOption("pre:reformulations prereformulations",
    "Kinds of reductions permitted during CPLEX presolve: \n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Reduce, values_prereformulations, 3);

  AddSolverOption("pre:relax prerelax",
    "Whether to use CPLEX's presolve on the initial LP relaxation of a MIP: \n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_Relax, values_autonoyes_, 0);


  AddSolverOption("pre:scale scale",
    "How to scale the problem:\n"
    "\n.. value-table::\n",
    CPXPARAM_Read_Scale, values_prescale, 0);

  AddSolverOption("pre:sos1enc presos1enc presos1reform",
    "Encoding used for SOS1 reformulation:\n"
    "\n.. value-table::\n",
    CPXPARAM_Preprocessing_SOS1Reform, values_presosenc, 0);

  AddSolverOption("pre:sos2enc presos2enc presos2reform",
    "Encoding used for SOS2 reformulation, see pre:sos1enc.",
    CPXPARAM_Preprocessing_SOS2Reform, -1, 1);

  AddStoredOption("tech:nosolve nosolve",
    "Stop after loading the problem and honoring any "
    "\"writeprob\" or \"writemipstart\" directives.",
    storedOptions_.noSolve_);


  
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
  if (objnparam_int_.rend() == it)
    throw std::runtime_error("Failed to find recorded option " +
      opt.wc_key_last__std_form());
  return it->second;
}
double CplexBackend::CplexGetObjDblParam(const SolverOption& opt) const {
  auto it = std::find_if(objnparam_dbl_.rbegin(), objnparam_dbl_.rend(),
    [&](const ObjNParam<int>& prm) {
      return prm.first == std::make_pair(opt.wc_tail(), opt.wc_keybody_last());
    });
  if (objnparam_dbl_.rend() == it)
    throw std::runtime_error("Failed to find recorded option " +
      opt.wc_key_last__std_form());
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

#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"

#include "gurobibackend.h"

extern "C" {
  #include "gurobi-ampls-c-api.h"    // Gurobi AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {

bool InterruptGurobi(void *model) {
  GRBterminate( static_cast<GRBmodel*>(model) );
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateGurobiBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::GurobiBackend()};
}

namespace mp {

/// Create Gurobi Model Manager
/// @param gc: the Gurobi Backend
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return GurobiModelMgr
std::unique_ptr<BasicModelManager>
CreateGurobiModelMgr(GurobiCommon& gc,
                     Env& e, pre::BasicPresolver*& pPre);


GurobiBackend::GurobiBackend() {
  pre::BasicPresolver* pPre;
  auto data = CreateGurobiModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetPresolver(pPre);
}

GurobiBackend::~GurobiBackend() {
  CloseGurobi();
}

std::string GurobiBackend::GetSolverVersion() {
  int a,b,c;
  GRBversion(&a, &b, &c);
  return fmt::format("{}.{}.{}", a, b, c);
}

void GurobiBackend::CloseGurobi() {
  /* Free the fixed model */
  if (has_model() && model() != model_fixed_) {
    assert(model_fixed_);
    GRBfreemodel(model_fixed_);
  }
  model_fixed_ = nullptr;

  /* Free model */
  if (has_model()) {
    GRBfreemodel(model());
    model_ref() = nullptr;
  }

  /* Free environment */
  if (has_env()) {
    GRBfreeenv(env());
    env_ref() = nullptr;
  }
}

void GurobiBackend::InitOptionParsing() {
  OpenGurobi();
}

void GurobiBackend::OpenGurobi() {
  // Typically try the registered function first;
  // if not available call the solver's API function to create
  // an empty environment. Will be started later.
  const auto& create_fn = GetCallbacks().cb_initsolver_;
  if (create_fn)
    set_env((GRBenv*)create_fn());
  else 
    GRB_CALL(GRBemptyenv(&env_ref()));
  /* Set default parameters */
  GRBsetintparam(env(), GRB_INT_PAR_OUTPUTFLAG, 0);
}

void GurobiBackend::OpenGurobiModel() {
  /* Create an empty model */
  GRB_CALL( GRBnewmodel(env(), &model_ref(), "amplgurobi", 0,
                        NULL, NULL, NULL, NULL, NULL) );
  assert(model());
  /* Init fixed model */
  model_fixed_ = model();

  /* Copy handlers to ModelAPI */
  copy_common_info_to_other();
}

void GurobiBackend::FinishOptionParsing() {
  if (servers().size()) {
    OpenGurobiComputeServer();
  }
  else if (cloudid().size() && cloudkey().size()) {
    OpenGurobiCloud();
  }
  else {
    // If a user defined function had been provided, the environment is assumed
    // as already started
    if (!GetCallbacks().cb_initsolver_)
    {
      int res = GRBstartenv(env_ref());
      if (res)
      {
        const auto& diag = GetCallbacks().cb_diagnostics_;
        if (diag)
          diag();
        else {
          MP_RAISE(
            fmt::format("Start environment failed with code {}, Gurobi message:\n{}",
              res, GRBgeterrormsg(env())));
        }
        exit(res);
      }
    }
  }
  OpenGurobiModel();
  if (paramfile_read().size())
    GRB_CALL(
          GRBreadparams(GRBgetenv(model()),
                        paramfile_read().c_str() ));
  if (paramfile_write().size())
    GRB_CALL(
          GRBwriteparams(GRBgetenv(model()),
                         paramfile_write().c_str() ));
  /// Tell the base class our verbosity
  set_verbose_mode(GrbGetIntParam(GRB_INT_PAR_OUTPUTFLAG));
}

void GurobiBackend::OpenGurobiComputeServer() {
  assert(servers().size());
  SetSolverOption(GRB_STR_PAR_COMPUTESERVER, servers().c_str());
  SetSolverOption(GRB_STR_PAR_SERVERPASSWORD, server_password().c_str());
  SetSolverOption(GRB_STR_PAR_CSROUTER, server_router().c_str());
  SetSolverOption(GRB_STR_PAR_CSGROUP, server_group().c_str());
  SetSolverOption(GRB_INT_PAR_CSTLSINSECURE, server_insecure());
  SetSolverOption(GRB_INT_PAR_CSPRIORITY, server_priority());
  SetSolverOption(GRB_INT_PAR_SERVERTIMEOUT, server_timeout());

  if (int i = GRBstartenv(env()))
    switch (i) {
    case GRB_ERROR_NETWORK:
      Abort(601, "Could not talk to Gurobi Compute Server(s).");
      break;
    case GRB_ERROR_JOB_REJECTED:
      Abort(602, "Job rejected by Gurobi Compute Server(s).");
      break;
    case GRB_ERROR_NO_LICENSE:
      Abort(603, "No license for specified Gurobi Compute Server(s).");
      break;
    default:
      Abort(604, fmt::format(
        "Surprise return {} while starting the compute server environment.", i));
    }
}


void GurobiBackend::OpenGurobiCloud() {
  assert(cloudid().size() && cloudkey().size());
  SetSolverOption(GRB_STR_PAR_CLOUDACCESSID, cloudid().c_str());
  SetSolverOption(GRB_STR_PAR_CLOUDSECRETKEY, cloudkey().c_str());
  SetSolverOption(GRB_STR_PAR_CLOUDPOOL, cloudpool().c_str());
  SetSolverOption(GRB_INT_PAR_CSPRIORITY, cloudpriority());

  if (int i = GRBstartenv(env())) {
    switch(i) {
    case GRB_ERROR_NETWORK:
      Abort(601, "Could not talk to Gurobi Instant Cloud.");
      break;
    case GRB_ERROR_JOB_REJECTED:
      Abort(602, "Job rejected by Gurobi Instant Cloud.");
      break;
    case GRB_ERROR_NO_LICENSE:
      Abort(603, "No license for specified Gurobi Instant Cloud.");
      break;
    case GRB_ERROR_CLOUD:
      Abort(605, "Bad value for cloudid or cloudkey, or Gurobi Cloud out of reach.");
      break;
    default:
      Abort(604, fmt::format(
              "Surprise return {} while starting the cloud environment", i));
    }
  }
}

bool GurobiBackend::IsMIP() const {
  return 1 == GrbGetIntAttr(GRB_INT_ATTR_IS_MIP);
}

bool GurobiBackend::IsQP() const {
  return 1 == GrbGetIntAttr(GRB_INT_ATTR_IS_QP);
}

bool GurobiBackend::IsQCP() const {
  return 1 == GrbGetIntAttr(GRB_INT_ATTR_IS_QCP);
}


std::vector<double> GurobiBackend::PrimalSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_X, NumVars());
}

Solution GurobiBackend::GetSolution() {
  auto mv = GetPresolver().PostsolveSolution(
        { PrimalSolution(), DualSolution() } );
  return { mv.GetVarValues()(), mv.GetConValues()(),
    GetObjectiveValues() };    // TODO postsolve obj values
}

pre::ValueMapDbl GurobiBackend::DualSolution() {
  return {{
      { CG_Linear, GurobiDualSolution_LP() },
      { CG_Quadratic, GurobiDualSolution_QCP() } }};
}

std::vector<double> GurobiBackend::GurobiDualSolution_LP() {
  return
    GrbGetDblAttrArray_VarCon(model_fixed_, GRB_DBL_ATTR_PI, 1);
}

std::vector<double> GurobiBackend::GurobiDualSolution_QCP() {
  return
    GrbGetDblAttrArray(model_fixed_, GRB_DBL_ATTR_QCPI, NumQPCons());
}

double GurobiBackend::ObjectiveValue() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_OBJVAL);
}

ArrayRef<double> GurobiBackend::GetObjectiveValues() {
  int no = NumObjs();
  if(no==0)
    return { };
  std::vector<double> objs(no,
                           std::numeric_limits<double>::quiet_NaN());
  bool fOk = true;
  if (NumObjs() == 1)
    objs[0] = GrbGetDblAttr(GRB_DBL_ATTR_OBJVAL, &fOk); // failsafe
  else {
    GRBenv* mdlenv = GRBgetenv(model());
    int objnumber = GrbGetIntParam(GRB_INT_PAR_OBJNUMBER);
    for (int i = 0; i < no; i++)
    {
      GRB_CALL( GRBsetintparam(mdlenv, GRB_INT_PAR_OBJNUMBER, i) );
      objs[i] = GrbGetDblAttr(GRB_DBL_ATTR_OBJNVAL, &fOk);
    }
    GRB_CALL( GRBsetintparam(mdlenv, GRB_INT_PAR_OBJNUMBER, objnumber) );
  }
  return objs;
}

ArrayRef<double> GurobiBackend::CurrentGrbPoolPrimalSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_XN, NumVars());
}

double GurobiBackend::CurrentGrbPoolObjectiveValue() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_POOLOBJVAL);
}


void GurobiBackend::MarkLazyOrUserCuts(ArrayRef<int> lazyVals) {
  auto vm = GetPresolver().PresolveLazyUserCutFlags({
                                                      {}, lazyVals
                                                    });
  std::vector<int> lv_lin = vm.GetConValues()(CG_Linear);
  if (lv_lin.size()) {
    for (size_t i=0; i<lv_lin.size(); ++i) {
      if (auto& val = lv_lin[i]) {
        if ((val>0 && !lazy_cuts()) ||
            (val<0 && !user_cuts()))
          val = 0;
      }
    }
    GRB_CALL( GRBsetintattrarray(model(),
                                 GRB_INT_ATTR_LAZY,
                                 0, (int)lv_lin.size(), (int*)lv_lin.data()) );
    // Testing
    ReportFirstLinearConstraintLazySuffix(lv_lin[0]);
  }
  // Testing
  const auto& lv_quad = vm.GetConValues()(CG_Quadratic);
  if (lv_quad.size())
    ReportFirstQCLazySuffix(lv_quad[0]);
}

SolutionBasis GurobiBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetPresolver().PostsolveBasis(
          { std::move(varstt),
            {{{ CG_Linear, std::move(constt) }}} } );
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
    assert(constt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void GurobiBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetPresolver().PresolveBasis(
        { basis.varstt, basis.constt } );
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

ArrayRef<int> GurobiBackend::VarStatii() {
  auto stt =
    GrbGetIntAttrArray(model_fixed_,
        GRB_INT_ATTR_VBASIS, NumVars());
  for (auto& s: stt) {
    switch (s) {
    case 0:
      s = (int)BasicStatus::bas;
      break;
    case -1:
      s = (int)BasicStatus::low;
      break;
    case -2:
      s = (int)BasicStatus::upp;
      break;
    case -3:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Gurobi VBasis value: {}", s));
    }
  }
  return stt;
}

ArrayRef<int> GurobiBackend::ConStatii() {
  auto stt =
    GrbGetIntAttrArray(model_fixed_,
        GRB_INT_ATTR_CBASIS, NumLinCons());
  for (auto& s: stt) {
    switch (s) {
    case 0:
      s = (int)BasicStatus::bas;
      break;
    case -1:
      s = (int)BasicStatus::sup;   // TODO exact value low/upp/equ??
      break;
    default:
      MP_RAISE(fmt::format("Unknown Gurobi CBasis value: {}", s));
    }
  }
  return stt;
}

void GurobiBackend::VarStatii(ArrayRef<int> vst) {
  std::vector<int> stt(vst.data(), vst.data()+vst.size());
  for (auto j=stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = 0;
      break;
    case BasicStatus::low:
    case BasicStatus::equ:
      s = -1;
      break;
    case BasicStatus::upp:
      s = -2;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = -3;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      if (!GRBgetdblattrelement(model(), GRB_DBL_ATTR_LB, j, &lb) &&
          !GRBgetdblattrelement(model(), GRB_DBL_ATTR_UB, j, &ub)) {
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
  GrbSetIntAttrArray(GRB_INT_ATTR_VBASIS, stt);
}

void GurobiBackend::ConStatii(ArrayRef<int> cst) {
  std::vector<int> stt(cst.data(), cst.data()+cst.size());
  for (auto& s: stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = 0;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status,
    case BasicStatus::low:    // as Gurobi 9.5 does not accept partial basis.
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = -1;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  GrbSetIntAttrArray(GRB_INT_ATTR_CBASIS, stt);
}

void GurobiBackend::AddPrimalDualStart(Solution sol0) {
  auto mv = GetPresolver().PresolveSolution(
        { sol0.primal, sol0.dual } );
  auto x0 = mv.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);
  GrbSetDblAttrArray(GRB_DBL_ATTR_PSTART, x0);
  GrbSetDblAttrArray(GRB_DBL_ATTR_DSTART, pi0);
}

void GurobiBackend::AddMIPStart(ArrayRef<double> x0) {
  switch (Gurobi_mipstart()) {
  case 0: break;
  case 1:
    GrbSetDblAttrArray(GRB_DBL_ATTR_START, x0);
    break;
  case 3:
    GrbSetIntAttrArray(GRB_INT_ATTR_VARHINTPRI,
                         ReadSuffix(sufHintPri));
    GrbSetDblAttrArray(GRB_DBL_ATTR_VARHINTVAL, x0);
    break;
  case 2:
    GrbSetDblAttrArray(GRB_DBL_ATTR_VARHINTVAL, x0);
    break;
  default:
    assert(0);
  }
}

void GurobiBackend::VarPriorities(ArrayRef<int> priority) {
  GrbSetIntAttrArray(GRB_INT_ATTR_BRANCHPRIORITY, priority);
}

void GurobiBackend::ObjPriorities(ArrayRef<int> priority) {
  for (int i=0; i<(int)priority.size(); ++i) {
    GrbSetIntParam(GRB_INT_PAR_OBJNUMBER, i);
    GrbSetIntAttr(GRB_INT_ATTR_OBJNPRIORITY, priority[i]);
  }
}

void GurobiBackend::ObjWeights(ArrayRef<double> val) {
  for (int i=0; i<(int)val.size(); ++i) {
    GrbSetIntParam(GRB_INT_PAR_OBJNUMBER, i);
    GrbSetDblAttr(GRB_DBL_ATTR_OBJNWEIGHT, val[i]);
  }
}

void GurobiBackend::ObjAbsTol(ArrayRef<double> val) {
  for (int i=0; i<(int)val.size(); ++i) {
    GrbSetIntParam(GRB_INT_PAR_OBJNUMBER, i);
    GrbSetDblAttr(GRB_DBL_ATTR_OBJNABSTOL, val[i]);
  }
}

void GurobiBackend::ObjRelTol(ArrayRef<double> val) {
  for (int i=0; i<(int)val.size(); ++i) {
    GrbSetIntParam(GRB_INT_PAR_OBJNUMBER, i);
    GrbSetDblAttr(GRB_DBL_ATTR_OBJNRELTOL, val[i]);
  }
}


ArrayRef<double> GurobiBackend::Ray() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_UNBDRAY, NumVars());
}

ArrayRef<double> GurobiBackend::DRay() {
  auto fd = GrbGetDblAttrArray(GRB_DBL_ATTR_FARKASDUAL, NumLinCons());
  auto vm = GetPresolver().PostsolveSolution({
                                               {},
                                               {{{CG_Linear, std::move(fd)}}}
                                             });
  return vm.GetConValues().MoveOut(0);   // need the vector itself
}

IIS GurobiBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();

  auto mv = GetPresolver().PostsolveIIS(
        { variis, coniis } );
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> GurobiBackend::VarsIIS() {
  auto iis_lb =
    GrbGetIntAttrArray(GRB_INT_ATTR_IIS_LB, NumVars());
  auto iis_ub =
    GrbGetIntAttrArray(GRB_INT_ATTR_IIS_UB, NumVars());
  for (size_t i = iis_lb.size(); i--; ) {
    if (iis_ub[i]) {
      if (iis_lb[i])
        iis_lb[i] = (int)IISStatus::fix;
      else
        iis_lb[i] = (int)IISStatus::upp;
    } else {
      if (iis_lb[i])
        iis_lb[i] = (int)IISStatus::low;
      else
        iis_lb[i] = (int)IISStatus::non;
    }
  }
  return iis_lb;
}

static void ConvertIIS2AMPL(std::vector<int>& ai) {
  for (int i=ai.size(); i--; ) {
    ai[i] = int(ai[i] ? IISStatus::mem : IISStatus::non);
  }
}

pre::ValueMapInt GurobiBackend::ConsIIS() {
  std::vector<int> iis_lincon =
      GrbGetIntAttrArray(GRB_INT_ATTR_IIS_CONSTR, NumLinCons());
  ConvertIIS2AMPL(iis_lincon);
  std::vector<int> iis_qc =
      GrbGetIntAttrArray(GRB_INT_ATTR_IIS_QCONSTR, NumQPCons());
  ConvertIIS2AMPL(iis_qc);
  std::vector<int> iis_soscon =
      GrbGetIntAttrArray(GRB_INT_ATTR_IIS_SOS, NumSOSCons());
  ConvertIIS2AMPL(iis_soscon);
  std::vector<int> iis_gencon =
      GrbGetIntAttrArray(GRB_INT_ATTR_IIS_GENCONSTR, NumGenCons());
  ConvertIIS2AMPL(iis_gencon);
  return {{{ CG_Linear, iis_lincon },
        { CG_Quadratic, iis_qc },
      { CG_SOS, iis_soscon },
      { CG_General, iis_gencon }}};
}

double GurobiBackend::MIPGap() {
  bool f;
  double g = GrbGetDblAttr(GRB_DBL_ATTR_MIPGAP, &f);
  return f ? g : Infinity();
}

double GurobiBackend::MIPGapAbs() {
    return std::fabs(
          ObjectiveValue() - BestDualBound() );
}

double GurobiBackend::BestDualBound() {
  bool f;
  double g = GrbGetDblAttr(GRB_DBL_ATTR_OBJBOUND, &f);
  return f ? g : -ModelSense() * Infinity();
}

double GurobiBackend::Kappa() {
  return GrbGetDblAttr(GRB_DBL_ATTR_KAPPA);
}


ArrayRef<double> GurobiBackend::Senslbhi() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SALBUp", 0); }
ArrayRef<double> GurobiBackend::Senslblo() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SALBLow", 0); }
ArrayRef<double> GurobiBackend::Sensobjhi() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SAObjUp", 0); }
ArrayRef<double> GurobiBackend::Sensobjlo() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SAObjLow", 0); }
ArrayRef<double> GurobiBackend::Sensrhshi() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SARHSUp", 1); }
ArrayRef<double> GurobiBackend::Sensrhslo() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SARHSLow", 1); }
ArrayRef<double> GurobiBackend::Sensubhi() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SAUBUp", 0); }
ArrayRef<double> GurobiBackend::Sensublo() const
{ return GrbGetDblAttrArray_VarCon(model_fixed_, "SAUBLow", 0); }



double GurobiBackend::NodeCount() const {
  bool f;
  return GrbGetDblAttr(GRB_DBL_ATTR_NODECOUNT, &f);
}

double GurobiBackend::SimplexIterations() const {
  bool f;
  return GrbGetDblAttr(GRB_DBL_ATTR_ITERCOUNT, &f);
}

int GurobiBackend::BarrierIterations() const {
  bool f;
  return GrbGetIntAttr(GRB_INT_ATTR_BARITERCOUNT, &f);
}

void GurobiBackend::ExportModel(const std::string &file) {
  GRB_CALL( GRBwrite(model(), file.c_str()) );
}

void GurobiBackend::InputExtras() {
  BaseBackend::InputExtras();
  InputGurobiExtras();
}

void GurobiBackend::InputGurobiExtras() {
  if (need_multiple_solutions())
    GrbSetIntParam(GRB_INT_PAR_POOLSEARCHMODE, storedOptions_.nPoolMode_);
  if (need_ray_primal() || need_ray_dual())
    GrbSetIntParam(GRB_INT_PAR_INFUNBDINFO, 1);
  GrbPlayObjNParams();
  if (feasrelax())
    DoGurobiFeasRelax();
  SetPartitionValues();
}

void GurobiBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptGurobi, model());
}


///////////////////////////////////// SOLVE /////////////////////////////////////////
void GurobiBackend::Solve() {
  PrepareGurobiSolve();

  GRB_CALL( GRBoptimize(model()) );

  WindupGurobiSolve();
}

void GurobiBackend::PrepareGurobiSolve() {
  /// After all attributes applied
  if (!storedOptions_.exportFile_.empty())
    ExportModel(storedOptions_.exportFile_);
  if (tunebase().size())
    DoGurobiTune();
}

void GurobiBackend::WindupGurobiSolve() { }

void GurobiBackend::ReportResults() {
  ReportGurobiResults();
  BaseBackend::ReportResults();
}

void GurobiBackend::ReportGurobiResults() {
  SetStatus( ConvertGurobiStatus() );
  AddGurobiMessage();
  if (need_multiple_solutions())
    ReportGurobiPool();
  if (need_fixed_MIP())
    ConsiderGurobiFixedModel();
}

void GurobiBackend::AddGurobiMessage() {
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

void GurobiBackend::DoGurobiTune() {
  assert(tunebase().size());
  GRB_CALL( GRBtunemodel(model()) );
//		TODO? solve_result_num = 532;
  auto n_results = GrbGetIntAttr(GRB_INT_ATTR_TUNE_RESULTCOUNT);
  if (n_results<=0)
    MP_RAISE("No tuning results!");
  auto tbc = tunebase();
  if (tbc.size()>=4 &&
      0==tbc.compare(tbc.size()-4, 4, ".prm"))
    tbc.resize(tbc.size()-4);
  tbc += "_{}_";
  tbc += ".prm";
  std::string tfn;
  for (int k=n_results; k--;) {
    GRB_CALL_MSG( GRBgettuneresult(model(), k),
      fmt::format(
        "Surprize return from GRBgettuneresult({})", k+1));
    tfn = fmt::format(tbc.c_str(), k+1);
    GRB_CALL_MSG( GRBwriteparams(GRBgetenv(model()), tfn.c_str()),
      fmt::format(
        "Surprize return from GRBwriteparams({})", tfn));
  }
  AddToSolverMessage(
        fmt::format("Tuning: wrote {} parameter files, best file: '{}'", n_results, tfn));
}

void GurobiBackend::ReportGurobiPool() {
  if (!IsMIP())         // Gurobi 9.1.2 returns 1 solution for LP
    return;             // but cannot retrieve its pool attributes
  int iPoolSolution = -1;
  while (++iPoolSolution < GrbGetIntAttr(GRB_INT_ATTR_SOLCOUNT)) {
    GrbSetIntParam(GRB_INT_PAR_SOLUTIONNUMBER, iPoolSolution);
    ReportIntermediateSolution(
          { CurrentGrbPoolPrimalSolution(),
            {}, { CurrentGrbPoolObjectiveValue() } });
  }
}

void GurobiBackend::ConsiderGurobiFixedModel() {
  if (!IsMIP())
    return;
  if (IsQCP()) {
    int i=0;
    if (GRBgetintparam(env(), GRB_INT_PAR_QCPDUAL, &i) || i == 0)
      return;
  }
  if (GRBmodel* mdl = GRBfixedmodel(model()))
    model_fixed_ = mdl;
  else
    return;
  auto msg = DoGurobiFixedModel();
  if (!msg.empty()) {
    AddToSolverMessage( msg +
                        " failed in DoGurobiFixedModel()." );
    GRBfreemodel(model_fixed_);
    model_fixed_ = model();
  }
}

std::string GurobiBackend::DoGurobiFixedModel() {
  GRBenv *env;
  if (!(env = GRBgetenv(model_fixed_)))
    return "GRBgetenv";
  if (GRBsetintparam(env, "Presolve", 0))     // why?
    return "GRBsetintparam(\"Presolve\")";
  int k = -12345;
  GRBgetintparam(env, GRB_INT_PAR_METHOD, &k);
  int& fixedmethod = storedOptions_.nFixedMethod_;
  if (fixedmethod < -1 || fixedmethod >= 5) { /* not specified or invalid */
    if (k >= 2 || k < 0)
      fixedmethod = 1;
    else
      fixedmethod = k;
    }
  if (fixedmethod != k)
    GRBsetintparam(env, GRB_INT_PAR_METHOD, k);
  if (!GRBgetintparam(env, GRB_INT_PAR_METHOD, &k) && k != fixedmethod)
    GRBsetintparam(env, GRB_INT_PAR_METHOD, fixedmethod);
  /// TODO output model(s)
  if (GRBoptimize(model_fixed_))
    return "optimize()";
  int i;
  if (GRBgetintattr(model_fixed_, GRB_INT_ATTR_STATUS, &i))
    return "getintattr()";
  static const char *statusname_from_infeasible[] = {
    "infeasible",
    "infeasible or unbounded",
    "unbounded",
    "cutoff",
    "iteration limit",
    "node limit",
    "time limit",
    "solution limit",
    "interrupted",
    "numeric difficulty",
    "suboptimal"
    };
  if (i != GRB_OPTIMAL) {
    if (i >= GRB_INFEASIBLE && i <= GRB_SUBOPTIMAL)
      return fmt::format(
            "Fixed model status: {}. GRBoptimize",
        statusname_from_infeasible[i-GRB_INFEASIBLE]);
    else
      return fmt::format(
            "Surprise status {} after GRBoptimize",
        i);
  }
  double f;
  if (!GRBgetdblattr(model_fixed_, GRB_DBL_ATTR_ITERCOUNT, &f)) {
    if (f)
      AddToSolverMessage( fmt::format(
            "Fixed MIP for mip:basis: {} simplex iteration{}",
            f, "s"[f == 1.]) );
  }
  return {};
}

void GurobiBackend::DoGurobiFeasRelax() {
  int reltype = feasrelax()-1,
      minrel = 0;
  if (reltype >= 3) {
    reltype -= 3;
    minrel = 1;
    feasrelax().flag_orig_obj_available();
  }
  /// Gurobi 9.5 only relaxes linear constraints
  auto mv = GetPresolver().PresolveSolution({
                                              {},
                                              feasrelax().rhspen()
                                            });
  const auto& rhspen = mv.GetConValues()(CG_Linear);
  /// Account for new variables (TODO: proper presolve)
  std::vector<double> lbpen = feasrelax().lbpen();
  if (lbpen.size() && lbpen.size()<(size_t)NumVars())
    lbpen.resize(NumVars());
  std::vector<double> ubpen = feasrelax().ubpen();
  if (ubpen.size() && ubpen.size()<(size_t)NumVars())
    ubpen.resize(NumVars());
  GRB_CALL( GRBfeasrelax(model(), reltype, minrel,
                         (double*)data_or_null( lbpen ),
                         (double*)data_or_null( ubpen ),
                         (double*)data_or_null( rhspen ),
                         &feasrelax().orig_obj_value()) );
}

void GurobiBackend::SetPartitionValues() {
  if (auto part = ReadIntSuffix( {"partition", suf::Kind::VAR} )) {
    GrbSetIntAttrArray(GRB_INT_ATTR_PARTITION, part);
  }
}


//////////////////////////////////////////////////////////////////////
////////////////////////// Solution Status ///////////////////////////
//////////////////////////////////////////////////////////////////////
std::pair<int, std::string> GurobiBackend::ConvertGurobiStatus() const {
  namespace sol = mp::sol;
  int optimstatus;
  GRB_CALL( GRBgetintattr(model(), GRB_INT_ATTR_STATUS, &optimstatus) );
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter()->Stop()) {
      return { sol::INTERRUPTED, "interrupted" };
    }
    int solcount;
    GRB_CALL( GRBgetintattr(model(), GRB_INT_ATTR_SOLCOUNT, &solcount) );
    if (solcount>0) {
      return { sol::UNCERTAIN, "feasible solution" };
    }
    return { sol::UNKNOWN, "unknown solution status" };
  case GRB_OPTIMAL:
    return { sol::SOLVED, "optimal solution" };
  case GRB_INFEASIBLE:
    return { sol::INFEASIBLE, "infeasible problem" };
  case GRB_INF_OR_UNBD:
    return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
  case GRB_UNBOUNDED:
    return { sol::UNBOUNDED, "unbounded problem" };
  case GRB_NUMERIC:
    return { sol::NUMERIC, "feasible or optimal but numeric issue" };
  case GRB_CUTOFF:
    return { sol::LIMIT, "objective cutoff" };
  case GRB_ITERATION_LIMIT:
    return { sol::LIMIT+1, "iteration limit" };
  case GRB_NODE_LIMIT:
    return { sol::LIMIT+2, "node limit" };
  case GRB_TIME_LIMIT:
    return { sol::LIMIT+3, "time limit" };
  case GRB_SOLUTION_LIMIT:
    return { sol::LIMIT+4, "solution limit" };
  case GRB_SUBOPTIMAL:
    return { sol::UNCERTAIN, "suboptimal" };
  case GRB_USER_OBJ_LIMIT:
    return { sol::UNCERTAIN+3, "bestobjstop or bestbndstop reached" };
  }
}


void GurobiBackend::ComputeIIS() {
  GRB_CALL(GRBcomputeIIS(model()));
  SetStatus( ConvertGurobiStatus() );   // could be new information
}



///////////////////////////////////////////////////////////////
////////////////////////// OPTIONS ////////////////////////////

// static possible values with descriptions

static const mp::OptionValueInfo values_barcrossover[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "None: return an interior solution", 0},
  { "1", "Push dual vars first, finish with primal simplex", 1},
  { "2", "Push dual vars first, finish with dual simplex", 2},
  { "3", "Push primal vars first, finish with primal simplex", 3},
  { "4", "Push primal vars first, finish with dual simplex.", 4}
};

static const mp::OptionValueInfo values_barcrossoverbasis[] = {
  { "0", "Favor speed (default)", 0},
  { "1", "Favor numerical stability.", 1}
};

static const mp::OptionValueInfo values_barhomogeneous[] = {
  {"-1", "Only when solving a MIP node relaxation (default)", -1},
  { "0", "Never", 0},
  { "1", "Always.", 1}
};

static const mp::OptionValueInfo values_barorder[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Approximate minimum degree", 0},
  { "1", "Nested dissection.", 1}
};

static const mp::OptionValueInfo values_bqpcuts[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Disallow BQP cuts", 0},
  { "1", "Enable moderate BQP cut generation", 1},
  { "2", "Enable aggressive BQP cut generation.", 2}
};

static const mp::OptionValueInfo values_branchdir[] = {
  {"-1", "Explore \"down\" branch first", -1},
  { "0", "Explore \"most promising\" branch first (default)", 0},
  { "1", "Explore \"up\" branch first.", 1}
};

static const mp::OptionValueInfo values_cuts[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "No cuts", 0},
  { "1", "Conservative cut generation", 1},
  { "2", "Aggressive cut generation", 2},
  { "3", "Very aggressive cut generation.", 3}
};
static constexpr int PrmCutsMin=-1, PrmCutsMax=3;
static const mp::OptionValueInfo values_cuts_upto2[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "No cuts", 0},
  { "1", "Conservative cut generation", 1},
  { "2", "Aggressive cut generation.", 2}
};

static const mp::OptionValueInfo values_disconnected[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "No", 0},
  { "1", "Moderate effort", 1},
  { "2", "Aggressive effort.", 2},
};

static const mp::OptionValueInfo values_iismethod[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Often faster than method 1", 0},
  { "1", "Can find a smaller IIS than method 0", 1},
  { "2", "Ignore the bound constraints.", 2},
};

static const mp::OptionValueInfo values_infproofcuts[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "No", 0},
  { "1", "Moderate cut generation", 1},
  { "2", "Aggressive cut generation.", 2},
};

static const mp::OptionValueInfo values_method[] = {
  { "-1", "Automatic (default): 3 for LP, 2 for QP, 1 for MIP", -1},
  { "0", "Primal simplex", 0},
  { "1", "Dual simplex", 1},
  { "2", "Barrier", 2},
  { "3", "Nondeterministic concurrent (several solves in parallel)", 3},
  { "4", "Deterministic concurrent", 4},
  { "5", "Deterministic concurrent simplex.", 5}
};

static const mp::OptionValueInfo values_mipfocus[] = {
  { "0", "Balance finding good feasible solutions and "
    "proving optimality (default)", 0},
  { "1", "Favor finding feasible solutions", 1},
  { "2", "Favor providing optimality", 2},
  { "3", "Focus on improving the best objective bound.", 3},
};

static const mp::OptionValueInfo values_mipstart_[4] = {
  {     "0", "No (overrides alg:start)", 0 },
  {     "1", "Yes (default)", 1},
  {     "2", "No, but use the incoming primal "
        "values as hints (VARHINTVAL), ignoring the .hintpri suffix", 2},
  {     "3", "Similar to 2, but use the .hintpri suffix on "
        "variables:  larger (integer) values give greater "
        "priority to the initial value of the associated "
        "variable.", 3}
};

static const mp::OptionValueInfo values_miqcpmethod[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Solve continuous QCP relaxations at each node", 0},
  { "1", "Use linearized outer approximations.", 1}
};

static const mp::OptionValueInfo values_multiobjmethod[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Primal simplex", 0},
  { "1", "Dual simplex", 1},
  { "2", "Ignore warm-start information; use the algorithm "
    "specified by the method keyword.", 2}
};

static const mp::OptionValueInfo values_multiobjpre[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Do not use Gurobi's presolve", 0},
  { "1", "Conservative presolve", 1},
  { "2", "Aggressive presolve, which may degrade lower priority objectives.", 2}
};

static const mp::OptionValueInfo values_nodemethod[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Primal simplex", 0},
  { "1", "Dual simplex", 1},
  { "2", "Barrier.", 2}
};

static const mp::OptionValueInfo values_nonconvex[] = {
  { "-1", "Default choice (currently the same as 1)", -1},
  { "0", "Complain about nonquadratic terms", 0},
  { "1", "Complain if Gurobi's presolve cannot discard or "
    "eliminate nonquadratic terms", 1},
  { "2", "Translate quadratic forms to bilinear form and use "
    "spatial branching.", 2}
};

static const mp::OptionValueInfo values_predeprow[] = {
  { "-1", "Only for continuous models (default)", -1},
  { "0", "Never", 0},
  { "1", "For all models.", 1}
};

static const mp::OptionValueInfo values_predual[] = {
  { "-1", "Automatic choice (default)", -1},
  { "0", "No", 0},
  { "1", "Yes", 1},
  { "2", "Form both primal and dual and use two threads to "
    "choose heuristically between them.", 2}
};

static const mp::OptionValueInfo values_premiqcpform[] = {
  { "-1", "Automatic choice (default)", -1},
  { "0", "Retain MIQCP form", 0},
  { "1", "Transform to second-order cone contraints", 1},
  { "2", "Transform to rotated cone constraints.", 2}
};

static const mp::OptionValueInfo values_preqlinearize[] = {
  { "-1", "Automatic choice (default)", -1},
  { "0", "Do not modify the quadratic part(s)\n"
    "\t\t\t 1 or 2 = try to linearize quadratic parts:", 0},
  { "1", "Focus on a strong LP relaxation", 1},
  { "2", "Focus on a compact LP relaxation.", 2}
};

static const mp::OptionValueInfo values_prescale[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "No", 0},
  { "1", "Yes", 1},
  { "2", "Yes, more aggressively", 2},
  { "3", "Yes, even more aggressively.", 3}
};
static const mp::OptionValueInfo values_pricing[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Partial pricing", 0},
  { "1", "Steepest edge", 1},
  { "2", "Devex", 2},
  { "3", "Quick-start steepest edge.", 3}
};

static const mp::OptionValueInfo values_pool_mode[] = {
  { "0", "Just collect solutions during normal solve, and sort them best-first", 0},
  { "1", "Make some effort at finding additional solutions", 1},
  { "2", "Seek \"poollimit\" best solutions (default)."
    "'Best solutions' are defined by the poolgap(abs) parameters.", 2}
};

static const mp::OptionValueInfo values_siftmethod_[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Primal simplex", 0},
  { "1", "Dual simplex", 1},
  { "2", "Barrier.", 2}
};

static const mp::OptionValueInfo values_tuneoutput_[] = {
  { "0", "None", 0},
  { "1", "Summarize each new best parameter set", 1},
  { "2", "Summarize each set tried (default)", 2},
  { "3", "Summary plus detailed solver output for each trial.", 3}
};

static const mp::OptionValueInfo values_varbranch[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Pseudo reduced - cost branching", 0},
  { "1", "Pseudo shadow - price branching", 1},
  { "2", "Maximum infeasibility branching", 2},
  { "3", "Strong branching.", 3}
};


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void GurobiBackend::InitCustomOptions() {

  set_option_header(
      fmt::format("{0} Optimizer Options for AMPL\n"
                  "-----------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``{1}_options``. For example::\n"
      "\n"
      "  ampl: option {1}_options 'opttol=1e-6';\n",
                  GetSolverName(), GetAMPLSolverName()).c_str());


  /// Options basis, sens are created internally if
  /// std features BASIS / SENSITIVITY are set.
  /// Change the help text
  AddToOptionDescription("alg:basis",
                         "Note that if you provide a valid starting extreme point, "
                         "either through primal/dual status, or through warmstart, "
                         "then Gurobi LP presolve will be disabled. For models where "
                         "presolve greatly reduces the problem size, "
                         "this might hurt performance."
      "\n\n"
      "For problems with "
      "integer variables or quadratic constraints, "
      "alg:basis=0 is assumed quietly unless mip:basis=1 or qcp:dual=1 "
      "is specified, respectively.");

  AddSolverOption("alg:cutoff cutoff",
    "If the optimal objective value is worse than cutoff, "
    "report \"objective cutoff\" and do not return a solution. "
    "Default: Infinity for minimizing, -Infinity for maximizing.",
    GRB_DBL_PAR_CUTOFF, MinusInfinity(), Infinity());

  AddToOptionDescription("alg:sens",
                         "For problems with integer variables or quadratic constraints, "
                         "alg:sens=0 is assumed quietly.");


  AddSolverOption("alg:method method lpmethod simplex",
    "Which algorithm to use for non-MIP problems or for the root node of MIP problems:\n"
    "\n.. value-table::\n", GRB_INT_PAR_METHOD, values_method, -1);


  AddSolverOption("alg:feasrelaxbigm feasrelaxbigm",
                  "Value of \"big-M\" sometimes used with constraints when doing "
                  "a feasibility relaxation.  Default = 1e6.",
                  GRB_DBL_PAR_FEASRELAXBIGM, 0.0, DBL_MAX);

  AddSolverOption("alg:feastol feastol",
                  "Primal feasibility tolerance (default 1e-6).",
                  GRB_DBL_PAR_FEASIBILITYTOL, 0.0, DBL_MAX);

  AddSolverOption("alg:numericfocus numericfocus",
                  "How much to try detecting and managing numerical issues:\n"
                  "\n"
      "| 0 - Automatic choice (default)\n"
      "| 1-3 - Increasing focus on more stable computations.",
                  GRB_INT_PAR_NUMERICFOCUS, 0, 3);

  AddSolverOption("bar:convtol barconvtol",
    "Tolerance on the relative difference between the primal and dual objectives "
    "for stopping the barrier algorithm "
    "(default 1e-8).", GRB_DBL_PAR_BARCONVTOL, 0.0, 1.0);


  AddSolverOption("bar:corr barcorrectors",
    "Limit on the number of central corrections done in each barrier iteration"
    "(default -1 = automatic choice).", GRB_INT_PAR_BARCORRECTORS, -1, GRB_MAXINT);

  AddSolverOption("bar:crossover crossover",
    "How to transform a barrier solution to a basic one:\n"
    "\n.. value-table::\n", GRB_INT_PAR_CROSSOVER, values_barcrossover, -1);

  AddSolverOption("bar:crossoverbasis crossoverbasis",
    "Strategy for initial basis construction during crossover:\n"
    "\n.. value-table::\n", GRB_INT_PAR_CROSSOVERBASIS, values_barcrossoverbasis, 0);

  AddSolverOption("bar:homog barhomogeneous",
    "Whether to use the homogeneous barrier algorithm (e.g., when method=2 is specified):\n"
    "\n.. value-table::\n"
    "The homogeneous barrier algorithm can detect infeasibility or unboundedness directly, "
    "without crossover, but is a bit slower than the nonhomogeneous barrier algorithm.",
    GRB_INT_PAR_BARHOMOGENEOUS, values_barhomogeneous, -1);


  AddSolverOption("bar:iterlim bariterlim",
    "Limit on the number of barrier iterations (default 1000).",
    GRB_INT_PAR_BARITERLIMIT, 0, GRB_MAXINT);

  AddSolverOption("bar:order barorder",
    "Ordering used to reduce fill in sparse-matrix factorizations "
    "during the barrier algorithm. Possible values:\n"
    "\n.. value-table::\n", GRB_INT_PAR_BARORDER, values_barorder, -1);

  AddSolverOption("bar:qcptol barqcptol",
    "Convergence tolerance on the relative difference between primal "
    "and dual objective values for barrier algorithms when solving problems "
    "with quadratic constraints (default 1e-6).", GRB_DBL_PAR_BARQCPCONVTOL,
    0.0, 1.0);


  /////////////////////// CUTS /////////////////////////
  AddSolverOption("cut:bqp bqpcuts",
    "Whether to enable Boolean Quadric Polytope cut generation:\n"
    "\n.. value-table::\n"
    "Overrides the \"cuts\" keyword.",
    GRB_INT_PAR_BQPCUTS, values_bqpcuts, -1);

  AddSolverOption("cut:clique cliquecuts",
    "Overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_CLIQUECUTS, PrmCutsMin, PrmCutsMax);

  AddSolverOption("cut:cover covercuts",
    "Overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_COVERCUTS, PrmCutsMin, PrmCutsMax);

  AddSolverOption("cut:agg cutagg cut:aggpasses",
    "Maximum number of constraint aggregation passes "
    "during cut generation (-1 = default = no limit); "
    "overrides \"cuts\".",
    GRB_INT_PAR_CUTAGGPASSES, -1, GRB_MAXINT);

  AddSolverOption("cut:passes cutpasses",
    "Maximum number of cutting-plane passes "
    "during root-cut generation; default = -1 ==> automatic choice.",
    GRB_INT_PAR_CUTPASSES, -1, GRB_MAXINT);

  AddSolverOption("cut:cuts cuts",
    "Global cut generation control, valid unless overridden "
    "by individual cut-type controls:"
    "\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_CUTS, values_cuts, -1);

  AddSolverOption("cut:flowcover flowcover",
    "Flowcover cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_FLOWCOVERCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:flowpath flowpath",
    "Overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_COVERCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:gomory gomory",
    "Maximum number of Gomory cut passes during cut generation "
        "(-1 = default = no limit); overrides \"cuts\".",
    GRB_INT_PAR_GOMORYPASSES, -1, GRB_MAXINT);
  AddSolverOption("cut:gubcover gubcover",
    "Overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_GUBCOVERCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:implied implied",
    "Implied cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_IMPLIEDCUTS, PrmCutsMin, PrmCutsMax);

  AddSolverOption("cut:infproof infproofcuts",
    "Whether to generate infeasibility proof cuts:"
    "\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_INFPROOFCUTS, values_infproofcuts, -1);

  AddSolverOption("cut:mipsep mipsep",
    "MIPsep cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_MIPSEPCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:mir mircuts",
    "MIR cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_MIRCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:modk modkcuts",
    "Mod-k cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_MODKCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:network networkcuts",
    "Network cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_NETWORKCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:relaxliftcuts relaxliftcuts",
    "Whether to enable relax-and-lift cut generation:\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_RELAXLIFTCUTS, values_cuts_upto2, -1);
  AddSolverOption("cut:rltcuts rltcuts",
    "Whether to enable generation of cuts by the Relaxation "
    "Linearization Technique (RLT):\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_RLTCUTS, values_cuts_upto2, -1);
  AddSolverOption("cut:submip submipcuts",
    "Sub-MIP cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_SUBMIPCUTS, PrmCutsMin, PrmCutsMax);
  AddSolverOption("cut:zerohalf zerohalfcuts",
    "Zero-half cuts: overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_ZEROHALFCUTS, PrmCutsMin, PrmCutsMax);




  AddSolverOption("lim:iter iterlim iterlimit",
    "Iteration limit (default: no limit).",
    GRB_DBL_PAR_ITERATIONLIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:minrelnodes minrelnodes",
    "Number of nodes for the Minimum Relaxation heuristic to "
    "explore at the MIP root node when a feasible solution has not "
    "been found by any other heuristic; default -1 ==> automatic choice.",
    GRB_INT_PAR_MINRELNODES, -1, GRB_MAXINT);

  AddSolverOption("lim:nodes nodelim nodelimit",
    "Maximum MIP nodes to explore (default: no limit).",
    GRB_DBL_PAR_NODELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:startnodes startnodelimit startnodes",
    "Limit on how many branch-and-bound nodes to explore when "
    "doing a partial MIP start:\n"
    "\n"
    "| -3 - Suppress MIP start processing\n"
    "| -2 - Only check full MIP starts for feasibility and "
            "ignore partial MIP starts\n"
    "| -1 - Use \"submipnodes\" (default)\n"
    "| >=0 - Specific node limit.",
    GRB_INT_PAR_STARTNODELIMIT, -3, GRB_MAXINT);

  AddSolverOption("lim:submipnodes submipnodes maxmipsub",
    "Limit on nodes explored by MIP-based heuristics, e.g., RINS. "
    "Default = 500.",
      GRB_INT_PAR_SUBMIPNODES, 0, GRB_MAXINT);

  AddSolverOption("lim:time timelim timelimit",
      "Limit on solve time (in seconds; default: no limit).",
      GRB_DBL_PAR_TIMELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:zeroobjnodes zeroobjnodes",
    "Number of nodes to explore in the zero objective heuristic. "
    "Note that this heuristic is only applied at the end of the "
    "MIP root, and only when no other root heuristic finds a "
    "feasible solution.\n"
    "\n"
    "This heuristic is quite expensive, and generally produces "
    "poor quality solutions. You should generally only use it if "
    "other means, including exploration of the tree with default "
    "settings, fail to produce a feasible solution.",
      GRB_INT_PAR_ZEROOBJNODES, -1, GRB_MAXINT);


  ////////////////////////// LP //////////////////////////

  AddSolverOption("lp:degenmoves degenmoves",
    "Limit on the number of degenerate simplex moves -- for use "
        "when too much time is taken after solving the initial root "
        "relaxation of a MIP problem and before cut generation or root "
        "heuristics have started.  Default -1 ==> automatic choice.",
    GRB_INT_PAR_DEGENMOVES, -1, GRB_MAXINT);

  AddSolverOption("lp:multprice_norm multprice_norm normadjust",
    "Choice of norm used in multiple pricing:\n"
    "\n"
      "|  -1 - Automatic choice (default)\n"
      "|  0, 1, 2, 3 - Specific choices:  hard to describe, "
        "but sometimes a specific choice will perform "
        "much better than the automatic choice.",
    GRB_INT_PAR_NORMADJUST, -1, 3);

  AddSolverOption("lp:perturb perturb",
    "Magnitude of simplex perturbation (when needed; default 2e-4).",
    GRB_DBL_PAR_PERTURBVALUE, 0.0, DBL_MAX);

  AddSolverOption("lp:pivtol pivtol markowitztol",
    "Markowitz pivot tolerance (default 7.8125e-3).",
    GRB_DBL_PAR_MARKOWITZTOL, 1e-4, 0.999);

  AddSolverOption("lp:pricing pricing",
    "Pricing strategy:\n"
    "\n.. value-table::",
    GRB_INT_PAR_SIMPLEXPRICING, values_pricing, -1);

  AddSolverOption("lp:quad quad",
    "Whether simplex should use quad-precision:\n"
    "\n.. value-table::",
    GRB_INT_PAR_QUAD, values_autonoyes_, -1);

  AddSolverOption("lp:sifting sifting",
    "Whether to use sifting within the dual simplex algorithm, "
    "which can be useful when there are many more variables than "
    "constraints:\n"
    "\n.. value-table::",
    GRB_INT_PAR_SIFTING, values_autonomodaggr_, -1);

  AddSolverOption("lp:siftmethod siftmethod",
    "Algorithm to use for sifting with the dual simplex method:\n"
    "\n.. value-table::",
    GRB_INT_PAR_SIFTMETHOD, values_siftmethod_, -1);


  ////////////////////////// MIP /////////////////////////

  AddSolverOption("mip:bestbndstop bestbndstop",
    "Stop once the best bound on the objective value "
    "is at least as good as this value.",
    GRB_DBL_PAR_BESTBDSTOP, MinusInfinity(), Infinity());

  AddSolverOption("mip:bestobjstop bestobjstop",
    "Stop after a feasible solution with objective value "
    "at least as good as this value has been found.",
    GRB_DBL_PAR_BESTOBJSTOP, MinusInfinity(), Infinity());

  AddSolverOption("mip:branchdir branchdir",
    "Which child node to explore first when branching:\n"
    "\n.. value-table::",
    GRB_INT_PAR_BRANCHDIR, values_branchdir, 0);


  AddSolverOption("mip:disconnected disconnected",
    "Whether to exploit independent MIP sub-models:"
    "\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_DISCONNECTED, values_disconnected, -1);

  AddStoredOption("mip:fixedmethod fixedmethod",
          "Value of \"method\" to use when seeking a basis for MIP problems "
          "when \"mip:basis=1\". Default: if \"method\" is 0 or 1 "
          "then \"method\" else 1.",
          storedOptions_.nFixedMethod_);


  AddSolverOption("mip:focus mipfocus",
    "MIP solution strategy:\n" "\n.. value-table::\n",
    GRB_INT_PAR_MIPFOCUS, values_mipfocus, 0);


  AddSolverOption("mip:gap mipgap",
    "Max. relative MIP optimality gap (default 1e-4).",
    GRB_DBL_PAR_MIPGAP, 1e-4, DBL_MAX);

  AddSolverOption("mip:gapabs mipgapabs",
    "Max. absolute MIP optimality gap (default 1e-10).",
    GRB_DBL_PAR_MIPGAPABS, 1e-10, DBL_MAX);

  AddSolverOption("mip:heurfrac heurfrac",
    "Fraction of time to spend in MIP heuristics (default 0.05).",
    GRB_DBL_PAR_HEURISTICS, 0.05, 1.0);

  AddSolverOption("alg:iismethod iismethod",
    "Which method to use when finding an IIS (irreducible infeasible "
    "set of constraints, including variable bounds):\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_IISMETHOD, values_iismethod, -1);


  AddSolverOption("mip:improvegap improvegap",
    "Optimality gap below which the MIP solver switches from "
        "trying to improve the best bound to trying to find better "
        "feasible solutions (default 0).",
    GRB_DBL_PAR_IMPROVESTARTGAP, 0.0, DBL_MAX);
  AddSolverOption("mip:improvetime improvetime",
    "Execution seconds after which the MIP solver switches from "
        "trying to improve the best bound to trying to find better "
        "feasible solutions (default Infinity).",
    GRB_DBL_PAR_IMPROVESTARTTIME, 0.0, DBL_MAX);
  AddSolverOption("mip:impstartnodes impstartnodes",
                  "Number of MIP nodes after which the solution strategy "
                      "will change from improving the best bound to finding better "
                      "feasible solutions (default Infinity).",
    GRB_DBL_PAR_IMPROVESTARTNODES, 0.0, DBL_MAX);

  AddSolverOption("mip:intfeastol intfeastol",
    "Feasibility tolerance for integer variables (default 1e-05).",
    GRB_DBL_PAR_INTFEASTOL, 0.0, DBL_MAX);


  AddSolverOption("mip:intfocus integralityfocus intfocus",
                  "Setting this parameter to 1 requests the solver to work "
                  "harder at finding solutions that are still (nearly) feasible "
                  "when all integer variables are rounded to exact integral "
                  "values to avoid numerical issues such as trickle flow:\n"
                  "\n.. value-table::\n",
                  GRB_INT_PAR_INTEGRALITYFOCUS, values_01_noyes_0default_, 0);
  AddToOptionDescription("mip:round",
                         "With Gurobi, for problems with numerical issues such as trickle flow, "
                         "option \"mip:intfocus\" can be more reliable.");


  AddToOptionDescription("mip:lazy",
                         "For Gurobi, lazy/user constraints are indicated with .lazy values of "
    "-1, "
    "1, 2, "
    "or 3 and are ignored until a solution feasible to the "
    "remaining constraints is found.  What happens next depends "
    "on the value of .lazy:\n"
    "\n"
    "|  -1 - Treat the constraint as a user cut; the "
            "constraint must be redundant with respect to the "
            "rest of the model, although it can cut off LP "
            "solutions;\n"
    "|  1  - The constraint may still be ignored if another "
            "lazy constraint cuts off the current solution;\n"
    "|  2  - The constraint will henceforth be enforced if it "
            "is violated by the current solution;\n"
    "|  3  - The constraint will henceforth be enforced."
  );



  AddSolverOption("mip:nodemethod nodemethod",
    "Algorithm used to solve relaxed MIP node problems:\n"
    "\n.. value-table::\n", GRB_INT_PAR_NODEMETHOD, values_nodemethod, -1);

  AddSolverOption("mip:norelheurtime norelheurtime",
    "Limits the amount of time (in seconds) spent in the NoRel heuristic; "
    "see the description of \"norelheurwork\" for details.  This "
    "parameter will introduce nondeterminism; use \"norelheurwork\" "
    "for deterministic results.  Default 0.",
    GRB_DBL_PAR_NORELHEURTIME, 0.0, DBL_MAX);

  AddSolverOption("mip:norelheurwork norelheurwork",
    "Limits the amount of work spent in the NoRel heuristic. "
    "This heuristic searches for high-quality feasible solutions "
    "before solving the root relaxation.  The work metrix is hard "
    "to define precisely, as it depends on the machine.  Default 0.",
    GRB_DBL_PAR_NORELHEURWORK, 0.0, DBL_MAX);



  AddSolverOption("mip:opttol opttol optimalitytolerance",
      "Dual feasibility tolerance.",
      GRB_DBL_PAR_OPTIMALITYTOL, 1e-9, 1e-2);

  AddSolverOption("mip:partition partitionplace",
      "Whether and how to use the .partition suffix on variables "
    "in the partition heuristic for MIP problems: sum of\n"
    "\n"
    "|      1 ==> When the branch-and-cut search ends\n"
    "|      2 ==> At nodes in the branch-and-cut search\n"
    "|      4 ==> After the root-cut loop\n"
    "|      8 ==> At the start of the root-cut loop\n"
    "|     16 ==> Before solving the root relaxation.\n"
    "\n"
    "Default = 15.  Values of .parition determine how variables "
    "participate in the partition heuristic.  Variables with\n"
    "\n"
    "|     .partition = -1 are always held fixed;\n"
    "|     .partition = 0 can vary in all sub-MIP models;\n"
    "|     .partition > 0 can vary only in in that sub-MIP model.\n"
    "\n"
    "The partition heuristic is only run when partition_place is "
    "between 1 and 31 and some variables have suitable .partition "
    "suffix values.",
      GRB_INT_PAR_PARTITIONPLACE, 0, 31);

  AddSolverOption("mip:pumppasses pumppasses",
    "Number of feasibility-pump passes to do after the MIP root "
    "when no other root heuristoc found a feasible solution. "
    "Default -1 = automatic choice.",
    GRB_INT_PAR_PUMPPASSES, -1, GRB_MAXINT);

  AddSolverOption("mip:rins rins",
    "How often to apply the RINS heuristic for MIP problems:\n"
    "\n"
    "| -1  - Automatic choice (default)\n"
    "| 0   - never\n"
    "| n > 0  - every n-th node.",
    GRB_INT_PAR_RINS, -1, GRB_MAXINT);

  AddStoredOption("mip:start mipstart intstart",
    "Whether to use initial guesses in problems with "
    "integer variables:\n"   "\n.. value-table::\n",
    storedOptions_.nMIPStart_, values_mipstart_);
  AddToOptionDescription("alg:start",
                         "Note that for LP, \"alg:basis\" is usually more efficient.\n"
                         "For Gurobi, "
                         "MIP-specific options can be tuned via \"mip:start\".");

  AddSolverOption("mip:symmetry symmetry",
    "MIP symmetry detection:\n"   "\n.. value-table::\n",
    GRB_INT_PAR_SYMMETRY, values_autonoconsaggr_, -1);

  AddSolverOption("mip:varbranch varbranch",
    "MIP branch variable selection strategy:\n"
    "\n.. value-table::\n", GRB_INT_PAR_VARBRANCH, values_varbranch, -1);


  AddSolverOption("miqcp:method miqcpmethod",
    "Method for solving mixed-integer quadratically constrained "
    "(MIQCP) problems:\n"  "\n.. value-table::\n",
    GRB_INT_PAR_MIQCPMETHOD, values_miqcpmethod, -1);


  /// Option "multiobj" is created internally if
  /// std feature MULTIOBJ is set.
  /// Change the help text
  ReplaceOptionDescription("obj:multi",
                           "0*/1: Whether to do multi-objective optimization.\n"
                           "When obj:multi = 1 and several objectives are present, suffixes "
                           ".objpriority, .objweight, .objreltol, and .objabstol on the "
                           "objectives are relevant.  Objectives with greater .objpriority "
                           "values (integer values) have higher priority.  Objectives with "
                           "the same .objpriority are weighted by .objweight.  Objectives "
                           "with positive .objabstol or .objreltol are allowed to be "
                           "degraded by lower priority objectives by amounts not exceeding "
                           "the .objabstol (absolute) and .objreltol (relative) limits. "
                           "The objectives must all be linear.  Objective-specific "
                           "convergence tolerances and method values may be assigned via "
                           "keywords of the form obj_n_<name>, such as obj_1_method for the "
                           "first objective.");

  AddSolverOption("obj:multiobjmethod multiobjmethod",
    "Choice of optimization algorithm for lower-priority objectives:\n"
    "\n.. value-table::\n"
    "The method keyword determines the algorithm to use for the highest "
    "priority objective.",
    GRB_INT_PAR_MULTIOBJMETHOD, values_multiobjmethod, -1);

  AddSolverOption("obj:multiobjpre multiobjpre",
    "How to apply Gurobi's presolve when doing multi-objective optimization:\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_MULTIOBJPRE, values_multiobjpre, -1);

  AddSolverOption("obj:scale objscale",
    "How to scale the objective:\n"
    "\n"
      "| 0 - Automatic choice (default)\n"
      "| -1..0 - Divide by max abs. coefficient "
             "raised to this power\n"
      "| >0 - Divide by this value.",
    GRB_DBL_PAR_OBJSCALE, -1.0, DBL_MAX);

  /// Wildcard options
  AddIntOption("obj:*:method obj_*_method", "Method for objective with index *",
            &GurobiBackend::GrbGetObjIntParam,
            &GurobiBackend::GrbSetObjIntParam);
  AddIntOption("obj:*:priority obj_*_priority", "Priority for objective with index *",
            &GurobiBackend::GrbGetObjIntParam,
            &GurobiBackend::GrbSetObjIntParam);
  AddDblOption("obj:*:weight obj_*_weight", "Weight for objective with index *",
            &GurobiBackend::GrbGetObjDblParam,
            &GurobiBackend::GrbSetObjDblParam);
  AddDblOption("obj:*:reltol obj_*_reltol", "Relative tolerance for objective with index *",
            &GurobiBackend::GrbGetObjDblParam,
            &GurobiBackend::GrbSetObjDblParam);
  AddDblOption("obj:*:abstol obj_*_abstol", "Absolute tolerance for objective with index *. "
               "Can only be applied on a multi-objective problem with obj:multi=1",
            &GurobiBackend::GrbGetObjDblParam,
            &GurobiBackend::GrbSetObjDblParam);


  AddSolverOption("pre:aggfill aggfill", "Amount of fill allowed during aggregation in presolve"
    "(default -1).", GRB_INT_PAR_AGGFILL, -1, GRB_MAXINT);


  AddSolverOption("pre:aggregate aggregate", "0/1*: whether to use aggregation in presolve."
    "Setting it to 0 can sometimes reduce numerical errors.", GRB_INT_PAR_AGGREGATE, 0, 1);

  AddSolverOption("pre:deprow predeprow",
    "Whether Gurobi's presolve should remove linearly dependent "
    "constraint-matrix rows:\n" "\n.. value-table::\n",
    GRB_INT_PAR_PREDEPROW, values_predeprow, -1);

  AddSolverOption("pre:dual predual",
    "Whether Gurobi's presolve should form the dual of a "
    "continuous model:\n" "\n.. value-table::\n",
    GRB_INT_PAR_PREDUAL, values_predual, -1);

  AddSolverOption("pre:dualreductions dualreductions",
    "Whether Gurobi's presolve should use dual reductions, which "
        "may be useful on a well-posed problem but can prevent "
        "distinguishing whether a problem is infeasible or unbounded:\n"
        "\n.. value-table::\n",
    GRB_INT_PAR_DUALREDUCTIONS, values_01_noyes_1default_, 1);

  AddSolverOption("pre:miqcpform premiqcpform",
    "For mixed-integer quadratically constrained (MIQCP) problems, "
    "how Gurobi should transform quadratic constraints:\n"
        "\n.. value-table::\n\n"
    "Choices 0 and 1 work with general quadratic constraints. "
    "Choices 1 and 2 only work with constraints of suitable forms.",
    GRB_INT_PAR_PREMIQCPFORM, values_premiqcpform, -1);

  AddSolverOption("pre:passes prepasses",
    "Limit on the number of Gurobi presolve passes:\n"
    "\n"
    "| -1 - Automatic choice (default)\n"
    "| n>=0 - At most n passes.",
    GRB_INT_PAR_PREPASSES, -1, GRB_MAXINT);

  AddSolverOption("pre:qlinearize preqlinearize preqlin",
    "How Gurobi's presolve should treat quadratic problems:\n"
        "\n.. value-table::\n",
    GRB_INT_PAR_PREQLINEARIZE, values_preqlinearize, -1);

  AddSolverOption("pre:scale scale",
    "Whether to scale the problem:\n"
        "\n.. value-table::\n"
    "Scaling typically reduces solution times, but it may lead "
    "to larger constraint violations in the original, unscaled "
    "model. Choosing a different scaling option can sometimes "
    "improve performance for particularly numerically difficult "
    "models.",
    GRB_INT_PAR_SCALEFLAG, values_prescale, -1);

  AddSolverOption("pre:solve presolve",
    "Whether to use Gurobi's presolve:\n"
        "\n.. value-table::\n",
    GRB_INT_PAR_PRESOLVE, values_autonoconsaggr_, -1);

  AddSolverOption("pre:sos1bigm presos1bigm",
    "Big-M for converting SOS1 constraints to binary form:\n"
    "\n"
    "| -1 - Automatic choice (default)\n"
    "|  0  - No conversion\n"
    "\n"
    "Large values (e.g., 1e8) may cause numeric trouble.",
    GRB_DBL_PAR_PRESOS1BIGM, -1.0, 1e10);

  AddSolverOption("pre:sos2bigm presos2bigm",
    "Big-M for converting SOS2 constraints to binary form:\n"
    "\n"
    "| -1 - Automatic choice (default)\n"
    "|  0  - No conversion\n"
    "\n"
    "Large values (e.g., 1e8) may cause numeric trouble.",
    GRB_DBL_PAR_PRESOS2BIGM, -1.0, 1e10);

  AddSolverOption("pre:sparsify presparsify",
    "Whether Gurobi's presolve should use its \"sparsify reduction\", "
    "which sometimes gives significant problem-size reductions:\n"
        "\n.. value-table::\n",
    GRB_INT_PAR_PRESPARSIFY, values_autonoyes_, -1);



  ///////////////////////// Q(C)P //////////////////////////////
  AddSolverOption("qp:nonconvex nonconvex",
                  "How to handle non-convex quadratic objectives "
                  "and constraints:\n"  "\n.. value-table::\n",
    GRB_INT_PAR_NONCONVEX, values_nonconvex, -1);

  AddSolverOption("qcp:dual qcpdual",
                  "Whether to compute dual variables when the problem "
                  "has quadratic constraints (which can be expensive):\n"
                  "\n.. value-table::\n",
    GRB_INT_PAR_QCPDUAL, values_01_noyes_0default_, 0);

  AddSolverOption("qp:psdtol psdtol",
                  "Maximum diagonal perturbation to correct indefiniteness "
                  "in quadratic objectives (default 1e-6).",
    GRB_DBL_PAR_PSDTOL, 0.0, DBL_MAX);


  /// Solution pool parameters
  /// Rely on MP's built-in options solutionstub or countsolutions
  AddSolverOption("sol:poolgap ams_eps poolgap",
      "Relative tolerance for reporting alternate MIP solutions "
      "		(default: Infinity, no limit).",
      GRB_DBL_PAR_POOLGAP, 0.0, DBL_MAX);
  AddSolverOption("sol:poolgapabs ams_epsabs poolgapabs",
      "Absolute tolerance for reporting alternate MIP solutions "
      "		(default: Infinity, no limit).",
      GRB_DBL_PAR_POOLGAPABS, 0.0, DBL_MAX);
  AddSolverOption("sol:poollimit ams_limit poollimit solnlimit",
      "Limit on the number of alternate MIP solutions written. Default: 10.",
      GRB_INT_PAR_POOLSOLUTIONS, 1, 2000000000);
  AddStoredOption("sol:poolmode ams_mode poolmode",
      "Search mode for MIP solutions when sol:stub/sol:count are specified "
                        "to request finding several alternative solutions:\n",
          storedOptions_.nPoolMode_, values_pool_mode);
  AddOptionSynonyms_Inline_Front("ams_stub", "sol:stub");

  /// Option "solutionstub" is created internally if
  /// std feature MULTISOL is set.
  /// Change the help text
  ReplaceOptionDescription("sol:stub",
                           "Stub for alternative MIP solutions, written to files with "
                        "names obtained by appending \"1.sol\", \"2.sol\", etc., to "
                        "<solutionstub>.  The number of such files written is affected "
                        "by the keywords poolgap, poolgapabs, poollimit, and poolmode. "
                        "The number of alternative MIP solution files written is "
                        "returned in suffix .nsol on the problem.");
  ReplaceOptionDescription("sol:count",
                           "0*/1: Whether to count the number of solutions "
                           "and return it in the ``.nsol`` problem suffix. "
                           "The number and kind of solutions are controlled by the "
                           "sol:pool... parameters. Value 1 implied by sol:stub.");


  AddStoredOption("tech:cloudid cloudid",
      "Use Gurobi Instant Cloud with this \"accessID\".",
          storedOptions_.cloudid_);
  AddStoredOption("tech:cloudkey cloudkey",
      "Use Gurobi Instant Cloud with this \"secretKey\". "
      "Both cloudid and cloudkey are required.",
          storedOptions_.cloudkey_);
  AddStoredOption("tech:cloudpool cloudpool",
      "Optional \"machine pool\" to use with Gurobi Instant Cloud.",
          storedOptions_.cloudpool_);
  AddStoredOption("tech:cloudpriority cloudpriority",
      "Priority of Cloud job, an integer >= -100 and <= 100. "
      "Default 0.  Jobs with priority 100 run immediately -- use "
      "caution when specifying this value.",
          storedOptions_.cloudpriority_, -100, 100);


  AddSolverOption("tech:concurrentmip concurrentmip",
      "How many independent MIP solves to allow at once when multiple "
      "threads are available. Optimization terminates when the first solve "
      "completes. The available threads are divided as "
      "evenly as possible among the concurrent solves.  Default = 1.\n"
      "\n"
      "See also \"tech:distmip\", \"tech:pooljobs\".",
          GRB_INT_PAR_CONCURRENTMIP, 1, GRB_MAXINT);

  AddSolverOption("tech:distmip pool_distmip distmip",
      "Enables distributed MIP. A value of n causes the MIP solver "
      "to divide the work of solving a MIP model among n machines. "
      "Use the \"tech:server\" parameter to indicate the name of the "
      "cluster where you would like your distributed MIP job to run "
      "(or use \"tech:workerpool\" if your client machine will act as manager "
      "and you just need a pool of workers). Default = 0.\n"
      "\n"
      "See also \"tech:concurrentmip\", \"tech:pooljobs\".",
          GRB_INT_PAR_DISTRIBUTEDMIPJOBS, 0, GRB_MAXINT);

  AddSolverOption("tech:pooljobs pool_jobs pooljobs",
      "Enables distributed concurrent optimization, which can be used "
      "to solve LP or MIP models on multiple machines. A value of n "
      "causes the solver to create n independent models, using different "
      "parameter settings for each. Each of these models is sent to a "
      "distributed worker for processing. Optimization terminates when "
      "the first solve completes. Use the \"tech:server\" parameter to "
      "indicate the name of the cluster where you would like your "
      "distributed concurrent job to run (or use \"tech:workerpool\" if your "
      "client machine will act as manager and you just need a pool of "
      "workers). Default = 0.\n"
      "\n"
      "See also \"tech:concurrentmip\", \"tech:distmip\".",
          GRB_INT_PAR_CONCURRENTJOBS, 0, GRB_MAXINT);


  AddSolverOption("tech:logfreq logfreq outfreq",
      "Interval in seconds between log lines (default 5).",
    GRB_INT_PAR_DISPLAYINTERVAL, 1, GRB_MAXINT);

  AddSolverOption("tech:logfile logfile",
      "Log file name.",
      GRB_STR_PAR_LOGFILE);

  AddSolverOption("tech:nodefiledir nodefiledir",
      "Directory where MIP tree nodes are written after memory "
      "for them exceeds \"nodefilestart\"; default \".\"",
    GRB_STR_PAR_NODEFILEDIR);

  AddSolverOption("tech:nodefilestart nodefilestart",
      "Gigabytes of memory to use for MIP tree nodes; "
      "default = Infinity (no limit, i.e., no node files written).",
    GRB_DBL_PAR_NODEFILESTART, 0.0, DBL_MAX);



  AddSolverOption("tech:outlev outlev",
      "0*/1: Whether to write gurobi log lines (chatter) to stdout and to file.",
    GRB_INT_PAR_OUTPUTFLAG, 0, 1);

  AddSolverOption("tech:param param",
                  "General way to specify values of both documented and "
    "undocumented Gurobi parameters; value should be a quoted "
    "string (delimited by ' or \") containing a parameter name, a "
    "space, and the value to be assigned to the parameter.  Can "
    "appear more than once.  Cannot be used to query current "
    "parameter values.",
                  "Dummy");

  AddStoredOption("tech:param:read param:read paramfile",
      "Name of Gurobi parameter file (surrounded by 'single' or "
      "\"double\" quotes if the name contains blanks). "
      "The suffix on a parameter file should be .prm, optionally followed "
      "by .zip, .gz, .bz2, or .7z.\n"
      "\n"
      "Lines that start with # are ignored.  Otherwise, each nonempty "
      "line should contain a name and a value, separated by a space.",
          storedOptions_.paramRead_);
  AddStoredOption("tech:param:write param:write",
      "Name of Gurobi parameter file (surrounded by 'single' or \"double\" quotes if the "
      "name contains blanks) to be written.",
          storedOptions_.paramWrite_);

  AddSolverOption("tech:resultfile resultfile",
      "Name of a file of extra information written after "
      "completion of optimization.  The name's suffix determines what "
      "is written:\n"
      "\n"
      "| .sol - Solution vector\n"
      "| .bas - Simplex basis\n"
      "| .mst - Integer variable solution vector\n"
      "| .ilp - IIS for an infeasible model\n"
      "| .mps, .rew, .lp, or .rlp - To capture the original model.\n"
      "\n"
      "The file suffix may optionally be followed by .gz, .bz2, or .7z, "
      "which produces a compressed result.",
      GRB_STR_PAR_RESULTFILE);


  AddSolverOption("tech:seed seed",
      "Random number seed (default 0), affecting perturbations that "
      "may influence the solution path.",
      GRB_INT_PAR_SEED, 0, GRB_MAXINT);


  AddStoredOption("tech:server server servers",
      "Comma-separated list of Gurobi Compute Servers, specified "
      "either by name or by IP address.  Default: run Gurobi locally "
      "(i.e., do not use a remote Gurobi server).",
          storedOptions_.servers_);
  AddOptionSynonyms_OutOfLine("tech:server_lic serverlic server_lic", "tech:optionfile");
  AddStoredOption("tech:server_group server_group",
      "Name of Compute Server Group, if any.",
          storedOptions_.server_group_);
  AddStoredOption("tech:server_router server_router",
      "Name or IP address of router for Compute Server, if any.",
          storedOptions_.server_router_);
  AddStoredOption("tech:server_insecure server_insecure",
      "Whether to user \"insecure mode\" with the Gurobi Compute "
      "Server.  Should be left at default value (0) unless an "
      "administrator specifies another value.",
          storedOptions_.server_insecure_, -GRB_MAXINT, GRB_MAXINT);
  AddStoredOption("tech:server_password server_password",
      "Password (if needed) for specified Gurobi Compute Server(s).",
          storedOptions_.server_password_);
  AddStoredOption("tech:server_priority server_priority",
      "Priority for Gurobi Compute Server(s).  Default = 0. "
      "Highest priority = 100.",
          storedOptions_.server_priority_, -100, 100);
  AddStoredOption("tech:server_timeout server_timeout",
      "Report job as rejected by Gurobi Compute Server if the "
      "job is not started within server_timeout seconds. "
      "Default = 10.",
          storedOptions_.server_timeout_, 1.0, DBL_MAX);


  AddSolverOption("tech:threads threads",
      "How many threads to use when using the barrier algorithm "
      "or solving MIP problems; default 0 ==> automatic choice.",
      GRB_INT_PAR_THREADS, 0, GRB_MAXINT);

  AddStoredOption("tech:tunebase tunebase",
                  "Base name for results of running Gurobi's search for best "
                  "parameter settings.  The search is run only when tunebase "
                  "is specified.  Results are written to files with names derived "
                  "from tunebase by appending \".prm\" if \".prm\" does not occur in "
                  "tunebase and inserting _1_, _2_, ... (for the first, second, "
                  "... set of parameter settings) before the right-most \".prm\". "
                  "The file with _1_ inserted is the best set and the solve "
                  "results returned are for this set.  In a subsequent \"solve;\", "
                  "you can use \"tech:param:read=...\" to apply the settings in results "
                  "file ... .",
                  storedOptions_.tunebase_);
  AddSolverOption("tech:tunejobs pool_tunejobs tunejobs",
      "Enables distributed parallel tuning, which can significantly "
      "increase the performance of the tuning tool. A value of n "
      "causes the tuning tool to distribute tuning work among n "
      "parallel jobs. These jobs are distributed among a set of "
      "machines. Use the \"tech:workerpool\" parameter to provide a distributed "
      "worker cluster. Default = 0.\n"
      "\n"
      "Note that distributed tuning is most effective when the worker "
      "machines have similar performance.",
          GRB_INT_PAR_TUNEJOBS, 0, GRB_MAXINT);
  AddSolverOption("tech:tuneoutput tuneoutput",
      "Amount of tuning output when tunebase is specified:\n"
                        "\n.. value-table::\n",
          GRB_INT_PAR_TUNEOUTPUT, values_tuneoutput_, 2);
  AddSolverOption("tech:tuneresults tuneresults",
      "Limit on the number of tuning result files to write "
      "when tunerbase is specified.  The default (-1) is to write "
      "results for all parameter sets on the efficient frontier.",
          GRB_INT_PAR_TUNERESULTS, -1, GRB_MAXINT);
  AddSolverOption("tech:tunetimelim tunetimelim lim:tunetime",
      "Time limit (in seconds) on tuning when tunebase "
      "is specified.  Default -1 ==> automatic choice of time limit.",
          GRB_DBL_PAR_TUNETIMELIMIT, -1.0, DBL_MAX);
  AddSolverOption("tech:tunetrials tunetrials",
      "Number of trials for each parameter set when tunebase "
      "is specified, each with a different random seed value. "
      "Default = 3.",
          GRB_INT_PAR_TUNETRIALS, 1, GRB_MAXINT);



  AddSolverOption("tech:workerpool pool_servers",
      "When using a distributed algorithm (distributed MIP, "
      "distributed concurrent, or distributed tuning), this "
      "parameter allows you to specify a Remote Services "
      "cluster that will provide distributed workers. You "
      "should also specify the access password for that "
      "cluster, if there is one, in the \"workerpassword\" "
      "parameter. Note that you don't need to set either "
      "of these parameters if your job is running on a "
      "Compute Server node and you want to use the same "
      "cluster for the distributed workers.\n"
      "\n"
      "You can provide a comma-separated list of machines "
      "for added robustness.",
      GRB_STR_PAR_WORKERPOOL);

  AddSolverOption("tech:workerpassword pool_password",
      "Password for the worker pool (if needed).",
      GRB_STR_PAR_WORKERPASSWORD);



  AddStoredOption("tech:exportfile writeprob writemodel exportfile",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);
}


void GurobiBackend::GrbSetObjIntParam(const SolverOption& opt, int val) {
  objnparam_int_.push_back( { {opt.wc_tail(), opt.wc_keybody_last()}, val } );
}
void GurobiBackend::GrbSetObjDblParam(const SolverOption& opt, double val) {
  objnparam_dbl_.push_back( { {opt.wc_tail(), opt.wc_keybody_last()}, val } );
}
int GurobiBackend::GrbGetObjIntParam(const SolverOption& opt) const {
  auto it = std::find_if(objnparam_int_.rbegin(), objnparam_int_.rend(),
                         [&](const ObjNParam<int>& prm){
    return prm.first == std::make_pair(opt.wc_tail(), opt.wc_keybody_last());
  });
  if (objnparam_int_.rend()==it)
    throw std::runtime_error("Failed to find recorded option " +
                             opt.wc_key_last__std_form());
  return it->second;
}
double GurobiBackend::GrbGetObjDblParam(const SolverOption& opt) const {
  auto it = std::find_if(objnparam_dbl_.rbegin(), objnparam_dbl_.rend(),
                         [&](const ObjNParam<int>& prm){
    return prm.first == std::make_pair(opt.wc_tail(), opt.wc_keybody_last());
  });
  if (objnparam_dbl_.rend()==it)
    throw std::runtime_error("Failed to find recorded option " +
                             opt.wc_key_last__std_form());
  return it->second;
}

/// What to do on certain "obj:*:..." option
/// The bool=true means get multiobj env
/// int is the obj number
static std::tuple<bool, int, const char*>
GrbGetObjParamAction(const GurobiBackend::ObjNParamKey& key) {
  int n;
  try {
    n = std::stoi(key.second)-1;  // subtract 1 for 0-based indexing
  }  catch (...) {
    throw std::runtime_error("Could not parse index '" + key.second +
                             "' of option 'obj:" + key.second + key.first + "'");
  }
  if (":method"==key.first)
    return {true, n, GRB_INT_PAR_METHOD};
  if (":priority"==key.first)
    return {false, n, GRB_INT_ATTR_OBJNPRIORITY};
  if (":weight"==key.first)
    return {false, n, GRB_DBL_ATTR_OBJNWEIGHT};
  if (":abstol"==key.first)
    return {false, n, GRB_DBL_ATTR_OBJNABSTOL};
  if (":reltol"==key.first)
    return {false, n, GRB_DBL_ATTR_OBJNRELTOL};
  throw std::runtime_error(
        "Unknown wildcard option 'obj:" + key.second + key.first + "'");
  return { false, -1, "" };
}

/// env() is only for error reporting
static void GrbDoSetObjParam(
    const GurobiBackend::ObjNParam<int>& prm,
    GRBmodel* model, GRBenv* env_) {
  auto action = GrbGetObjParamAction(prm.first);
  auto iobj = std::get<1>(action);
  auto prm_attr = std::get<2>(action);
  auto env = [env_]() { return env_; };            // for GRB_CALL
  if (std::get<0>(action)) {
    auto envN = GRBgetmultiobjenv(model, iobj);
    GRB_CALL( GRBsetintparam(envN, prm_attr, prm.second) );
  } else {
    GRB_CALL( GRBsetintparam(GRBgetenv(model),
                             GRB_INT_PAR_OBJNUMBER, iobj) );
    GRB_CALL( GRBsetintattr(model, prm_attr, prm.second) );
  }
}
/// env() is only for error reporting
static void GrbDoSetObjParam(
    const GurobiBackend::ObjNParam<double>& prm,
    GRBmodel* model, GRBenv* env_) {
  auto action = GrbGetObjParamAction(prm.first);
  auto iobj = std::get<1>(action);
  auto prm_attr = std::get<2>(action);
  auto env = [env_]() { return env_; };            // for GRB_CALL
  if (std::get<0>(action)) {
    auto envN = GRBgetmultiobjenv(model, iobj);
    GRB_CALL( GRBsetdblparam(envN, prm_attr, prm.second) );
  } else {
    GRB_CALL( GRBsetintparam(GRBgetenv(model),
                             GRB_INT_PAR_OBJNUMBER, iobj) );
    GRB_CALL( GRBsetdblattr(model, prm_attr, prm.second) );
  }
}

template <class T>
static void DoPlayGrbObjNParams(
    const std::vector< GurobiBackend::ObjNParam<T> >& objnp,
    GRBmodel* model, GRBenv* env) {
  for (const auto& p: objnp)
    GrbDoSetObjParam(p, model, env);
}

void GurobiBackend::GrbPlayObjNParams() {
  int nobj;
  GRB_CALL( GRBgetintattr(model(), GRB_INT_ATTR_NUMOBJ, &nobj) );
  if (debug_mode())
    printf("Number of objectives: %d\n", nobj);
  DoPlayGrbObjNParams(objnparam_int_, model(), env());
  DoPlayGrbObjNParams(objnparam_dbl_, model(), env());
}


} // namespace mp


AMPLS_MP_Solver* AMPLSOpenGurobi(const char* slv_opt) {
  AMPLS_MP_Solver* slv =
      AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::GurobiBackend()},
        slv_opt);
  return slv;
}

void AMPLSCloseGurobi(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

GRBmodel* GetGRBmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::GurobiBackend*>(AMPLSGetBackend(slv))->model();
}

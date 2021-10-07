#include <vector>
#include <climits>
#include <cfloat>

#include "gurobibackend.h"

#define GRB_CALL( call ) do { if (int e=call) MP_RAISE( \
    fmt::format("  Call failed: '{}' with code {}, message: {}", #call, \
        e, GRBgeterrormsg(env_) ) \
  ); } while (0)

namespace {

bool InterruptGurobi(void *model) {
  GRBterminate( static_cast<GRBmodel*>(model) );
  return true;
}

}  // namespace

namespace mp {

const char* GurobiBackend::GetSolverInvocationName() { return "gurobidirect"; }
const char* GurobiBackend::GetBackendName() { return "GurobiBackend"; }

std::string GurobiBackend::GetSolverVersion() {
  int a,b,c;
  GRBversion(&a, &b, &c);
  return fmt::format("{}.{}.{}", a, b, c);
}

GurobiBackend::GurobiBackend() {}

GurobiBackend::~GurobiBackend() {
  CloseGurobi();
}

void GurobiBackend::CloseGurobi() {
  /* Free the fixed model */
  if (model_ != model_fixed_) {
    assert(model_);
    assert(model_fixed_);
    GRBfreemodel(model_fixed_);
  }
  model_fixed_ = nullptr;

  /* Free model */
  if (model_) {
    GRBfreemodel(model_);
    model_ = nullptr;
  }

  /* Free environment */
  if (env_) {
    GRBfreeenv(env_);
    env_ = nullptr;
  }
}

void GurobiBackend::InitOptionParsing() {
  OpenGurobi();
}

void GurobiBackend::OpenGurobi() {
  GRB_CALL( GRBloadenv(&env_, NULL) );
  OpenGurobiModel();
}

void GurobiBackend::OpenGurobiModel() {
  /* Set default parameters */
  GRBsetintparam(env_, GRB_INT_PAR_OUTPUTFLAG, 0);
  /* Create an empty model */
  GRB_CALL( GRBnewmodel(env_, &model_, "amplgurobi", 0,
                        NULL, NULL, NULL, NULL, NULL) );
  assert(model_);
  /* Init fixed model */
  model_fixed_ = model_;
}

void GurobiBackend::FinishOptionParsing() {
  if (cloudid().size() && cloudkey().size()) {
    OpenGurobiCloud();
  }
}

void GurobiBackend::OpenGurobiCloud() {
  assert(cloudid().size() && cloudkey().size());
  auto logf = GrbGetStrParam(GRB_STR_PAR_LOGFILE);
  if (env_) {
    CloseGurobi();
  }
  if (int i = GRBloadcloudenv(&env_, logf.c_str(),
                              cloudid().c_str(), cloudkey().c_str(),
                              cloudpool().c_str(), cloudpriority()
                              )) {
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
              "Surprise return {} from GRBloadcloudenv().", i));
    }
  }
  OpenGurobiModel();
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

int GurobiBackend::NumberOfConstraints() const {
  return GrbGetIntAttr(GRB_INT_ATTR_NUMCONSTRS);
}

int GurobiBackend::NumberOfVariables() const {
  return GrbGetIntAttr(GRB_INT_ATTR_NUMVARS);
}

int GurobiBackend::NumberOfObjectives() const {
  return GrbGetIntAttr(GRB_INT_ATTR_NUMOBJ);
}

int GurobiBackend::ModelSense() const {
  return GrbGetIntAttr(GRB_INT_ATTR_MODELSENSE);
}

ArrayRef<double> GurobiBackend::PrimalSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_X, NumberOfVariables());
}

ArrayRef<double> GurobiBackend::DualSolution() {
  return MakeDualsFromLPAndQCPDuals(
        GurobiDualSolution_LP(), GurobiDualSolution_QCP());
}

std::vector<double> GurobiBackend::GurobiDualSolution_LP() {
  return
    GrbGetDblAttrArray_VarCon(model_fixed_, GRB_DBL_ATTR_PI, 1);
}

std::vector<double> GurobiBackend::GurobiDualSolution_QCP() {
  int nqc;
  GRBgetintattr(model_fixed_, GRB_INT_ATTR_NUMQCONSTRS, &nqc);
  return
    GrbGetDblAttrArray(model_fixed_, GRB_DBL_ATTR_QCPI, nqc);
}

double GurobiBackend::ObjectiveValue() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_OBJVAL);
}

ArrayRef<double> GurobiBackend::ObjectiveValues() const {
  int no = NumberOfObjectives();
  if(no==0)
    return std::vector<double>();
  std::vector<double> objs(no, std::numeric_limits<double>::quiet_NaN());

  if (NumberOfObjectives() == 1)
    objs[0] = GrbGetDblAttr(GRB_DBL_ATTR_OBJVAL);
  else {
    GRBenv* env = GRBgetenv(model_);
    int objnumber = GrbGetIntParam(GRB_INT_PAR_OBJNUMBER);
    for (int i = 0; i < no; i++)
    {
      GRBsetintparam(env, GRB_INT_PAR_OBJNUMBER, i);
      objs[i] = GrbGetDblAttr(GRB_DBL_ATTR_OBJNVAL);
    }
    GRBsetintparam(env, GRB_INT_PAR_OBJNUMBER, objnumber);
  }
  return objs;
}

ArrayRef<double> GurobiBackend::CurrentGrbPoolPrimalSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_XN, NumberOfVariables());
}

double GurobiBackend::CurrentGrbPoolObjectiveValue() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_POOLOBJVAL);
}


ArrayRef<int> GurobiBackend::VarStatii() {
  auto stt =
    GrbGetIntAttrArray(model_fixed_,
        GRB_INT_ATTR_VBASIS, NumberOfVariables());
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
        GRB_INT_ATTR_CBASIS, NumberOfConstraints());
  for (auto& s: stt) {
    switch (s) {
    case 0:
      s = (int)BasicStatus::bas;
      break;
    case -1:
      s = (int)BasicStatus::none;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Gurobi CBasis value: {}", s));
    }
  }
  return stt;
}

void GurobiBackend::VarStatii(ArrayRef<int> vst) {
  std::vector<int> stt(vst.data(), vst.data()+vst.size());
  for (auto& s: stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = 0;
      break;
    case BasicStatus::low:
      s = -1;
      break;
    case BasicStatus::upp:
      s = -2;
      break;
    case BasicStatus::sup:
      s = -3;
      break;
    case BasicStatus::none:
    case BasicStatus::equ:
    case BasicStatus::btw:
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
    case BasicStatus::none:
      s = -1;
      break;
    case BasicStatus::upp:
    case BasicStatus::sup:
    case BasicStatus::low:
    case BasicStatus::equ:
    case BasicStatus::btw:
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  GrbSetIntAttrArray(GRB_INT_ATTR_CBASIS, stt);
}

void GurobiBackend::InputPrimalDualStart(ArrayRef<double> x0, ArrayRef<double> pi0) {
  GrbSetDblAttrArray(GRB_DBL_ATTR_PSTART, x0);
  if (!IsQCP())
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
    GrbGetDblAttrArray(GRB_DBL_ATTR_UNBDRAY, NumberOfVariables());
}

ArrayRef<double> GurobiBackend::DRay() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_FARKASDUAL, NumberOfConstraints());
}


ArrayRef<int> GurobiBackend::VarsIIS() {
  auto iis_lb =
    GrbGetIntAttrArray(GRB_INT_ATTR_IIS_LB, NumberOfVariables());
  auto iis_ub =
    GrbGetIntAttrArray(GRB_INT_ATTR_IIS_UB, NumberOfVariables());
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

ArrayRef<int> GurobiBackend::ConsIIS() {
  // Adjust for non linear constraints, which always come
  // after the linear ones in the NL file
  int nl = GrbGetIntAttr(GRB_INT_ATTR_NUMSOS) +
    GrbGetIntAttr(GRB_INT_ATTR_NUMQCONSTRS) +
    GrbGetIntAttr(GRB_INT_ATTR_NUMGENCONSTRS);
  auto iis_con =
    GrbGetIntAttrArray(GRB_INT_ATTR_IIS_CONSTR,
      (std::size_t)NumberOfConstraints() + nl, nl);
  for (int i=iis_con.size(); i--; ) {
    iis_con[i] = int(iis_con[i] ? IISStatus::mem : IISStatus::non);
  }
  return iis_con;
}
double GurobiBackend::MIPGap() const {
  bool f;
  double g = GrbGetDblAttr(GRB_DBL_ATTR_MIPGAP, &f);
  return f ? g : Infinity();
}

double GurobiBackend::BestDualBound() const {
  bool f;
  double g = GrbGetDblAttr(GRB_DBL_ATTR_OBJBOUND, &f);
  return f ? g : -ModelSense() * Infinity();
}

double GurobiBackend::Kappa() const {
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

double GurobiBackend::NumberOfIterations() const {
  bool f;
  return GrbGetDblAttr(GRB_DBL_ATTR_ITERCOUNT, &f);
}

void GurobiBackend::ExportModel(const std::string &file) {
  GRB_CALL( GRBwrite(model_, file.c_str()) );
}


void GurobiBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptGurobi, model_);
}


///////////////////////////////////// SOLVE /////////////////////////////////////////
void GurobiBackend::SolveAndReportIntermediateResults() {
  PrepareGurobiSolve();

  GRB_CALL( GRBoptimize(model_) );

  if (need_multiple_solutions())
    ReportGurobiPool();
  if (need_fixed_MIP())
    ConsiderGurobiFixedModel();
}

void GurobiBackend::PrepareGurobiSolve() {
  if (need_multiple_solutions())
    GrbSetIntParam(GRB_INT_PAR_POOLSEARCHMODE, storedOptions_.nPoolMode_);
  if (need_ray_primal() || need_ray_dual())
    GrbSetIntParam(GRB_INT_PAR_INFUNBDINFO, 1);
  if (feasrelax_IOdata())
    DoGurobiFeasRelax();
}

void GurobiBackend::ReportGurobiPool() {
  if (!IsMIP())         // Gurobi 9.1.2 returns 1 solution for LP
    return;             // but cannot retrieve its pool attributes
  int iPoolSolution = -1;
  while (++iPoolSolution < GrbGetIntAttr(GRB_INT_ATTR_SOLCOUNT)) {
    GrbSetIntParam(GRB_INT_PAR_SOLUTIONNUMBER, iPoolSolution);
    ReportIntermediateSolution(
          CurrentGrbPoolObjectiveValue(),
          CurrentGrbPoolPrimalSolution(),
          {});
  }
}

void GurobiBackend::ConsiderGurobiFixedModel() {
  if (!IsMIP())
    return;
  if (IsQCP()) {
    int i=0;
    if (GRBgetintparam(env_, GRB_INT_PAR_QCPDUAL, &i) || i == 0)
      return;
  }
  if (GRBmodel* mdl = GRBfixedmodel(model_))
    model_fixed_ = mdl;
  else
    return;
  auto msg = DoGurobiFixedModel();
  if (!msg.empty()) {
    AddToSolverMessage( msg +
                        " failed in DoGurobiFixedModel()." );
    GRBfreemodel(model_fixed_);
    model_fixed_ = model_;
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
            f, "s" + (f == 1.)) );
  }
  return {};
}

void GurobiBackend::DoGurobiFeasRelax() {
  int reltype = feasrelax_IOdata().mode()-1,
      minrel = 0;
  if (reltype >= 3) {
    reltype -= 3;
    minrel = 1;
    feasrelax_IOdata().origObjAvailable_ = true;
  }
  GRB_CALL( GRBfeasrelax(model_, reltype, minrel,
                         (double*)data_or_null( feasrelax_IOdata().lbpen ),
                         (double*)data_or_null( feasrelax_IOdata().ubpen ),
                         (double*)data_or_null( feasrelax_IOdata().rhspen ),
                         &feasrelax_IOdata().origObjValue_) );
}


//////////////////////////////////////////////////////////////////////
////////////////////////// Solution Status ///////////////////////////
//////////////////////////////////////////////////////////////////////
std::string GurobiBackend::ConvertSolutionStatus(
    const mp::Interrupter &interrupter, int &solve_code) {
  namespace sol = mp::sol;
  int optimstatus;
  GRB_CALL( GRBgetintattr(model_, GRB_INT_ATTR_STATUS, &optimstatus) );
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter.Stop()) {
      solve_code = sol::INTERRUPTED;
      return "interrupted";
    }
    int solcount;
    GRB_CALL( GRBgetintattr(model_, GRB_INT_ATTR_SOLCOUNT, &solcount) );
    if (solcount>0) {
      solve_code = sol::UNCERTAIN;
      return "feasible solution";
    }
    solve_code = sol::FAILURE + 1;
    return "unknown solution status";
  case GRB_OPTIMAL:
    solve_code = sol::SOLVED;
    return "optimal solution";
  case GRB_INFEASIBLE:
    solve_code = sol::INFEASIBLE;
    return "infeasible problem";
  case GRB_UNBOUNDED:
    solve_code = sol::UNBOUNDED;
    return "unbounded problem";
  case GRB_INF_OR_UNBD:
    solve_code = sol::INFEASIBLE + 1;
    return "infeasible or unbounded problem";
  case GRB_NUMERIC:
    solve_code = sol::FAILURE;
    return "error";
  }
}


void GurobiBackend::ComputeIIS() {
  GRB_CALL(GRBcomputeIIS(model_));
}


//////////////////////////////////////////////////////////////////////////
/////////////////////////// Modeling interface ///////////////////////////
//////////////////////////////////////////////////////////////////////////
void GurobiBackend::InitProblemModificationPhase() {
  stats.time = steady_clock::now();
}

void GurobiBackend::AddVariable(Variable var) {
  char vtype = var::Type::CONTINUOUS==var.type() ?
        GRB_CONTINUOUS : GRB_INTEGER;
  auto lb=var.lb(), ub=var.ub();
  GRB_CALL( GRBaddvars(model_, 1, 0,
                       NULL, NULL, NULL, NULL,                  // placeholders, no matrix here
                       &lb, &ub, &vtype, NULL) );
}

void GurobiBackend::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (1>iobj) {
    GrbSetIntAttr( GRB_INT_ATTR_MODELSENSE,
                  obj::Type::MAX==lo.obj_sense() ? GRB_MAXIMIZE : GRB_MINIMIZE);
    NoteGurobiMainObjSense(lo.obj_sense());
    GrbSetDblAttrList( GRB_DBL_ATTR_OBJ, lo.vars(), lo.coefs() );
  } else {
    GRB_CALL( GRBsetobjectiven(model_, iobj, 0,           // default priority 0
                               /// Gurobi allows opposite sense by weight sign
                               lo.obj_sense()==GetGurobiMainObjSense() ? 1.0 : -1.0,
                               0.0, 0.0, nullptr,
                               0.0, lo.num_terms(),
                               (int*)lo.vars().data(), (double*)lo.coefs().data()) );
  }
}

void GurobiBackend::SetQuadraticObjective(int iobj, const QuadraticObjective &qo) {
  if (1>iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    GRB_CALL( GRBaddqpterms(model_, qt.num_terms(),
                                (int*)qt.vars1(), (int*)qt.vars2(),
                            (double*)qt.coefs()) );
  } else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void GurobiBackend::AddConstraint( const LinearConstraint& lc ) {
  GRB_CALL( GRBaddrangeconstr(model_, lc.nnz(),
                              (int*)lc.pvars(), (double*)lc.pcoefs(),
                              lc.lb(), lc.ub(), NULL) );
}

void GurobiBackend::AddConstraint( const QuadraticConstraint& qc ) {
  const auto& qt = qc.GetQPTerms();
  if (qc.lb()==qc.ub())
    GRB_CALL( GRBaddqconstr(model_, qc.nnz(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                            qt.num_terms(), (int*)qt.vars1(), (int*)qt.vars2(),
                            (double*)qt.coefs(), GRB_EQUAL, qc.lb(), NULL) );
  else {            // Let solver deal with lb>~ub etc.
    if (qc.lb()>MinusInfinity()) {
      GRB_CALL( GRBaddqconstr(model_, qc.nnz(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                              qt.num_terms(), (int*)qt.vars1(), (int*)qt.vars2(),
                              (double*)qt.coefs(), GRB_GREATER_EQUAL, qc.lb(), NULL) );
    }
    if (qc.ub()<Infinity()) {
      GRB_CALL( GRBaddqconstr(model_, qc.nnz(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                              qt.num_terms(), (int*)qt.vars1(), (int*)qt.vars2(),
                              (double*)qt.coefs(), GRB_LESS_EQUAL, qc.ub(), NULL) );
    }
  }
}

void GurobiBackend::AddConstraint(const MaximumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMax(model_, NULL,
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               MinusInfinity()) );
}

void GurobiBackend::AddConstraint(const MinimumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMin(model_, NULL,
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               Infinity()) );
}

void GurobiBackend::AddConstraint(const AbsConstraint &absc)  {
  const auto& args = absc.GetArguments();
  GRB_CALL( GRBaddgenconstrAbs(model_, NULL,
                               absc.GetResultVar(), args[0]) );
}

void GurobiBackend::AddConstraint(const ConjunctionConstraint &cc)  {
  const auto& args = cc.GetArguments();
  GRB_CALL( GRBaddgenconstrAnd(model_, NULL,
                               cc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiBackend::AddConstraint(const DisjunctionConstraint &dc)  {
  const auto& args = dc.GetArguments();
  GRB_CALL( GRBaddgenconstrOr(model_, NULL,
                               dc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiBackend::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model_, NULL,
                               ic.b_, ic.bv_, (int)ic.c_.size(),
                               ic.v_.data(), ic.c_.data(), GRB_LESS_EQUAL, ic.rhs_ ) );
}

//////////////////// General constraints /////////////////////
void GurobiBackend::AddConstraint(const SOS1Constraint &sos)  {
  int type = GRB_SOS_TYPE1;
  int beg = 0;
  GRB_CALL( GRBaddsos(model_, 1, sos.size(), &type, &beg,
              (int*)sos.get_vars().data(),
                      (double*)sos.get_weights().data()) );
}

void GurobiBackend::AddConstraint(const SOS2Constraint &sos)  {
  int type = GRB_SOS_TYPE2;
  int beg = 0;
  GRB_CALL( GRBaddsos(model_, 1, sos.size(), &type, &beg,
              (int*)sos.get_vars().data(),
                      (double*)sos.get_weights().data()) );
}

void GurobiBackend::AddConstraint(const ExpConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExp(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const ExpAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExpA(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiBackend::AddConstraint(const LogConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLog(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const LogAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLogA(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiBackend::AddConstraint(const PowConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrPow(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiBackend::AddConstraint(const SinConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrSin(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const CosConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrCos(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const TanConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrTan(model_, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const PLConstraint& plc) {
  PLPoints plp(plc.GetParameters());
  GRB_CALL( GRBaddgenconstrPWL(model_, NULL,
              plc.GetArguments()[0], plc.GetResultVar(),
              plp.x_.size(), plp.x_.data(), plp.y_.data()) );
}


///////////////////////////////////////////////////////
void GurobiBackend::FinishProblemModificationPhase() {
  // Update before adding statuses etc
  GRB_CALL( GRBupdatemodel(model_) );

  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
}


///////////////////////////////////////////////////////////////
////////////////////////// OPTIONS ////////////////////////////

// static possible values with descriptions

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

static const mp::OptionValueInfo values_iismethod[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Often faster than method 1", 0},
    { "1", "Can find a smaller IIS than method 0", 1},
    { "2", "Ignore the bound constraints.", 2},
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

static const mp::OptionValueInfo values_multiobjmethod[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Primal simplex", 0},
    { "1", "Dual simplex", 1},
    {"2", "Ignore warm-start information; use the algorithm "
        "specified by the method keyword.", 2}
};

static const mp::OptionValueInfo values_multiobjpre[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Do not use Gurobi's presolve", 0},
    { "1", "Conservative presolve", 1},
    {"2", "Aggressive presolve, which may degrade lower priority objectives.", 2}
};

static const mp::OptionValueInfo values_nodemethod[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Primal simplex", 0},
    { "1", "Dual simplex", 1},
    {"2", "Barrier.", 2}
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

static const mp::OptionValueInfo values_pool_mode[] = {
    {"0", "Just collect solutions during normal solve, and sort them best-first", 0},
    { "1", "Make some effort at finding additional solutions", 1},
    { "2", "Seek \"poollimit\" best solutions (default)."
      "'Best solutions' are defined by the poolgap(abs) parameters.", 2}
};

static const mp::OptionValueInfo values_varbranch[] = {
    {"-1", "Automatic choice (default)",-1},
    { "0", "Pseudo reduced - cost branching",0},
    { "1", "Pseudo shadow - price branching",1},
    { "2", "Maximum infeasibility branching",2},
    { "3", "Strong branching.",3}
};

static constexpr int PrmCutsMin=-1, PrmCutsMax=3;


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void GurobiBackend::InitCustomOptions() {

  set_option_header(
      fmt::format("Gurobi Optimizer Options for AMPL\n"
                  "---------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``{0}_options``. For example::\n"
      "\n"
      "  ampl: option {0}_options 'opttol=1e-6';\n",
                  GetSolverInvocationName()).c_str());


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
      "integer variables and quadratic constraints, "
      "alg:basis=0 is assumed quietly unless qcp:dual=1 "
      "is specified.");

  AddToOptionDescription("alg:sens",
                         "For problems with both integer variables and quadratic constraints, "
                         "alg:sens=0 is assumed quietly.");


  AddSolverOption("alg:method method lpmethod simplex",
    "Which algorithm to use for non-MIP problems or for the root node of MIP problems:\n"
    "\n.. value-table::\n", GRB_INT_PAR_METHOD, values_method, -1);


  AddSolverOption("alg:feasrelaxbigm feasrelaxbigm",
                  "Value of \"big-M\" sometimes used with constraints when doing "
                  "a feasibility relaxation.  Default = 1e6.",
                  GRB_DBL_PAR_FEASRELAXBIGM, 0.0, DBL_MAX);

  AddSolverOption("bar:convtol barconvtol",
    "Tolerance on the relative difference between the primal and dual objectives "
    "for stopping the barrier algorithm "
    "(default 1e-8).", GRB_DBL_PAR_BARCONVTOL, 0.0, 1.0);


  AddSolverOption("bar:corr barcorrectors",
    "Limit on the number of central corrections done in each barrier iteration"
    "(default -1 = automatic choice).", GRB_INT_PAR_BARCORRECTORS, -1, GRB_MAXINT);

  AddSolverOption("bar:homog barhomogeneous",
    "Whether to use the homogeneous barrier algorithm (e.g., when method=2 is specified):\n"
    "\n.. value-table::\n"
    "The homogeneous barrier algorithm can detect infeasibility or unboundedness directly, "
    "without crossover, but is a bit slower than the nonhomogeneous barrier algorithm.",
    GRB_INT_PAR_BARHOMOGENEOUS, values_barhomogeneous, -1);


  AddSolverOption("bar:iterlim bariterlim",
    "Limit on the number of barrier iterations (default 1000).",
    GRB_INT_PAR_BARITERLIMIT, 0, GRB_MAXINT);

  AddSolverOption("bar:order barorder", "Ordering used to reduce fill in sparse-matrix factorizations during the barrier algorithm. Possible values:\n"
    "\n.. value-table::\n", GRB_INT_PAR_BARORDER, values_barorder, -1);

  AddSolverOption("bar:qcptol barqcptol",
    "Convergence tolerance on the relative difference between primal and dual objective values for barrier algorithms when solving problems "
    "with quadratic constraints (default 1e-6).", GRB_DBL_PAR_BARQCPCONVTOL,
    0.0, 1.0);



  AddSolverOption("mip:bestbndstop bestbndstop",
    "Stop once the best bound on the objective value "
    "is at least as good as this value.",
    GRB_DBL_PAR_BESTBDSTOP, MinusInfinity(), Infinity());

  AddSolverOption("mip:bestobjstop bestobjstop",
    "Stop after a feasible solution with objective value "
    "at least as good as this value has been found.",
    GRB_DBL_PAR_BESTOBJSTOP, MinusInfinity(), Infinity());

  AddSolverOption("mip:bqpcuts bqpcuts",
    "Whether to enable Boolean Quadric Polytope cut generation:\n"
    "\n.. value-table::\n"
    "Overrides the \"cuts\" keyword.",
    GRB_INT_PAR_BQPCUTS, values_bqpcuts, -1);

  AddSolverOption("mip:branchdir branchdir",
    "Which child node to explore first when branching:\n"
    "\n.. value-table::",
    GRB_INT_PAR_BRANCHDIR, values_branchdir, 0);

  AddSolverOption("mip:cliquecuts cliquecuts",
    "Overrides \"cuts\"; choices as for \"cuts\".",
    GRB_INT_PAR_CLIQUECUTS, PrmCutsMin, PrmCutsMax);


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


  AddSolverOption("mip:focus mipfocus",
    "MIP solution strategy:\n" "\n.. value-table::\n",
    GRB_INT_PAR_MIPFOCUS, values_mipfocus, 0);

  AddSolverOption("mip:heurfrac heurfrac",
    "Fraction of time to spend in MIP heuristics (default 0.05).",
    GRB_DBL_PAR_HEURISTICS, 0.05, 1.0);

  AddSolverOption("alg:iismethod iismethod",
    "Which method to use when finding an IIS (irreducible infeasible "
    "set of constraints, including variable bounds):\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_IISMETHOD, values_iismethod, -1);

  AddStoredOption("mip:start mipstart intstart",
    "Whether to use initial guesses in problems with "
    "integer variables:\n"   "\n.. value-table::\n",
    storedOptions_.nMIPStart_, values_mipstart_);

  AddToOptionDescription("alg:start",
                         "MIP-specific options can be tuned via mip:start.");

  AddSolverOption("mip:maxmipsub maxmipsub",
    "Maximum number of nodes for RIMS heuristic to explore on MIP problems (default 500).",
      GRB_INT_PAR_SUBMIPNODES, 500, GRB_MAXINT);

  AddSolverOption("mip:gap mipgap",
    "Max relative MIP optimality gap (default 1e-4)",
    GRB_DBL_PAR_MIPGAP, 1e-4, DBL_MAX);

  AddSolverOption("mip:gapabs mipgapabs",
    "Max absolute MIP optimality gap (default 1e-10).",
    GRB_DBL_PAR_MIPGAPABS, 1e-10, DBL_MAX);

  AddSolverOption("mip:opttol opttol optimalitytolerance",
      "Dual feasibility tolerance.",
      GRB_DBL_PAR_OPTIMALITYTOL, 1e-9, 1e-2);

  AddSolverOption("mip:intfocus integralityfocus intfocus",
                  "Setting this parameter to 1 requests the solver to work "
                  "harder at finding solutions that are still (nearly) feasible "
                  "when all integer variables are rounded to exact integral "
                  "values to avoid numerical issues such as trickle flow:\n"
                  "\n.. value-table::\n",
      GRB_INT_PAR_INTEGRALITYFOCUS, values_01_noyes_0default_, 0);

  AddToOptionDescription("mip:round",
                         "For problems with numerical issues such as trickle flow, "
                         "option \"mip:intfocus\" can be more reliable.");



  AddStoredOption("mip:fixedmethod fixedmethod",
          "Value of \"method\" to use when seeking a basis for MIP problems "
          "when \"mip:basis=1\". Default: if \"method\" is 0 or 1 "
          "then \"method\" else 1.",
          storedOptions_.nFixedMethod_);

  AddSolverOption("mip:minrelnodes minrelnodes",
    "Number of nodes for the Minimum Relaxation heuristic to "
    "explore at the MIP root node when a feasible solution has not "
    "been found by any other heuristic; default -1 ==> automatic choice.",
    GRB_INT_PAR_MINRELNODES, -1, GRB_MAXINT);

  AddSolverOption("mip:nodemethod nodemethod",
    "Algorithm used to solve relaxed MIP node problems:\n"
    "\n.. value-table::\n", GRB_INT_PAR_NODEMETHOD, values_nodemethod, -1);

  AddSolverOption("mip:varbranch varbranch",
    "MIP branch variable selection strategy:\n"
    "\n.. value-table::\n", GRB_INT_PAR_VARBRANCH, values_varbranch, -1);



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
                           "keywords of the form obj_n_name, such as obj_1_method for the "
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
    GRB_INT_PAR_MULTIOBJPRE, values_multiobjmethod, -1);



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


  AddSolverOption("qp:nonconvex nonconvex",
                  "How to handle non-convex quadratic objectives and constraints:\n"
                  "\n.. value-table::\n",
    GRB_INT_PAR_NONCONVEX, values_nonconvex, -1);

  AddSolverOption("qcp:dual qcpdual",
                  "Whether to compute dual variables when the problem "
                  "has quadratic constraints (which can be expensive):\n"
                  "\n.. value-table::\n",
    GRB_INT_PAR_QCPDUAL, values_01_noyes_0default_, 0);



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
                        "to request finding several alternative solutions:\n"
                        "\n.. value-table::\n",
          storedOptions_.nPoolMode_, values_pool_mode);
  AddOptionSynonymsFront("ams_stub", "sol:stub");

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



  AddSolverOption("lim:iter iterlim iterlimit",
    "Iteration limit (default: no limit).",
    GRB_DBL_PAR_ITERATIONLIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:nodes nodelim nodelimit",
    "Maximum MIP nodes to explore (default: no limit).",
    GRB_DBL_PAR_NODELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:time timelim timelimit",
      "Limit on solve time (in seconds; default: no limit).",
      GRB_DBL_PAR_TIMELIMIT, 0.0, DBL_MAX);


  AddSolverOption("tech:outlev outlev",
      "0*/1: Whether to write gurobi log lines (chatter) to stdout and to file.",
    GRB_INT_PAR_OUTPUTFLAG, 0, 1);

  AddSolverOption("tech:logfreq logfreq outfreq",
      "Interval in seconds between log lines (default 5).",
    GRB_INT_PAR_DISPLAYINTERVAL, 1, GRB_MAXINT);

  AddSolverOption("tech:logfile logfile",
      "Log file name.",
      GRB_STR_PAR_LOGFILE);

  AddSolverOption("tech:threads threads",
      "How many threads to use when using the barrier algorithm\n"
      "or solving MIP problems; default 0 ==> automatic choice.",
      GRB_INT_PAR_THREADS, 0, GRB_MAXINT);

  AddStoredOption("tech:writeprob writeprob exportfile",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);
}

int GurobiBackend::GrbGetIntParam(const char *key) const {
  int v;
  GetSolverOption(key, v);
  return v;
}
double GurobiBackend::GrbGetDblParam(const char *key) const {
  double v;
  GetSolverOption(key, v);
  return v;
}
std::string GurobiBackend::GrbGetStrParam(const char *key) const {
  std::string v;
  GetSolverOption(key, v);
  return v;
}
void GurobiBackend::GrbSetIntParam(const char *key, int value) {
  SetSolverOption(key, value);
}
void GurobiBackend::GrbSetDblParam(const char *key, double value) {
  SetSolverOption(key, value);
}
void GurobiBackend::GrbSetStrParam(const char *key, const std::string& value) {
  SetSolverOption(key, value);
}

void GurobiBackend::GetSolverOption(const char *key, int &value) const {
  GRB_CALL( GRBgetintparam(GRBgetenv(model_), key, &value) );
}

void GurobiBackend::SetSolverOption(const char *key, int value) {
  GRB_CALL( GRBsetintparam(GRBgetenv(model_), key, value) );
}

void GurobiBackend::GetSolverOption(const char *key, double &value) const {
  GRB_CALL( GRBgetdblparam(GRBgetenv(model_), key, &value) );
}

void GurobiBackend::SetSolverOption(const char *key, double value) {
  GRB_CALL( GRBsetdblparam(GRBgetenv(model_), key, value) );
}

void GurobiBackend::GetSolverOption(const char *key, std::string &value) const {
  char buffer[GRB_MAX_STRLEN];
  GRB_CALL( GRBgetstrparam(GRBgetenv(model_), key, buffer) );
  value = buffer;
}

void GurobiBackend::SetSolverOption(const char *key, const std::string& value) {
  GRB_CALL( GRBsetstrparam(GRBgetenv(model_), key, value.c_str()) );
}


/// Shortcuts for attributes
int GurobiBackend::GrbGetIntAttr(const char* attr_id, bool *flag) const {
  int tmp=0;
  int error = GRBgetintattr(model_, attr_id, &tmp);
  if (flag)
    *flag = (0==error);
  else if (error)
    MP_RAISE( fmt::format("Failed to obtain attribute {}, error code {}",
                       attr_id, error ) );
  return tmp;
}

double GurobiBackend::GrbGetDblAttr(const char* attr_id, bool *flag) const {
  double tmp=0.0;
  int error = GRBgetdblattr(model_, attr_id, &tmp);
  if (flag)
    *flag = (0==error);
  else if (error)
    MP_RAISE( fmt::format("Failed to obtain attribute {}, error code {}",
                       attr_id, error ) );
  return tmp;
}

void GurobiBackend::GrbSetIntAttr(
    const char *attr_id, int val) {
  GRB_CALL( GRBsetintattr(model_, attr_id, val) );
}

void GurobiBackend::GrbSetDblAttr(
    const char *attr_id, double val) {
  GRB_CALL( GRBsetdblattr(model_, attr_id, val) );
}

std::vector<int> GurobiBackend::GrbGetIntAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset) const
{ return GrbGetIntAttrArray(model_, attr_id, size, offset); }

std::vector<int> GurobiBackend::GrbGetIntAttrArray(
    GRBmodel* mdl, const char* attr_id,
    std::size_t size, std::size_t offset) const {
  std::vector<int> res(size);
  int error = GRBgetintattrarray(mdl, attr_id,
    0, (int)(size-offset), res.data()+offset);
  if (error)
    res.clear();
  return res;
}

std::vector<double> GurobiBackend::GrbGetDblAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset) const
{ return GrbGetDblAttrArray(model_, attr_id, size, offset); }

std::vector<double> GurobiBackend::GrbGetDblAttrArray(
    GRBmodel* mdl, const char* attr_id,
    std::size_t size, std::size_t offset ) const {
  std::vector<double> res(size);
  int error = GRBgetdblattrarray(mdl, attr_id,
    0, (int)(size-offset), res.data()+offset);
  if (error)
    res.clear();
  return res;
}

std::vector<double> GurobiBackend::GrbGetDblAttrArray_VarCon(
    const char* attr, int varcon) const
{ return GrbGetDblAttrArray_VarCon(model_, attr, varcon); }

std::vector<double> GurobiBackend::GrbGetDblAttrArray_VarCon(
    GRBmodel* mdl, const char* attr, int varcon) const {
  return GrbGetDblAttrArray(mdl, attr,
                            varcon ? NumberOfConstraints() :
                                     NumberOfVariables());
}


void GurobiBackend::GrbSetIntAttrArray(
    const char *attr_id, ArrayRef<int> values, std::size_t start) {
  if (values)
    GRB_CALL( GRBsetintattrarray(model_, attr_id,
              (int)start, (int)values.size(), (int*)values.data()) );
}

void GurobiBackend::GrbSetDblAttrArray(
    const char *attr_id, ArrayRef<double> values, std::size_t start) {
  if (values)
    GRB_CALL( GRBsetdblattrarray(model_, attr_id,
              (int)start, (int)values.size(), (double*)values.data()) );
}

void GurobiBackend::GrbSetIntAttrList(const char *attr_id,
                                      const std::vector<int> &idx,
                                      const std::vector<int> &val) {
  assert(idx.size()==val.size());
  if (idx.size())
    GRB_CALL( GRBsetintattrlist(model_, attr_id,
                (int)idx.size(), (int*)idx.data(), (int*)val.data()) );
}

void GurobiBackend::GrbSetDblAttrList(const char *attr_id,
                                      const std::vector<int> &idx,
                                      const std::vector<double> &val) {
  assert(idx.size()==val.size());
  if (idx.size())
    GRB_CALL( GRBsetdblattrlist(model_, attr_id,
                                (int)idx.size(),
                                (int*)idx.data(), (double*)val.data()) );
}

void GurobiBackend::NoteGurobiMainObjSense(obj::Type s) { main_obj_sense_ = s; }

obj::Type GurobiBackend::GetGurobiMainObjSense() const { return main_obj_sense_; }

} // namespace mp

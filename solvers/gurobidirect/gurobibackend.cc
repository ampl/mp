#include <vector>
#include <climits>
#include <cfloat>

#include "gurobibackend.h"

#define GRB_CALL( call ) do { if (int e=call) RAISE( \
    fmt::format("  Call failed: '{}' with code {}", #call, e ) \
  ); } while (0)
#define RAISE(msg) throw std::runtime_error(msg)

namespace {

bool InterruptGurobi(void *model) {
  GRBterminate( static_cast<GRBmodel*>(model) );
  return true;
}

}  // namespace

namespace mp {

GurobiBackend::GurobiBackend() {
  OpenSolver();
}
GurobiBackend::~GurobiBackend() {
  CloseSolver();
}

const char* GurobiBackend::GetSolverInvocationName() { return "gurobidirect"; }
const char* GurobiBackend::GetBackendName() { return "GurobiBackend"; }

std::string GurobiBackend::GetSolverVersion() {
  int a,b,c;
  GRBversion(&a, &b, &c);
  return fmt::format("{}.{}.{}", a, b, c);
}

void GurobiBackend::OpenSolver() {
  
  GRB_CALL( GRBloadenv(&env, NULL) );

  /* Create an empty model */
  GRB_CALL( GRBnewmodel(env, &model, "amplgurobidirectmodel", 0, NULL, NULL, NULL, NULL, NULL) );

}

void GurobiBackend::CloseSolver() {
  /* Free model */
  GRBfreemodel(model);

  /* Free environment */
  GRBfreeenv(env);

}


bool GurobiBackend::IsMIP() const {
  return 1 == GrbGetIntAttr(GRB_INT_ATTR_IS_MIP);
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

std::vector<double> GurobiBackend::PrimalSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_X, NumberOfVariables());
}

std::vector<double> GurobiBackend::DualSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_PI, NumberOfConstraints());
}

double GurobiBackend::ObjectiveValue() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_OBJVAL);
}

std::vector<double> GurobiBackend::ObjectiveValues() const {
  int no = NumberOfObjectives();
  if(no==0)
    return std::vector<double>();
  std::vector<double> objs(no, std::numeric_limits<double>::quiet_NaN());

  if (NumberOfObjectives() == 1)
    objs[0] = GrbGetDblAttr(GRB_DBL_ATTR_OBJVAL);
  else {
    GRBenv* env = GRBgetenv(model);
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

void GurobiBackend::StartPoolSolutions() {
  assert(-2==iPoolSolution);
  iPoolSolution = -1;
}

bool GurobiBackend::SelectNextPoolSolution() {
  assert(iPoolSolution>=-1);
  if (!IsMIP())         // Gurobi 9.1.2 returns 1 solution for LP
    return false;       // but cannot retrieve its pool attributes
  ++iPoolSolution;
  if (iPoolSolution < GrbGetIntAttr(GRB_INT_ATTR_SOLCOUNT)) {
    GrbSetIntParam(GRB_INT_PAR_SOLUTIONNUMBER, iPoolSolution);
    return true;
  }
  return false;
}

void GurobiBackend::EndPoolSolutions() {
  assert(iPoolSolution>=-1);
  iPoolSolution = -2;
}

std::vector<double> GurobiBackend::CurrentPoolPrimalSolution() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_XN, NumberOfVariables());
}

double GurobiBackend::CurrentPoolObjectiveValue() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_POOLOBJVAL);
}


std::vector<int> GurobiBackend::VarStatii() {
  auto stt =
    GrbGetIntAttrArray(GRB_INT_ATTR_VBASIS, NumberOfVariables());
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
      RAISE(fmt::format("Unknown Gurobi VBasis value: {}", s));
    }
  }
  return stt;
}

std::vector<int> GurobiBackend::ConStatii() {
  auto stt =
    GrbGetIntAttrArray(GRB_INT_ATTR_CBASIS, NumberOfConstraints());
  for (auto& s: stt) {
    switch (s) {
    case 0:
      s = (int)BasicStatus::bas;
      break;
    case -1:
      s = (int)BasicStatus::none;
      break;
    default:
      RAISE(fmt::format("Unknown Gurobi CBasis value: {}", s));
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
      RAISE(fmt::format("Unknown AMPL var status value: {}", s));
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
      RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  GrbSetIntAttrArray(GRB_INT_ATTR_CBASIS, stt);
}

void GurobiBackend::InputSimplexStart(ArrayRef<double> x0, ArrayRef<double> pi0) {
  GrbSetDblAttrArray(GRB_DBL_ATTR_PSTART, x0);
  GrbSetDblAttrArray(GRB_DBL_ATTR_DSTART, pi0);
}

void GurobiBackend::InputMIPStart(ArrayRef<double> x0) {
  switch (mipstart()) {
  case 0: break;
  case 1:
    GrbSetDblAttrArray(GRB_DBL_ATTR_START, x0);
    break;
  case 3:
    GrbSetIntAttrArray(GRB_INT_ATTR_VARHINTPRI,
                         ReadSuffix(sufHintPri));
    /// fall through
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


std::vector<double> GurobiBackend::Ray() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_UNBDRAY, NumberOfVariables());
}

std::vector<double> GurobiBackend::DRay() {
  return
    GrbGetDblAttrArray(GRB_DBL_ATTR_FARKASDUAL, NumberOfConstraints());
}


std::vector<int> GurobiBackend::VarsIIS() {
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

std::vector<int> GurobiBackend::ConsIIS() {
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

double GurobiBackend::NodeCount() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_NODECOUNT);
}

double GurobiBackend::Niterations() const {
  return GrbGetDblAttr(GRB_DBL_ATTR_ITERCOUNT);
}

void GurobiBackend::ExportModel(const std::string &file) {
  GRB_CALL( GRBwrite(model, file.c_str()) );
}


void GurobiBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptGurobi, model);
}

void GurobiBackend::SolveAndReportIntermediateResults() {
  PrepareParameters();

  GRB_CALL( GRBoptimize(model) );
}

void GurobiBackend::PrepareParameters() {
  if (need_multiple_solutions())
    GrbSetIntParam(GRB_INT_PAR_POOLSEARCHMODE, storedOptions_.nPoolMode_);
  if (need_ray_primal() || need_ray_dual())
    GrbSetIntParam(GRB_INT_PAR_INFUNBDINFO, 1);
}

std::string GurobiBackend::ConvertSolutionStatus(
    const mp::Interrupter &interrupter, int &solve_code) {
  namespace sol = mp::sol;
  int optimstatus;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus) );
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter.Stop()) {
      solve_code = sol::INTERRUPTED;
      return "interrupted";
    }
    int solcount;
    GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solcount) );
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
  GRB_CALL(GRBcomputeIIS(model));
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
  GRB_CALL( GRBaddvars(model, 1, 0,
                       NULL, NULL, NULL, NULL,                  // placeholders, no matrix here
                       &lb, &ub, &vtype, NULL) );
}

void GurobiBackend::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (1>iobj) {
    GrbSetIntAttr( GRB_INT_ATTR_MODELSENSE,
                  obj::Type::MAX==lo.obj_sense() ? GRB_MAXIMIZE : GRB_MINIMIZE);
    SetMainObjSense(lo.obj_sense());
    GrbSetDblAttrList( GRB_DBL_ATTR_OBJ, lo.vars(), lo.coefs() );
  } else {
    GRB_CALL( GRBsetobjectiven(model, iobj, 0,           // default priority 0
                               /// Gurobi allows opposite sense by weight sign
                               lo.obj_sense()==GetMainObjSense() ? 1.0 : -1.0,
                               0.0, 0.0, nullptr,
                               0.0, lo.num_terms(),
                               (int*)lo.vars().data(), (double*)lo.coefs().data()) );
  }
}

void GurobiBackend::SetQuadraticObjective(int iobj, const QuadraticObjective &qo) {
  if (1>iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    GRB_CALL( GRBaddqpterms(model, qt.num_terms(),
                                (int*)qt.vars1(), (int*)qt.vars2(),
                            (double*)qt.coefs()) );
  } else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void GurobiBackend::AddConstraint( const LinearConstraint& lc ) {
  GRB_CALL( GRBaddrangeconstr(model, lc.nnz(),
                              (int*)lc.pvars(), (double*)lc.pcoefs(),
                              lc.lb(), lc.ub(), NULL) );
}

void GurobiBackend::AddConstraint( const QuadraticConstraint& qc ) {
  const auto& qt = qc.GetQPTerms();
  if (qc.lb()==qc.ub())
    GRB_CALL( GRBaddqconstr(model, qc.nnz(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                            qt.num_terms(), (int*)qt.vars1(), (int*)qt.vars2(),
                            (double*)qt.coefs(), GRB_EQUAL, qc.lb(), NULL) );
  else {            // Let solver deal with lb>~ub etc.
    if (qc.lb()>MinusInfinity()) {
      GRB_CALL( GRBaddqconstr(model, qc.nnz(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                              qt.num_terms(), (int*)qt.vars1(), (int*)qt.vars2(),
                              (double*)qt.coefs(), GRB_GREATER_EQUAL, qc.lb(), NULL) );
    }
    if (qc.ub()<Infinity()) {
      GRB_CALL( GRBaddqconstr(model, qc.nnz(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                              qt.num_terms(), (int*)qt.vars1(), (int*)qt.vars2(),
                              (double*)qt.coefs(), GRB_LESS_EQUAL, qc.ub(), NULL) );
    }
  }
}

void GurobiBackend::AddConstraint(const MaximumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMax(model, NULL,
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               MinusInfinity()) );
}

void GurobiBackend::AddConstraint(const MinimumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMin(model, NULL,
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               Infinity()) );
}

void GurobiBackend::AddConstraint(const AbsConstraint &absc)  {
  const auto& args = absc.GetArguments();
  GRB_CALL( GRBaddgenconstrAbs(model, NULL,
                               absc.GetResultVar(), args[0]) );
}

void GurobiBackend::AddConstraint(const ConjunctionConstraint &cc)  {
  const auto& args = cc.GetArguments();
  GRB_CALL( GRBaddgenconstrAnd(model, NULL,
                               cc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiBackend::AddConstraint(const DisjunctionConstraint &dc)  {
  const auto& args = dc.GetArguments();
  GRB_CALL( GRBaddgenconstrOr(model, NULL,
                               dc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiBackend::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model, NULL,
                               ic.b_, ic.bv_, (int)ic.c_.size(),
                               ic.v_.data(), ic.c_.data(), GRB_LESS_EQUAL, ic.rhs_ ) );
}

//////////////////// Nonlinear /////////////////////
void GurobiBackend::AddConstraint(const ExpConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExp(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const ExpAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExpA(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiBackend::AddConstraint(const LogConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLog(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const LogAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLogA(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiBackend::AddConstraint(const PowConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrPow(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiBackend::AddConstraint(const SinConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrSin(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const CosConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrCos(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const TanConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrTan(model, NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiBackend::AddConstraint(const PLConstraint& plc) {
  PLPoints plp(plc.GetParameters());
  GRB_CALL( GRBaddgenconstrPWL(model, NULL,
              plc.GetArguments()[0], plc.GetResultVar(),
              plp.x_.size(), plp.x_.data(), plp.y_.data()) );
}


///////////////////////////////////////////////////////
void GurobiBackend::FinishProblemModificationPhase() {
  // Update before adding statuses etc
  GRB_CALL( GRBupdatemodel(model) );

  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
}


///////////////////////////////////////////////////////////////
////////////////////////// OPTIONS ////////////////////////////

// static possible values with descriptions

const mp::OptionValueInfo values_barhomogeneous[] = {
    {"-1", "Only when solving a MIP node relaxation (default)", -1},
    { "0", "Never", 0},
    { "1", "Always", 1}
};

const mp::OptionValueInfo values_barorder[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Approximate minimum degree", 0},
    { "1", "Nested dissection", 1}
};

const mp::OptionValueInfo values_iismethod[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Often faster than method 1", 0},
    { "1", "Can find a smaller IIS than method 0", 1},
    { "2", "Ignore the bound constraints", 2},
};

const mp::OptionValueInfo values_method[] = {
    { "-1", "Automatic (default): 3 for LP, 2 for !P, 1 for MIP", -1},
    { "0", "Primal simplex", 0},
    { "1", "Dual simplex", 1},
    { "2", "Barrier", 2},
    { "3", "Nondeterministic concurrent (several solves in parallel)", 3},
    { "4", "Deterministic concurrent", 4},
    { "5", "Deterministic concurrent simplex", 5}
};

const mp::OptionValueInfo values_mipfocus[] = {
    { "0", "Balance finding good feasible solutions and "
			    "proving optimality (default)", 0},
    { "1", "Favor finding feasible solutions", 1},
    { "2", "Favor providing optimality", 2},
    { "3", "Focus on improving the best objective bound", 3},
};

const mp::OptionValueInfo values_mipstart_[4] = {
    {     "0", "No (overrides mp:warmstart)", 0 },
    {     "1", "Yes (default)", 1},
    {     "2", "No, but use the incoming primal "
          "values as hints (VARHINTVAL), ignoring the .hintpri suffix", 2},
    {     "3", "Similar to 2, but use the .hintpri suffix on "
          "variables:  larger (integer) values give greater "
          "priority to the initial value of the associated "
          "variable.", 3}
};

const mp::OptionValueInfo values_multiobjmethod[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Primal simplex", 0},
    { "1", "Dual simplex", 1},
    {"2", "Ignore warm-start information; use the algorithm "
				"specified by the method keyword", 2}
};

const mp::OptionValueInfo values_multiobjpre[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Do not use Gurobi's presolve", 0},
    { "1", "Conservative presolve", 1},
    {"2", "Aggressive presolve, which may degrade lower priority objectives", 2}
};

const mp::OptionValueInfo values_nodemethod[] = {
    {"-1", "Automatic choice (default)", -1},
    { "0", "Primal simplex", 0},
    { "1", "Dual simplex", 1},
    {"2", "Barrier", 2}
};

const mp::OptionValueInfo values_predeprow[] = {
    { "-1", "Only for continuous models (default)", -1},
    { "0", "Never", 0},
    { "1", "For all models", 1}
};

const mp::OptionValueInfo values_predual[] = {
    { "-1", "Automatic choice (default)", -1},
    { "0", "No", 0},
    { "1", "Yes", 1},
    { "2", "Form both primal and dual and use two threads to "
			     "choose heuristically between them.", 2}
};

const mp::OptionValueInfo values_pool_mode[] = {
    {"0", "Just collect solutions during normal solve, and sort them best-first", 0},
    { "1", "Make some effort at finding additional solutions", 1},
    { "2", "Seek \"poollimit\" best solutions (default)."
      "'Best solutions' are defined by the poolgap(abs) parameters.", 2}
};

const mp::OptionValueInfo values_varbranch[] = {
    {"-1", "Automatic choice (default)",-1},
    { "0", "Pseudo reduced - cost branching",0},
    { "1", "Pseudo shadow - price branching",1},
    { "2", "Maximum infeasibility branching",2},
    { "3", "Strong branching.",3}
};



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

  AddSolverOption("pre:aggfill aggfill", "Amount of fill allowed during aggregation in presolve"
    "(default -1).", GRB_INT_PAR_AGGFILL, -1, INT_MAX);

  
  AddSolverOption("pre:aggregate aggregate", "0/1*: whether to use aggregation in presolve."
    "Setting it to 0 can sometimes reduce numerical errors.", GRB_INT_PAR_AGGREGATE, 0, 1);

  AddSolverOption("pre:deprow predeprow",
    "Whether Gurobi's presolve should remove linearly dependent "
    "constraint-matrix rows:\n" "\n.. value-table::\n",
    GRB_INT_PAR_PREDEPROW, values_predeprow, -1);

  AddSolverOption("pre:dual predual",
    "Whether gurobi's presolve should form the dual of a "
    "continuous model:\n" "\n.. value-table::\n",
    GRB_INT_PAR_PREDUAL, values_predual, -1);
  
  AddSolverOption("bar:convtol barconvtol",
    "Tolerance on the relative difference between the primal and dual objectives "
    "for stopping the barrier algorithm "
    "(default 1e-8).", GRB_DBL_PAR_BARCONVTOL, 0.0, 1.0);


  AddSolverOption("bar:corr barcorrectors",
    "Limit on the number of central corrections done in each barrier iteration"
    "(default -1 = automatic choice).", GRB_INT_PAR_BARCORRECTORS, -1, INT_MAX);

  AddSolverOption("bar:homog barhomogeneous",
    "Whether to use the homogeneous barrier algorithm (e.g., when method=2 is specified):\n"
    "\n.. value-table::\n"
    "The homogeneous barrier algorithm can detect infeasibility or unboundedness directly, "
    "without crossover, but is a bit slower than the nonhomogeneous barrier algorithm.",
    GRB_INT_PAR_BARHOMOGENEOUS, values_barhomogeneous, -1);


  AddSolverOption("bar:iterlim bariterlim",
    "Limit on the number of barrier iterations (default 1000).", 
    GRB_INT_PAR_BARITERLIMIT, 0, INT_MAX);

  AddSolverOption("bar:order barorder", "Ordering used to reduce fill in sparse-matrix factorizations during the barrier algorithm. Possible values:\n"
    "\n.. value-table::\n", GRB_INT_PAR_AGGREGATE, values_barorder, -1);

  AddSolverOption("bar:qcptol barqcptol",
    "Convergence tolerance on the relative difference between primal and dual objective values for barrier algorithms when solving problems "
    "with quadratic constraints (default 1e-6).", GRB_DBL_PAR_BARQCPCONVTOL, 
    0.0, 1.0);

  AddSolverOption("log:file logfile",
      "Log file name.",
      GRB_STR_PAR_LOGFILE);




  /// Option "multiobj" is created internally if
  /// ThisBackend::IfMultipleObj() returns true.
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

  AddSolverOption("mip:focus mipfocus",
    "MIP solution strategy:\n" "\n.. value-table::\n",
    GRB_INT_PAR_MIPFOCUS, values_mipfocus, 0);

  AddSolverOption("mip:heurfrac heurfrac",
    "Fraction of time to spend in MIP heuristics (default 0.05)",
    GRB_DBL_PAR_HEURISTICS, 0.05, 1.0);

  AddSolverOption("mip:iismethod iismethod",
    "Which method to use when finding an IIS (irreducible infeasible "
    "set of constraints, including variable bounds):\n"
    "\n.. value-table::\n",
    GRB_INT_PAR_IISMETHOD, values_iismethod, -1);

  AddStoredOption("mip:start mipstart",
    "Whether to use initial guesses in problems with "
    "integer variables:\n"   "\n.. value-table::\n",
    storedOptions_.nMIPStart_, values_mipstart_);

  AddSolverOption("mip:maxmipsub maxmipsub",
    "Maximum number of nodes for RIMS heuristic to explore on MIP problems (default 500).",
      GRB_INT_PAR_SUBMIPNODES, 500, INT_MAX);

  AddSolverOption("mip:gap mipgap",
    "Max relative MIP optimality gap (default 1e-4)",
    GRB_DBL_PAR_MIPGAP, 1e-4, DBL_MAX);

  AddSolverOption("mip:gapabs mipgapabs",
    "Max absolute MIP optimality gap (default 1e-10)",
    GRB_DBL_PAR_MIPGAPABS, 1e-10, DBL_MAX);

  AddSolverOption("mip:opttol optimalitytolerance opttol",
      "Dual feasibility tolerance.",
      GRB_DBL_PAR_OPTIMALITYTOL, 1e-9, 1e-2);

  AddSolverOption("mip:minrelnodes minrelnodes",
    "Number of nodes for the Minimum Relaxation heuristic to "
		"explore at the MIP root node when a feasible solution has not "
		"been found by any other heuristic; default -1 ==> automatic choice.",
    GRB_INT_PAR_MINRELNODES, -1, INT_MAX);

  AddSolverOption("mip:nodemethod nodemethod",
    "Algorithm used to solve relaxed MIP node problems:\n"
    "\n.. value-table::\n", GRB_INT_PAR_NODEMETHOD, values_nodemethod, -1);

  AddSolverOption("mip:varbranch varbranch",
    "MIP branch variable selection strategy:\n"
    "\n.. value-table::\n", GRB_INT_PAR_VARBRANCH, values_varbranch, -1);



  AddSolverOption("log:lev outlev",
      "0*/1: Whether to write gurobi log lines (chatter) to stdout and to file.", 
    GRB_INT_PAR_OUTPUTFLAG, 0, 1);
  SetSolverOption(GRB_INT_PAR_OUTPUTFLAG, 0);

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
  AddSolverOption("sol:poollimit ams_limit poollimit",
      "Limit on the number of alternate MIP solutions written. Default: 10.",
      GRB_INT_PAR_POOLSOLUTIONS, 1, 2000000000);
  AddStoredOption("sol:poolmode ams_mode poolmode",
      "Search mode for MIP solutions when sol:stub/sol:count are specified "
                        "to request finding several alternative solutions:\n"
                        "\n.. value-table::\n",
          storedOptions_.nPoolMode_, values_pool_mode);
  AddOptionSynonymsFront("ams_stub", "sol:stub");

  /// Option "solutionstub" is created internally if
  /// ThisBackend::IfMultipleSol() returns true.
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


  AddSolverOption("gen:method method lpmethod",
    "Which algorithm to use for non-MIP problems or for the root node of MIP problems:\n"
    "\n.. value-table::\n", GRB_INT_PAR_METHOD, values_method, -1);

  AddSolverOption("gen:threads threads",
      "How many threads to use when using the barrier algorithm\n"
      "or solving MIP problems; default 0 ==> automatic choice.",
      GRB_INT_PAR_THREADS, 0, INT_MAX);


  AddSolverOption("lim:iterlimit iterlimit",
    "Iteration limit (default: no limit).",
    GRB_DBL_PAR_ITERATIONLIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:node nodelim",
    "Maximum MIP nodes to explore (default: no limit).",
    GRB_DBL_PAR_NODELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:time timelim",
      "Limit on solve time (in seconds; default: no limit).",
      GRB_DBL_PAR_TIMELIMIT, 0.0, DBL_MAX);



  AddStoredOption("gen:writeprob writeprob exportfile",
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
  GRB_CALL( GRBgetintparam(GRBgetenv(model), key, &value) );
}

void GurobiBackend::SetSolverOption(const char *key, int value) {
  GRB_CALL( GRBsetintparam(GRBgetenv(model), key, value) );
}

void GurobiBackend::GetSolverOption(const char *key, double &value) const {
  GRB_CALL( GRBgetdblparam(GRBgetenv(model), key, &value) );
}

void GurobiBackend::SetSolverOption(const char *key, double value) {
  GRB_CALL( GRBsetdblparam(GRBgetenv(model), key, value) );
}

void GurobiBackend::GetSolverOption(const char *key, std::string &value) const {
  char buffer[GRB_MAX_STRLEN];
  GRB_CALL( GRBgetstrparam(GRBgetenv(model), key, buffer) );
  value = buffer;
}

void GurobiBackend::SetSolverOption(const char *key, const std::string& value) {
  GRB_CALL( GRBsetstrparam(GRBgetenv(model), key, value.c_str()) );
}


/// Shortcuts for attributes
int GurobiBackend::GrbGetIntAttr(const char* attr_id, bool *flag) const {
  int tmp=INT_MIN;
  int error = GRBgetintattr(model, attr_id, &tmp);
  if (flag)
    *flag = (0==error);
  else if (error)
    RAISE( fmt::format("Failed to obtain attribute {}, error code {}",
                       attr_id, error ) );
  return tmp;
}

double GurobiBackend::GrbGetDblAttr(const char* attr_id, bool *flag) const {
  double tmp=0.0;
  int error = GRBgetdblattr(model, attr_id, &tmp);
  if (flag)
    *flag = (0==error);
  else if (error)
    RAISE( fmt::format("Failed to obtain attribute {}, error code {}",
                       attr_id, error ) );
  return tmp;
}

void GurobiBackend::GrbSetIntAttr(
    const char *attr_id, int val) {
  GRB_CALL( GRBsetintattr(model, attr_id, val) );
}

void GurobiBackend::GrbSetDblAttr(
    const char *attr_id, double val) {
  GRB_CALL( GRBsetdblattr(model, attr_id, val) );
}

std::vector<int> GurobiBackend::GrbGetIntAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset) const {
  std::vector<int> res(size);
  int error = GRBgetintattrarray(model, attr_id,
    0, (int)(size-offset), res.data()+offset);
  if (error)
    res.clear();
  return res;
}

std::vector<double> GurobiBackend::GrbGetDblAttrArray(const char* attr_id,
  std::size_t size, std::size_t offset ) const {
  std::vector<double> res(size);
  int error = GRBgetdblattrarray(model, attr_id,
    0, (int)(size-offset), res.data()+offset);
  if (error)
    res.clear();
  return res;
}

void GurobiBackend::GrbSetIntAttrArray(
    const char *attr_id, ArrayRef<int> values, std::size_t start) {
  if (values)
    GRB_CALL( GRBsetintattrarray(model, attr_id,
              (int)start, (int)values.size(), (int*)values.data()) );
}

void GurobiBackend::GrbSetDblAttrArray(
    const char *attr_id, ArrayRef<double> values, std::size_t start) {
  if (values)
    GRB_CALL( GRBsetdblattrarray(model, attr_id,
              (int)start, (int)values.size(), (double*)values.data()) );
}

void GurobiBackend::GrbSetIntAttrList(const char *attr_id,
                                      const std::vector<int> &idx,
                                      const std::vector<int> &val) {
  assert(idx.size()==val.size());
  if (idx.size())
    GRB_CALL( GRBsetintattrlist(model, attr_id,
                (int)idx.size(), (int*)idx.data(), (int*)val.data()) );
}

void GurobiBackend::GrbSetDblAttrList(const char *attr_id,
                                      const std::vector<int> &idx,
                                      const std::vector<double> &val) {
  assert(idx.size()==val.size());
  if (idx.size())
    GRB_CALL( GRBsetdblattrlist(model, attr_id,
                                (int)idx.size(),
                                (int*)idx.data(), (double*)val.data()) );
}

void GurobiBackend::SetMainObjSense(obj::Type s) { main_obj_sense_ = s; }

obj::Type GurobiBackend::GetMainObjSense() const { return main_obj_sense_; }

} // namespace mp

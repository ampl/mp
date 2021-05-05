#include <vector>
#include <climits>
#include <cfloat>

#include "gurobibackend.h"

#define GRB_CALL( call ) do { if (int e=call) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)


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
  return 1 == GetGrbIntAttribute(GRB_INT_ATTR_IS_MIP);
}

bool GurobiBackend::IsQCP() const {
  return 1 == GetGrbIntAttribute(GRB_INT_ATTR_IS_QCP);
}

int GurobiBackend::NumberOfConstraints() const {
  return GetGrbIntAttribute(GRB_INT_ATTR_NUMCONSTRS);
}

int GurobiBackend::NumberOfVariables() const {
  return GetGrbIntAttribute(GRB_INT_ATTR_NUMVARS);
}

int GurobiBackend::NumberOfObjectives() const {
  return GetGrbIntAttribute(GRB_INT_ATTR_NUMOBJ);
}

std::vector<double> GurobiBackend::PrimalSolution() {
  return
    GetGrbDblArrayAttribute(GRB_DBL_ATTR_X, NumberOfVariables());
}

std::vector<double> GurobiBackend::DualSolution() {
  return
    GetGrbDblArrayAttribute(GRB_DBL_ATTR_PI, NumberOfConstraints());
}

double GurobiBackend::ObjectiveValue() const {
  return GetGrbDblAttribute(GRB_DBL_ATTR_OBJVAL);
}

std::vector<int> GurobiBackend::VarStatii() {
  return
    GetGrbIntArrayAttribute(GRB_INT_ATTR_VBASIS, NumberOfVariables());
}

std::vector<int> GurobiBackend::ConStatii() {
  return
    GetGrbIntArrayAttribute(GRB_INT_ATTR_CBASIS, NumberOfConstraints());
}

std::vector<int> GurobiBackend::VarsIIS() {
  return
    GetGrbIntArrayAttribute(GRB_INT_ATTR_IIS_LB, NumberOfVariables());
}

std::vector<int> GurobiBackend::ConsIIS() {
  // Adjust for non linear constraints, which always come
  // after the linear ones in the NL file
  int nl = GetGrbIntAttribute(GRB_INT_ATTR_NUMSOS) +
    GetGrbIntAttribute(GRB_INT_ATTR_NUMQCONSTRS) +
    GetGrbIntAttribute(GRB_INT_ATTR_NUMGENCONSTRS);
  return
    GetGrbIntArrayAttribute(GRB_INT_ATTR_IIS_CONSTR,
      (std::size_t)NumberOfConstraints() + nl, nl);
}
double GurobiBackend::MIPGap() {
  return GetGrbDblAttribute(GRB_DBL_ATTR_MIPGAP);
}

double GurobiBackend::NodeCount() const {
  return GetGrbDblAttribute(GRB_DBL_ATTR_NODECOUNT);
}

double GurobiBackend::Niterations() const {
  return GetGrbDblAttribute(GRB_DBL_ATTR_ITERCOUNT);
}

void GurobiBackend::ExportModel(const std::string &file) {
  GRB_CALL( GRBwrite(model, file.c_str()) );
}


void GurobiBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptGurobi, model);
}

void GurobiBackend::SolveAndReportIntermediateResults() {
  GRB_CALL( GRBoptimize(model) );
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
void GurobiBackend::AddLinearObjective( const LinearObjective& lo ) {
  if (1>=NumberOfObjectives()) {
    GRB_CALL( GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE,
                  obj::Type::MAX==lo.get_sense() ? GRB_MAXIMIZE : GRB_MINIMIZE) );
    GRB_CALL( GRBsetdblattrlist(model, GRB_DBL_ATTR_OBJ,
                                lo.get_num_terms(),
                                (int*)lo.get_vars().data(),
                                (double*)lo.get_coefs().data()) );
  } else {
    throw std::runtime_error("Multiple objectives not supported");
//    TODO
//    GRB_CALL( GRBsetobjectiven(model, 0, 1, 0.0, 0.0, 0.0, "primary",
//                               0.0, nnz, (int*)v, (double*)c) );
  }
}
void GurobiBackend::AddQuadraticObjective(const QuadraticObjective &qo) {
  if (1>=NumberOfObjectives()) {
    AddLinearObjective(qo);                         // add the linear part
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
                               args.size(), args.data(),
                               MinusInfinity()) );
}

void GurobiBackend::AddConstraint(const MinimumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMin(model, NULL,
                               mc.GetResultVar(),
                               args.size(), args.data(),
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
                               args.size(), args.data()) );
}

void GurobiBackend::AddConstraint(const DisjunctionConstraint &dc)  {
  const auto& args = dc.GetArguments();
  GRB_CALL( GRBaddgenconstrOr(model, NULL,
                               dc.GetResultVar(),
                               args.size(), args.data()) );
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



///////////////////////////////////////////////////////
void GurobiBackend::FinishProblemModificationPhase() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
}


///////////////////////////////////////////////////////////////
////////////////////////// OPTIONS ////////////////////////////
void GurobiBackend::InitOptions() {


  set_option_header(
      fmt::format("Gurobi Optimizer Options for AMPL\n"
                  "---------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``{0}_options``. For example::\n"
      "\n"
      "  ampl: option {0}_options 'optimalitytolerance=1e-6';\n",
                  GetSolverInvocationName()).c_str());

  AddSolverOption("outlev",
      "1: output logging (console and file). "
      "Default = 0 (no logging).", GRB_INT_PAR_OUTPUTFLAG, 0, 1);
  SetSolverOption(GRB_INT_PAR_OUTPUTFLAG, 0);

  AddSolverOption("logfile",
      "Log file name.",
      GRB_STR_PAR_LOGFILE);

  AddStoredOption("exportfile",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  AddSolverOption("optimalitytolerance",
      "Dual feasibility tolerance.",
      GRB_DBL_PAR_OPTIMALITYTOL, 1e-9, 1e-2);

  AddSolverOption("threads",
      "How many threads to use when using the barrier algorithm\n"
      "or solving MIP problems; default 0 ==> automatic choice.",
      GRB_INT_PAR_THREADS, 0, INT_MAX);

  AddSolverOption("timelim",
      "limit on solve time (in seconds; default: no limit).",
      GRB_DBL_PAR_TIMELIMIT, 0.0, DBL_MAX);


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
int GurobiBackend::GetGrbIntAttribute(const char* attr_id) const {
  int tmp;
  GRB_CALL( GRBgetintattr(model, attr_id, &tmp) );
  return tmp;
}
double GurobiBackend::GetGrbDblAttribute(const char* attr_id) const {
  double tmp;
  GRB_CALL( GRBgetdblattr(model, attr_id, &tmp) );
  return tmp;
}

std::vector<int> GurobiBackend::GetGrbIntArrayAttribute(const char* attr_id,
  std::size_t size, std::size_t offset) const {
  std::vector<int> res(size);
  int error = GRBgetintattrarray(model, attr_id,
    0, size-offset, res.data()+offset);
  if (error)
    res.clear();
  return res;
}
std::vector<double> GurobiBackend::GetGrbDblArrayAttribute(const char* attr_id,
  std::size_t size, std::size_t offset ) const {
  std::vector<double> res(size);
  int error = GRBgetdblattrarray(model, attr_id,
    0, size-offset, res.data()+offset);
  if (error)
    res.clear();
  return res;
}



} // namespace mp

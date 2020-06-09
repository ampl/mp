#include <vector>

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

GurobiBackend::GurobiBackend() :
   BaseBackend("gurobidirect", 0, 0, MULTIPLE_SOL | MULTIPLE_OBJ)
   {
  InitBackend();
  InitOptions();
}

GurobiBackend::~GurobiBackend() {
  CloseBackend();
}

void GurobiBackend::InitBackend() {
  GRB_CALL( GRBemptyenv(&env) );

  GRB_CALL( GRBsetstrparam(env, "LogFile", "gurobi.log") );

  GRB_CALL( GRBstartenv(env) );

  /* Create an empty model */
  GRB_CALL( GRBnewmodel(env, &model, "amplgurobidirectmodel", 0, NULL, NULL, NULL, NULL, NULL) );

}

void GurobiBackend::CloseBackend() {
  /* Free model */
  GRBfreemodel(model);

  /* Free environment */
  GRBfreeenv(env);
}

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

void GurobiBackend::PrimalSolution(std::vector<double> &x) {
  int num_vars = NumberOfVariables();
  x.resize(num_vars);
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, num_vars, x.data());
}

void GurobiBackend::DualSolution(std::vector<double> &pi) {
  int num_cons = NumberOfConstraints();
  pi.resize(num_cons);
  GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, num_cons, pi.data());
}

double GurobiBackend::ObjectiveValue() const {
  return GetGrbDblAttribute(GRB_DBL_ATTR_OBJVAL);
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

void GurobiBackend::DoOptimize() {
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
      solve_code = 600;
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

void GurobiBackend::InitProblemModificationPhase(const Problem &p) {
  stats.time = steady_clock::now();
}

void GurobiBackend::AddVariables(int n, double *lbs, double *ubs, var::Type *types) {
  std::vector<char> vtypes(n, GRB_CONTINUOUS);
  for (int var = 0; var < n; ++var) {
    if (types[var]!=var::Type::CONTINUOUS)
      vtypes[var] = GRB_INTEGER;
  }
  GRB_CALL( GRBaddvars(model, n, 0,
                       NULL, NULL, NULL, NULL,                  // placeholders, no matrix here
                       lbs, ubs, vtypes.data(), NULL) );
}
void GurobiBackend::AddLinearObjective( obj::Type sense, int nnz,
                         const double* c, const int* v) {
  if (1>=NumberOfObjectives()) {
    GRB_CALL( GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE,
                          obj::Type::MAX==sense ? GRB_MAXIMIZE : GRB_MINIMIZE) );
    GRB_CALL( GRBsetdblattrlist(model, GRB_DBL_ATTR_OBJ, nnz, (int*)v, (double*)c) );
  } else {
//    TODO
//    GRB_CALL( GRBsetobjectiven(model, 0, 1, 0.0, 0.0, 0.0, "primary",
//                               0.0, nnz, (int*)v, (double*)c) );
  }
}
void GurobiBackend::AddLinearConstraint(int nnz, const double* c, const int* v,
                         double lb, double ub) {
//  this->Print( "  ADD LIN CONSTR:  {} <= ", lb);
//  for (int i=0; i<nnz; ++i) {
//    this->Print( "{}*[{}] ", c[i], v[i] );
//    if (i<nnz-1 && c[i+1]>=0.0)
//      this->Print( "+ " );
//  }
//  this->Print( "<= {}\n", ub );
  if (lb==ub)
    GRB_CALL( GRBaddconstr(model, nnz, (int*)v, (double*)c, GRB_EQUAL, lb, NULL) );
  else {            // Let solver deal with lb>~ub etc.
    if (lb>MinusInfinity()) {
      GRB_CALL( GRBaddconstr(model, nnz, (int*)v, (double*)c, GRB_GREATER_EQUAL, lb, NULL) );
    }
    if (ub<Infinity()) {
      GRB_CALL( GRBaddconstr(model, nnz, (int*)v, (double*)c, GRB_LESS_EQUAL, ub, NULL) );
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
  const auto& exportFile = strOptMgr_.get(EXPORT_FILE);
  if (0<exportFile.size()) {
    ExportModel(exportFile);
  }
}


///////////////////////////////////////////////////////////////
////////////////////////// OPTIONS ////////////////////////////
void GurobiBackend::InitOptions() {

  int a,b,c;
  GRBversion(&a, &b, &c);
  set_long_name(fmt::format("Gurobi {}.{}.{}", a, b, c));
  set_version(fmt::format("AMPL/Gurobi Optimizer [{}.{}.{}]", a,b,c));

  AddSuffix("priority", 0, suf::VAR);

  set_option_header(
      "Gurobi Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``gurobidirect_options``. For example::\n"
      "\n"
      "  ampl: option gurobidirect_options 'optimalitytolerance=1e-6';\n");

  AddOption("outputflag",
      "1: output logging (console and file). "
      "Default = 0 (no logging).",
      slvOptInt_, GRB_INT_PAR_OUTPUTFLAG);
  SetSolverOption(GRB_INT_PAR_OUTPUTFLAG, 0);

  AddOption("logfile",
      "Log file name.",
      slvOptString_, GRB_STR_PAR_LOGFILE);

  AddOption("exportfile",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      strOptMgr_, EXPORT_FILE);

  AddOption("optimalitytolerance",
      "Dual feasibility tolerance.",
      slvOptDouble_, GRB_DBL_PAR_OPTIMALITYTOL);

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

void GurobiBackend::SetSolverOption(const char *key, fmt::StringRef value) {
  GRB_CALL( GRBsetstrparam(GRBgetenv(model), key, value.to_string().c_str()) );
}

/////////////////////////////////////////////////////////////////////////////////////////
SolverPtr create_gurobidirect(const char *) { return SolverPtr(new GurobiBackend()); }

} // namespace mp

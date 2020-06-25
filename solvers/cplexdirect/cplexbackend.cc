#include <vector>

#include "cplexbackend.h"

#define CPLEX_CALL( call ) do { if (int e=call) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

namespace {

volatile int terminate_flag = 0;
bool InterruptCplex(void *) {
  terminate_flag = 1;
  return true;
}

}  // namespace

namespace mp {

const char* CplexBackend::GetBackendName() { return "CplexBackend"; }

CplexBackend::CplexBackend() :
   BaseBackend("cplexdirect", 0, 0, MULTIPLE_SOL | MULTIPLE_OBJ)
   {
  InitBackend();
  InitOptions();
}

CplexBackend::~CplexBackend() {
  CloseBackend();
}

void CplexBackend::InitBackend() {
  int status;
  env = CPXopenCPLEX (&status);
  if ( env == NULL ) {
     char  errmsg[CPXMESSAGEBUFSIZE];
     CPXgeterrorstring (env, status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open CPLEX environment.\n{}", errmsg ) );
  }

  CPLEX_CALL( CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON) );

  /* Create an empty model */
  lp = CPXcreateprob (env, &status, "amplcplexdirectmodel");
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create LP, error code {}.", status ) );
}

void CplexBackend::CloseBackend() {
  if ( lp != NULL ) {
     CPLEX_CALL( CPXfreeprob (env, &lp) );
  }
  /* Free up the CPLEX environment, if necessary */
  if ( env != NULL ) {
     CPLEX_CALL( CPXcloseCPLEX (&env) );
  }
}

bool CplexBackend::IsMIP() const {
  int probtype = CPXgetprobtype (env, lp);
  return
      CPXPROB_MILP == probtype ||
      CPXPROB_MIQP == probtype ||
      CPXPROB_MIQCP == probtype;
}

bool CplexBackend::IsQCP() const {
  int probtype = CPXgetprobtype (env, lp);
  return probtype >= 5;
}

int CplexBackend::NumberOfConstraints() const {
  return CPXgetnumrows (env, lp);
}

int CplexBackend::NumberOfVariables() const {
  return CPXgetnumcols (env, lp);
}

int CplexBackend::NumberOfObjectives() const {
  return CPXgetnumobjs (env, lp);
}

void CplexBackend::PrimalSolution(std::vector<double> &x) {
  int num_vars = NumberOfVariables();
  x.resize(num_vars);
  CPLEX_CALL( CPXgetx (env, lp, x.data(), 0, num_vars-1) );
}

void CplexBackend::DualSolution(std::vector<double> &pi) {
  int num_cons = NumberOfConstraints();
  pi.resize(num_cons);
  CPLEX_CALL( CPXgetpi (env, lp, pi.data(), 0, num_cons-1) );
}

double CplexBackend::ObjectiveValue() const {
  double objval;
  CPLEX_CALL( CPXgetobjval (env, lp, &objval ) );
  return objval;
}

double CplexBackend::NodeCount() const {
  return CPXgetnodecnt (env, lp);
}

double CplexBackend::Niterations() const {
  return CPXgetmipitcnt (env, lp);
}

void CplexBackend::ExportModel(const std::string &file) {
  CPLEX_CALL( CPXwriteprob (env, lp, file.c_str(), NULL) );
}


void CplexBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCplex, nullptr);
  CPLEX_CALL( CPXsetterminate (env, &terminate_flag) );
}

void CplexBackend::DoOptimize() {
  CPLEX_CALL( CPXmipopt(env, lp) );
}

std::string CplexBackend::ConvertSolutionStatus(
    const mp::Interrupter &interrupter, int &solve_code) {
  namespace sol = mp::sol;
  int optimstatus = CPXgetstat(env, lp);
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter.Stop()) {
      solve_code = 600;
      return "interrupted";
    }
    int solcount;
    solcount = CPXgetsolnpoolnumsolns (env, lp);  // Can we use it without CPXpopulate?
    if (solcount>0) {
      solve_code = sol::UNCERTAIN;
      return "feasible solution";
    }
    solve_code = sol::FAILURE + 1;
    return "unknown solution status";
  case CPX_STAT_OPTIMAL:
  case CPXMIP_OPTIMAL:
  case CPX_STAT_MULTIOBJ_OPTIMAL:
    solve_code = sol::SOLVED;
    return "optimal solution";
  case CPX_STAT_INFEASIBLE:
  case CPXMIP_INFEASIBLE:
  case CPX_STAT_MULTIOBJ_INFEASIBLE:
    solve_code = sol::INFEASIBLE;
    return "infeasible problem";
  case CPX_STAT_UNBOUNDED:
  case CPXMIP_UNBOUNDED:
  case CPX_STAT_MULTIOBJ_UNBOUNDED:
    solve_code = sol::UNBOUNDED;
    return "unbounded problem";
  case CPX_STAT_INForUNBD:
  case CPXMIP_INForUNBD:
  case CPX_STAT_MULTIOBJ_INForUNBD:
    solve_code = sol::INFEASIBLE + 1;
    return "infeasible or unbounded problem";
  }
}


void CplexBackend::InitProblemModificationPhase() {
  stats.time = steady_clock::now();
}

void CplexBackend::AddVariable(Variable var) {
  char vtype = var::Type::CONTINUOUS==var.type() ?
        CPX_CONTINUOUS : CPX_INTEGER;
  auto lb=var.lb(), ub=var.ub();
  CPLEX_CALL( CPXnewcols (env, lp, 1, NULL, &lb, &ub, &vtype, NULL) );
}

void CplexBackend::AddLinearObjective( obj::Type sense, int nnz,
                         const double* c, const int* v) {
  if (1>=NumberOfObjectives()) {
    CPLEX_CALL( CPXchgobjsen (env, lp,
                          obj::Type::MAX==sense ? CPX_MAX : CPX_MIN) );
    CPLEX_CALL( CPXchgobj (env, lp, nnz, v, c) );
  } else {
//    TODO
  }
}
void CplexBackend::AddLinearConstraint(int nnz, const double* c, const int* v,
                         double lb, double ub) {
  char sense = 'E';                     // good to initialize
  double rhs = lb;
  if (lb==ub)
    sense = 'E';
  else {            // Let solver deal with lb>~ub etc.
    if (lb>MinusInfinity()) {
      sense = 'G';
    }
    if (ub<Infinity()) {
      if ('G'==sense)
        sense = 'R';
      else {
        sense = 'L';
        rhs = ub;
      }
    }
  }
  int rmatbeg[] = { 0 };
  CPLEX_CALL( CPXaddrows (env, lp, 0, 1, nnz, &rhs,
                          &sense, rmatbeg, v, c,
                          NULL, NULL) );
  if ('R'==sense) {
    int indices = NumberOfConstraints()-1;
    double range = ub-lb;
    CPLEX_CALL( CPXchgrngval (env, lp, 1, &indices, &range) );
  }
}


void CplexBackend::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  CPLEX_CALL( CPXaddindconstr (env, lp, ic.b_, !ic.bv_, (int)ic.c_.size(),
                               ic.rhs_, 'L', ic.v_.data(), ic.c_.data(), NULL) );
}

void CplexBackend::FinishProblemModificationPhase() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
}


////////////////////////////// OPTIONS /////////////////////////////////

void CplexBackend::InitOptions() {

  int version;
  CPXversionnumber(env, &version);
  set_long_name(fmt::format("IBM ILOG CPLEX {}", version));
  set_version(fmt::format("AMPL/CPLEX Optimizer [{}]", version));

  set_option_header(
      "IBM ILOG CPLEX Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cplexdirect_options``. For example::\n"
      "\n"
      "  ampl: option cplexdirect_options 'mipgap=1e-6';\n");

  AddOption("mipdisplay",
      "0-5: output logging verbosity. "
      "Default = 0 (no logging).",
      slvOptInt_, CPXPARAM_MIP_Display);
  SetSolverOption(CPXPARAM_MIP_Display, 0);

  AddOption("exportfile",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  AddOption("mipgap",
      "Relative optimality gap |bestbound-bestinteger|/(1e-10+|bestinteger|).",
      slvOptDouble_, CPXPARAM_MIP_Tolerances_MIPGap);

}

void CplexBackend::GetSolverOption(int key, int &value) const {
  CPLEX_CALL( CPXgetintparam(env, key, &value) );
}

void CplexBackend::SetSolverOption(int key, int value) {
  CPLEX_CALL( CPXsetintparam(env, key, value) );
}

void CplexBackend::GetSolverOption(int key, double &value) const {
  CPLEX_CALL( CPXgetdblparam(env, key, &value) );
}

void CplexBackend::SetSolverOption(int key, double value) {
  CPLEX_CALL( CPXsetdblparam(env, key, value) );
}

void CplexBackend::GetSolverOption(int key, std::string &value) const {
  char buffer[CPX_STR_PARAM_MAX];
  CPLEX_CALL( CPXgetstrparam(env, key, buffer) );
  value = buffer;
}

void CplexBackend::SetSolverOption(int key, const std::string& value) {
  CPLEX_CALL( CPXsetstrparam(env, key, value.c_str()) );
}


//////////////////////////////////////////////////////////////////////////////////////
SolverPtr create_cplexdirect(const char *) { return SolverPtr(new CplexBackend()); }

} // namespace mp

#include <vector>
#include <climits>
#include <cfloat>

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

CplexBackend::CplexBackend() {
  OpenSolver();
}

CplexBackend::~CplexBackend() {
  CloseSolver();
}

const char* CplexBackend::GetBackendName()
  { return "CplexBackend"; }

std::string CplexBackend::GetSolverVersion() {
  int version;
  CPXversionnumber(env, &version);
  return fmt::format("{}", version);
}

void CplexBackend::OpenSolver() {
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

void CplexBackend::CloseSolver() {
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

int CplexBackend::NumLinCons() const {
  return CPXgetnumrows (env, lp);
}

int CplexBackend::NumVars() const {
  return CPXgetnumcols (env, lp);
}

int CplexBackend::NumObjs() const {
  return CPXgetnumobjs (env, lp);
}

Solution CplexBackend::GetSolution() {
  auto mv = GetPresolver().PostsolveSolution(
        { PrimalSolution(), DualSolution() } );
  return { mv.GetVarValues()(), mv.GetConValues()(),
    GetObjectiveValues() };   // TODO postsolve obj values
}

ArrayRef<double> CplexBackend::PrimalSolution() {
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  int error = CPXgetx (env, lp, x.data(), 0, num_vars-1);
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
  int error = CPXgetpi (env, lp, pi.data(), 0, num_cons-1);
  if (error)
    pi.clear();
  return pi;
}

double CplexBackend::ObjectiveValue() const {
  double objval = -1e308;
  CPXgetobjval (env, lp, &objval );   // failsafe
  return objval;
}

double CplexBackend::NodeCount() const {
  return CPXgetnodecnt (env, lp);
}

double CplexBackend::SimplexIterations() const {
  return std::max(
        CPXgetmipitcnt (env, lp), CPXgetitcnt (env, lp));
}

int CplexBackend::BarrierIterations() const {
  return CPXgetbaritcnt (env, lp);
}

void CplexBackend::ExportModel(const std::string &file) {
  CPLEX_CALL( CPXwriteprob (env, lp, file.c_str(), NULL) );
}


void CplexBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCplex, nullptr);
  CPLEX_CALL( CPXsetterminate (env, &terminate_flag) );
}

void CplexBackend::SolveAndReportIntermediateResults() {
  CPLEX_CALL( CPXmipopt(env, lp) );

  WindupCPLEXSolve();
}

void CplexBackend::WindupCPLEXSolve() {
  SetStatus( ConvertCPLEXStatus() );
  AddCPLEXMessages();
}

void CplexBackend::AddCPLEXMessages() {
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

std::pair<int, std::string> CplexBackend::ConvertCPLEXStatus() {
  namespace sol = mp::sol;
  int optimstatus = CPXgetstat(env, lp);
  switch (optimstatus) {
  default:
    // Fall through.
    if (interrupter()->Stop()) {
      return { sol::INTERRUPTED, "interrupted" };
    }
    int solcount;
    solcount = CPXgetsolnpoolnumsolns (env, lp);  // Can we use it without CPXpopulate?
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


void CplexBackend::InitProblemModificationPhase() { }

void CplexBackend::AddVariables(const VarArrayDef& v) {
  std::vector<char> vtypes(v.size());
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS==v.ptype()[i] ?
          CPX_CONTINUOUS : CPX_INTEGER;
  CPLEX_CALL( CPXnewcols (env, lp, (int)v.size(), NULL,
                          v.plb(), v.pub(), vtypes.data(), NULL) );
}

void CplexBackend::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (1>iobj) {
    CPLEX_CALL( CPXchgobjsen (env, lp,
                    obj::Type::MAX==lo.obj_sense() ? CPX_MAX : CPX_MIN) );
    CPLEX_CALL( CPXchgobj (env, lp, lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) );
  } else {
//    TODO
  }
}
void CplexBackend::AddConstraint(const RangeLinCon& lc) {
  char sense = 'E';                     // good to initialize things
  double rhs = lc.lb();
  if (lc.lb()==lc.ub())
    sense = 'E';
  else {                                // Let solver deal with lb>~ub etc.
    if (lc.lb()>MinusInfinity()) {
      sense = 'G';
    }
    if (lc.ub()<Infinity()) {
      if ('G'==sense)
        sense = 'R';
      else {
        sense = 'L';
        rhs = lc.ub();
      }
    }
  }
  int rmatbeg[] = { 0 };
  CPLEX_CALL( CPXaddrows (env, lp, 0, 1, lc.size(), &rhs,
                          &sense, rmatbeg, lc.pvars(), lc.pcoefs(),
                          NULL, NULL) );
  if ('R'==sense) {
    int indices = NumLinCons()-1;
    double range = lc.ub()-lc.lb();
    CPLEX_CALL( CPXchgrngval (env, lp, 1, &indices, &range) );
  }
}


void CplexBackend::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  CPLEX_CALL( CPXaddindconstr (env, lp,
                               ic.get_binary_var(), !ic.get_binary_value(),
                               (int)ic.get_constraint().size(),
                               ic.get_constraint().rhs(), 'L',
                               ic.get_constraint().pvars(),
                               ic.get_constraint().pcoefs(), NULL) );
}
void CplexBackend::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  CPLEX_CALL( CPXaddindconstr (env, lp,
                               ic.get_binary_var(), !ic.get_binary_value(),
                               (int)ic.get_constraint().size(),
                               ic.get_constraint().rhs(), 'E',
                               ic.get_constraint().pvars(),
                               ic.get_constraint().pcoefs(), NULL) );
}

void CplexBackend::FinishProblemModificationPhase() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
}


////////////////////////////// OPTIONS /////////////////////////////////

void CplexBackend::InitCustomOptions() {

  set_option_header(
      "IBM ILOG CPLEX Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cplexdirect_options``. For example::\n"
      "\n"
      "  ampl: option cplexdirect_options 'mipgap=1e-6';\n");

  AddSolverOption("tech:outlev",
      "0-5: output logging verbosity. "
      "Default = 0 (no logging).",
      CPXPARAM_MIP_Display, 0, 5);
  SetSolverOption(CPXPARAM_MIP_Display, 0);

  AddStoredOption("tech:exportfile writeprob",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

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


} // namespace mp

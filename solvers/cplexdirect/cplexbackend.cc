#include "cplexbackend.h"

#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>

#include <iostream>

using std::strcmp;
using std::vector;

namespace {


const mp::OptionValueInfo OPTIMIZERS[] = {
  {
    "auto",
    "CP Optimizer if the problem has nonlinear objective/constraints "
    "or logical constraints, CPLEX otherwise", 0
  },
  {
    "cp",
    "CP Optimizer", 0
  },
  {
    "cplex",
    "CPLEX Optimizer", 0
  }
};


mp::OptionError GetOptionValueError(
    const mp::SolverOption &opt, fmt::StringRef message) {
  throw mp::OptionError(fmt::format(
      "Can't get value of option {}: {}", opt.name(), message));
}


bool HasNonlinearObj(const mp::Problem &p) {
  if (p.num_objs() == 0)
    return false;
  mp::NumericExpr expr = p.obj(0).nonlinear_expr();
  return expr && !mp::Cast<mp::NumericConstant>(expr);
}

std::string ConvertSolutionStatus( CPXENVptr env,
    CPXLPptr lp, const mp::Interrupter &interrupter, int &solve_code) {
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

volatile int terminate_flag = 0;

bool InterruptCplex(void *) {
  terminate_flag = 1;
  return true;
}

}  // namespace

namespace mp {

CplexBackend::CplexBackend() :
   BaseSolverImpl("cplexdirect", 0, 0, MULTIPLE_SOL | MULTIPLE_OBJ)
   {
  InitBackend();

  options_[DEBUGEXPR] = false;
  options_[USENUMBEROF] = true;
  options_[SOLUTION_LIMIT] = -1;

  int version;
  CPXversionnumber(env, &version);
  set_long_name(fmt::format("IBM ILOG CPLEX {}", version));
  set_version(fmt::format("AMPL/CPLEX Optimizer [{}]", version));

  AddSuffix("priority", 0, suf::VAR);

  set_option_header(
      "IBM ILOG CPLEX Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``ilogcp_options``. For example::\n"
      "\n"
      "  ampl: option ilogcp_options 'optimalitytolerance=1e-6 "
      "searchtype=restart';\n");

  AddStrOption("optimizer",
      "Specifies which optimizer to use. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``auto``.",
      &CplexBackend::GetOptimizer, &CplexBackend::SetOptimizer, OPTIMIZERS);

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

std::string CplexBackend::GetOptimizer(const SolverOption &) const {
  switch (optimizer_) {
  default:
    assert(false);
    // Fall through.
  case AUTO:  return "auto";
  case CP:    return "cp";
  case CPLEX: return "cplex";
  }
}

void CplexBackend::SetOptimizer(const SolverOption &opt, fmt::StringRef value) {
  if (value == "auto")
    optimizer_ = AUTO;
  else if (value == "cp")
    optimizer_ = CP;
  else if (value == "cplex")
    optimizer_ = CPLEX;
  else
    throw InvalidOptionValue(opt, value);
}

void CplexBackend::SetBoolOption(
    const SolverOption &opt, int value, Option id) {
  if (value != 0 && value != 1)
    throw InvalidOptionValue(opt, value);
  options_[id] = value;
}

void CplexBackend::DoSetIntOption(
    const SolverOption &opt, int value, Option id) {
  if (value < 0)
    throw InvalidOptionValue(opt, value);
  options_[id] = value;
}


void CplexBackend::SolveWithCplex(
    Problem &p,
    Stats &stats, SolutionHandler &sh) {
  interrupter()->SetHandler(InterruptCplex, nullptr);
  CPLEX_CALL( CPXsetterminate (env, &terminate_flag) );

  std::cout << "Exporting model" << std::endl;
  ExportModel("model_cplex.lp");

  stats.setup_time = GetTimeAndReset(stats.time);
  CPLEX_CALL( CPXmipopt(env, lp) );
  stats.solution_time = GetTimeAndReset(stats.time);

  // Convert solution status.
  int solve_code = 0;
  std::string status =
      ConvertSolutionStatus(env, lp, *interrupter(), solve_code);

  fmt::MemoryWriter writer;
  writer.write("{}: {}\n", long_name(), status);
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  vector<double> solution, dual_solution;
  if (solve_code < sol::INFEASIBLE) {
    PrimalSolution(solution);

    if (IsMIP()) {
      writer << NodeCount() << " nodes, ";
    } else {                                    // Also for QCP
      DualSolution(dual_solution);
    }
    writer << Niterations() << " iterations";

    if (NumberOfObjectives() > 0) {
      writer.write(", objective {}", FormatObjValue(ObjectiveValue()));
    }
  }
  sh.HandleSolution(solve_code, writer.c_str(),
      solution.empty() ? 0 : solution.data(),
      dual_solution.empty() ? 0 : dual_solution.data(), obj_value);
}

void CplexBackend::Solve(Problem &p, SolutionHandler &sh) {
  Resolve(p, sh);
}

void CplexBackend::InitProblemModificationPhase(const Problem &p) {
  stats.time = steady_clock::now();

  optimizer = optimizer_;
  if (optimizer == AUTO) {
    if (p.num_logical_cons() != 0 || p.has_nonlinear_cons() ||
        HasNonlinearObj(p)) {
      optimizer = CP;
    } else {
      optimizer = CPLEX;
    }
  }

}

void CplexBackend::AddVariables(int n, double *lbs, double *ubs, var::Type *types) {
  std::vector<char> vtypes(n, CPX_CONTINUOUS);
  for (int var = 0; var < n; ++var) {
    if (types[var]!=var::Type::CONTINUOUS)
      vtypes[var] = CPX_INTEGER;
  }
  CPLEX_CALL( CPXnewcols (env, lp, n, NULL, lbs, ubs, vtypes.data(), NULL) );
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
}

void CplexBackend::Convert(Problem &p) {
  InitProblemModificationPhase(p);
}

void CplexBackend::Resolve(Problem &p, SolutionHandler &sh) {

  SolveWithCplex(p, stats, sh);
  double output_time = GetTimeAndReset(stats.time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          stats.setup_time, stats.solution_time, output_time);
  }
}


SolverPtr create_ilogcp(const char *) { return SolverPtr(new CplexBackend()); }
}

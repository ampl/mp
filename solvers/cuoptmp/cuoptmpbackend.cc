#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cuoptmpbackend.h"

extern "C" {
  #include "cuoptmp-ampls-c-api.h"    // cuoptmp AMPLS C API
}
#include "mp/ampls-cpp-api.h"

#include "cuopt/version_config.hpp"

namespace {


bool InterruptCuoptMp(void* prob) {
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCuoptmpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CuoptmpBackend()};
}


namespace mp {

/// Create Cuoptmp Model Manager
/// @param gc: the Cuoptmp common handle
/// @param e: environment
/// @param pre: presolver to be names[Solver::returned]= "";
/// need it to convert solution data
/// @return CuoptmptModelMgr
std::unique_ptr<BasicModelManager>
CreateCuoptmpModelMgr(CuoptmpCommon&, Env&, pre::BasicValuePresolver*&);


CuoptmpBackend::CuoptmpBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateCuoptmpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

CuoptmpBackend::~CuoptmpBackend() {
  CloseSolver();
}

void CuoptmpBackend::OpenSolver() {
  lp_ = new ProblemData();
  lp_->num_constraints = 0;
  lp_->num_variables = 0;
  lp_->constraint_matrix_row_offsets.clear();
  lp_->constraint_matrix_coefficients.clear();
  lp_->constraint_matrix_column_indices.clear();
  lp_->constraint_sense.clear();
  lp_->rhs.clear();
  lp_->lower_bounds.clear();
  lp_->upper_bounds.clear();
  lp_->variable_types.clear();
  lp_->objective_coefficients.clear();
  lp_->objective_sense = CUOPT_MINIMIZE;
  lp_->objective_offset = 0.0;
  int status = 0;
  status = cuOptCreateSolverSettings(&lp_->settings);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Failed to create solver settings, error code {}.", status));
  }
}

void CuoptmpBackend::CloseSolver() {
  /* Cleanup: close problem and environment */
  if ( lp() != NULL ) {
    cuOptDestroyProblem(&lp_->problem);
    cuOptDestroySolverSettings(&lp_->settings);
    cuOptDestroySolution(&lp_->solution);
    delete lp_;
  }
}

const char* CuoptmpBackend::GetBackendName()
  { return "CuoptmpBackend"; }

std::string CuoptmpBackend::GetSolverVersion() {
  int32_t major = CUOPT_VERSION_MAJOR, minor = CUOPT_VERSION_MINOR, patch = CUOPT_VERSION_PATCH;
  return fmt::format("{}.{}.{}", major, minor, patch);
}


bool CuoptmpBackend::IsMIP() const {
  ProblemData* problem_data = lp();
  cuopt_int_t is_mip;
  cuopt_int_t status = cuOptIsMIP(problem_data->problem, &is_mip);
  if (status != CUOPT_SUCCESS) {
    fmt::print("Error checking if problem is MIP\n");
  }
  return is_mip == 1;
 }



ArrayRef<double> CuoptmpBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error = 1;
  std::vector<double> x(num_vars);
  // We should always return a solution when available,
  // even if infeasible/suboptimal etc.
  // User decides on it using solve_result.
  ProblemData* problem_data = lp();
  cuopt_int_t status = cuOptGetPrimalSolution(problem_data->solution, x.data());
  if (status != CUOPT_SUCCESS) {
    x.clear();
  }

  for (int j = 0; j < num_vars; j++) {
    if (problem_data->variable_types[j] == CUOPT_INTEGER) {
      x[j] = std::round(x[j]);
    }
  }

  return x;
}

pre::ValueMapDbl CuoptmpBackend::DualSolution() {
  return { {
    { CG_Linear, DualSolution_LP() },
  } };
}

ArrayRef<double> CuoptmpBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  ProblemData* problem_data = lp();
  cuopt_int_t status = cuOptGetDualSolution(problem_data->solution, pi.data());
  if (status != CUOPT_SUCCESS) {
    pi.clear();
  }
  return pi;
}

double CuoptmpBackend::ObjectiveValue() const {
  ProblemData* problem_data = lp();
  cuopt_float_t objective_value;
  cuopt_int_t status = cuOptGetObjectiveValue(problem_data->solution, &objective_value);
  if (status != CUOPT_SUCCESS) {
    fmt::print("Error getting objective value\n");
  }
  return objective_value;
}

double CuoptmpBackend::NodeCount() const {
  return 0;
}

double CuoptmpBackend::SimplexIterations() const {
  return 0;
//  return getIntAttr(CUOPTLP_INTATTR_SIMPLEXITER);
}

int CuoptmpBackend::BarrierIterations() const {
  return 0;
}


void CuoptmpBackend::SetInterrupter(mp::Interrupter *inter) {
 // inter->SetHandler(InterruptCuoptmp, lp());
  // TODO Check interrupter
  //CUOPTLP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CuoptmpBackend::Solve() {
  ProblemData* problem_data = lp();

  cuopt_int_t is_mip;
  cuopt_int_t status = cuOptIsMIP(problem_data->problem, &is_mip);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Failed to check if problem is MIP, error code {}.", status));
  }
  status = cuOptSolve(problem_data->problem, problem_data->settings, &problem_data->solution);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Failed to solve problem, error code {}.", status));
  }


  WindupCUOPTLPSolve();
}

void CuoptmpBackend::WindupCUOPTLPSolve() {
}

void CuoptmpBackend::ReportResults() {
  ReportCUOPTLPResults();
  BaseBackend::ReportResults();
}

void CuoptmpBackend::ReportCUOPTLPResults() {
  SetStatus( GetSolveResult() );
  AddCUOPTLPMessages();
}

void CuoptmpBackend::AddCUOPTLPMessages() {
  if(auto si = SimplexIterations())
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", si));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> CuoptmpBackend::GetSolveResult() {
  namespace sol = mp::sol;
  /*
   * TODO.
   * Follow guidelines from mp::sol::Status.
     * Keep new result codes added
     * in AddOptions() via AddSolveResults().
     */
  ProblemData* problem_data = lp();
  cuopt_int_t termination_status;
  cuopt_int_t status = cuOptGetTerminationStatus(problem_data->solution, &termination_status);
  if (status != CUOPT_SUCCESS) {
    fmt::print("Error getting termination status\n");
  }
  switch (termination_status) {
    case CUOPT_TERIMINATION_STATUS_OPTIMAL:
      return { sol::SOLVED, "optimal" };
    case CUOPT_TERIMINATION_STATUS_INFEASIBLE:
      return { sol::INFEASIBLE_NO_IIS, "infeasible" };
    case CUOPT_TERIMINATION_STATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded" };
    case CUOPT_TERIMINATION_STATUS_ITERATION_LIMIT:
      return { sol::LIMIT_NO_FEAS_ITER, "iteration limit" }; // TODO: check if solution is feasible
    case CUOPT_TERIMINATION_STATUS_TIME_LIMIT:
      return { sol::LIMIT_NO_FEAS_TIME, "time limit" };  // TODO: check if solution is feasible
    case CUOPT_TERIMINATION_STATUS_NUMERICAL_ERROR:
      return { sol::NUMERIC, "numerical error" };
    case CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE:
      return { sol::LIMIT_FEAS, "primal feasible" };
    case CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND:
      return { sol::LIMIT_FEAS_TIME, "feasible found" };
    case CUOPT_TERIMINATION_STATUS_CONCURRENT_LIMIT:
      return { sol::LIMIT_FEAS_WORK, "concurrent limit" };
    case CUOPT_TERIMINATION_STATUS_NO_TERMINATION:
      return { sol::UNKNOWN, "not solved" };
    default:
      return { sol::UNKNOWN, "unknown termination status" };
  }
}


void CuoptmpBackend::FinishOptionParsing() {}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo verbosity_values_[] = {
  { "0", "Only statistics", 0},
  { "1", "All info", 1}
};

static const mp::OptionValueInfo values_method[] = {
    { "0", "Concurrent (default)", 0},
    { "1", "Pdlp", 1},
    { "2", "Dual simplex", 2}
};

static const mp::OptionValueInfo values_pdlp_solver_mode[] = {
    { "0", "Stable1", 0},
    { "1", "Stable2 (default)", 1},
    { "2", "Methodical1", 2},
    { "3", "Fast1", 3}
};

static const mp::OptionValueInfo values_bool[] = {
    { "true", "True", 0},
    { "false", "False", 1}
};

void CuoptmpBackend::InitCustomOptions() {

  set_option_header(
      "CUOPTLP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cuoptmp_options``. For example::\n"
      "\n"
      "  ampl: option cuoptmp_options 'mipgap=1e-6';\n");

  AddSolverOption("lim:timelimit timelimit",
    "Time limit in seconds after which the solver will stop and return the current solution",
    CUOPT_TIME_LIMIT, 0.0, DBL_MAX);

  AddSolverOption("timelim timelim",
    "Time limit in seconds after which the solver will stop and return the current solution",
    CUOPT_TIME_LIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:ncputhreads ncputhreads",
    "Number of CPU threads used in the LP and MIP solvers",
    CUOPT_NUM_CPU_THREADS, -1, INT_MAX);

  AddSolverOption("threads threads",
    "Number of CPU threads used in the LP and MIP solvers",
    CUOPT_NUM_CPU_THREADS, -1, INT_MAX);

  // LP Solver Options

  AddSolverOption("alg:method method",
    "Designate the method to solve the LP problem:\n"
    "\n.. value-table::\n",
    CUOPT_METHOD,
    values_method, CUOPT_METHOD_CONCURRENT);

  AddSolverOption("alg:solver_mode solver_mode",
    "Designate the solver mode used by PDLP to solve the problem:\n"
    "\n.. value-table::\n",
    CUOPT_PDLP_SOLVER_MODE,
    values_pdlp_solver_mode, CUOPT_PDLP_SOLVER_MODE_STABLE2);

  AddSolverOption("lim:iterations iterations",
    "Iteration limit after which the solver will stop and return the current solution",
    CUOPT_ITERATION_LIMIT, 0, INT_MAX);

  AddSolverOption("lp:infeasdetect infeasdetect",
    "Detect infeasibility PDLP",
    CUOPT_INFEASIBILITY_DETECTION,
    values_bool, "false");

  AddSolverOption("lp:strictinfeas strictinfeas",
    "Stop if current or the average solution is detected as infeasible",
    CUOPT_STRICT_INFEASIBILITY,
    values_bool, "false");

  AddSolverOption("lp:crossover crossover",
    "Crossover to a basic solution after a optimal solution is found",
    CUOPT_CROSSOVER,
    values_bool, "false");

  AddSolverOption("lp:savebestprimal savebestprimal",
    "Save the best primal solution so far",
    CUOPT_SAVE_BEST_PRIMAL_SO_FAR,
    values_bool, "false");

  AddSolverOption("lp:firstprimalfeas firstprimalfeas",
    "Stop when the first primal feasible solution is found",
    CUOPT_FIRST_PRIMAL_FEASIBLE ,
    values_bool, "false");

  AddSolverOption("lp:perconsres perconsres",
    "Compute the primal & dual residual per constraint instead of globally",
    CUOPT_PER_CONSTRAINT_RESIDUAL,
    values_bool, "false");

  AddSolverOption("lp:absprimaltol absprimaltol",
    "Absolute primal tolerance used in PDLP's primal feasibility check",
    CUOPT_ABSOLUTE_PRIMAL_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("lp:absdualtol absdualtol",
    "Absolute dual tolerance used in PDLP's dual feasibility check",
    CUOPT_ABSOLUTE_DUAL_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("lp:absgaptol absgaptol",
    "Absolute gap tolerance used in PDLP's duality gap check",
    CUOPT_ABSOLUTE_GAP_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("lp:relprimaltol relprimaltol",
    "Relative primal tolerance used in PDLP's primal feasibility check",
    CUOPT_RELATIVE_PRIMAL_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("lp:reldualtol reldualtol",
    "Relative dual tolerance used in PDLP's dual feasibility check",
    CUOPT_RELATIVE_DUAL_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("lp:relgaptol relgaptol",
    "Relative gap tolerance used in PDLP's duality gap check",
    CUOPT_RELATIVE_GAP_TOLERANCE, 0.0, 1e-1);

  // MIP Solver Options

  AddSolverOption("mip:abstol abstol",
    "Absolute tolerance used in mip",
    CUOPT_MIP_ABSOLUTE_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("mip:reltol reltol",
    "Relative tolerance used in mip",
    CUOPT_MIP_RELATIVE_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("mip:inttol inttol",
    "Integrality tolerance used in mip",
    CUOPT_MIP_INTEGRALITY_TOLERANCE, 0.0, 1e-1);

  AddSolverOption("mip:absgap absgap",
    "Absolute tolerance used to terminate the MIP solve",
    CUOPT_MIP_ABSOLUTE_GAP, 0.0, 1e-1);

  AddSolverOption("lp:relgap relgap",
    "Relative tolerance used to terminate the MIP solve",
    CUOPT_MIP_RELATIVE_GAP, 0.0, 1e-1);

  AddSolverOption("mip:hueristicsonly hueristicsonly",
    "Run only the GPU heuristics",
    CUOPT_MIP_HEURISTICS_ONLY,
    values_bool, "false");

  AddSolverOption("mip:scale scale",
    "Apply Scaling to MIP problems",
    CUOPT_MIP_SCALING,
    values_bool, "true");

  // Logging Options

  AddSolverOption("tech:consolelog consolelog",
    "Log information to the console during a solve",
    CUOPT_LOG_TO_CONSOLE,
    values_bool, "true");

  AddSolverOption("tech:outlev outlev",
    "Whether to log information to the console",
    CUOPT_LOG_TO_CONSOLE,
    values_bool, "true");
}


double CuoptmpBackend::MIPGap() {
  ProblemData* problem_data = lp();
  if (!IsMIP()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  cuopt_float_t mip_gap;
  cuopt_int_t status = cuOptGetMIPGap(problem_data->solution, &mip_gap);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Error getting MIP gap, error code {}.", status));
  }
  return mip_gap;
}

double CuoptmpBackend::BestDualBound() {
  ProblemData* problem_data = lp();
  if (!IsMIP()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  cuopt_float_t best_dual_bound;
  cuopt_int_t status = cuOptGetSolutionBound(problem_data->solution, &best_dual_bound);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Error getting best dual bound, error code {}.", status));
  }
  return best_dual_bound;
}

double CuoptmpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


} // namespace mp


// AMPLs
void* AMPLSOpenCuoptmp(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::CuoptmpBackend()},
        cb);
}

void AMPLSCloseCuoptmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetCuoptmpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CuoptmpBackend*>(AMPLSGetBackend(slv))->lp();
}

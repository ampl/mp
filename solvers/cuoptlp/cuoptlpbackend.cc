#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cuoptlpbackend.h"

extern "C" {
  #include "cuoptlp-ampls-c-api.h"    // Cuoptlp AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptCuoptlp(void* prob) {
  //return CUOPTLP_Interrupt((cuoptlp_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCuoptlpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CuoptlpBackend()};
}


namespace mp {

/// Create Cuoptlp Model Manager
/// @param gc: the Cuoptlp common handle
/// @param e: environment
/// @param pre: presolver to be names[Solver::returned]= "";
/// need it to convert solution data
/// @return CuoptlpModelMgr
std::unique_ptr<BasicModelManager>
CreateCuoptlpModelMgr(CuoptlpCommon&, Env&, pre::BasicValuePresolver*&);


CuoptlpBackend::CuoptlpBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateCuoptlpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

CuoptlpBackend::~CuoptlpBackend() {
  CloseSolver();
}

void CuoptlpBackend::OpenSolver() {
  fmt::print("opening solver\n");

  lp_ = new ProblemData();
  printf("lp_ %p this %p\n", lp_, this);

  int status = 0;
  status = cuOptCreateSolverSettings(&lp_->settings);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Failed to create solver settings, error code {}.", status));
  }

  // TODO Typically this function creates an instance of the solver environment
  // and an empty model
  void* env_p;
  // Typically try the registered function first;
  // if not available call the solver's API function directly
  /*
  const auto& create_fn = GetCallbacks().cb_initsolver_;
  if (create_fn)
    set_env((GRBenv*)create_fn());
  else
    status = createEnv(&env_p);
    */
  // set_env(env_p);

  /* Todo catch errors
  if ( env() == NULL ) {
    // char  errmsg[CPXMESSAGEBUFSIZE];
    // CPXgeterrorstring (env(), status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open CUOPTLP environment.\n{}", status) );
  }
  */


  /* TODO Create problem instance
  cuoptlp_prob* prob;
  status = CUOPTLP_CreateProb(env_p, &prob);
 */
  //Solver::SolverModel* prob = Solver::CreateSolverModel();
  //set_lp(prob); // Assign it
  //if (status)
  //  throw std::runtime_error( fmt::format(
  //        "Failed to create problem, error code {}.", status ) );
  /* TODO Typically check call */
  /// Turn off verbosity by default
  // CUOPTLP_CCALL(CUOPTLP_SetIntParam(prob, "Logging", 0));

}

void CuoptlpBackend::CloseSolver() {
  fmt::print("closing solver\n");
  /* Cleanup: close problem and environment */
  if ( lp() != NULL ) {
    cuOptDestroyProblem(&lp_->problem);
    cuOptDestroySolverSettings(&lp_->settings);
    cuOptDestroySolution(&lp_->solution);
    delete lp_;
  }
}

const char* CuoptlpBackend::GetBackendName()
  { return "CuoptlpBackend"; }

std::string CuoptlpBackend::GetSolverVersion() {
  // TODO Return version from solver API
  int32_t major, minor, patch;
  return fmt::format("{}.{}.{}", major, minor, patch);
}


bool CuoptlpBackend::IsMIP() const {
  ProblemData* problem_data = lp();
  cuopt_int_t is_mip;
  cuopt_int_t status = cuOptIsMIP(problem_data->problem, &is_mip);
  if (status != CUOPT_SUCCESS) {
    fmt::print("Error checking if problem is MIP\n");
  }
  return is_mip == 1;
 }



ArrayRef<double> CuoptlpBackend::PrimalSolution() {
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
  return x;
}

pre::ValueMapDbl CuoptlpBackend::DualSolution() {
  return { {
    { CG_Linear, DualSolution_LP() },
    { CG_Quadratic, DualSolution_QP() } } };
}

ArrayRef<double> CuoptlpBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  ProblemData* problem_data = lp();
  cuopt_int_t status = cuOptGetDualSolution(problem_data->solution, pi.data());
  if (status != CUOPT_SUCCESS) {
    pi.clear();
  }
  return pi;
}
ArrayRef<double> CuoptlpBackend::DualSolution_QP() {
  std::vector<double> pi(1);
  pi.clear();
  return pi;
}

double CuoptlpBackend::ObjectiveValue() const {
  ProblemData* problem_data = lp();
  cuopt_float_t objective_value;
  cuopt_int_t status = cuOptGetObjectiveValue(problem_data->solution, &objective_value);
  if (status != CUOPT_SUCCESS) {
    fmt::print("Error getting objective value\n");
  }
  return objective_value;
}

double CuoptlpBackend::NodeCount() const {
  return 0;
//  return getIntAttr(CUOPTLP_INTATTR_NODECNT);
}

double CuoptlpBackend::SimplexIterations() const {
  return 0;
//  return getIntAttr(CUOPTLP_INTATTR_SIMPLEXITER);
}

int CuoptlpBackend::BarrierIterations() const {
  return 0;
//  return getIntAttr(CUOPTLP_INTATTR_BARRIERITER);
}


void CuoptlpBackend::SetInterrupter(mp::Interrupter *inter) {
 // inter->SetHandler(InterruptCuoptlp, lp());
  // TODO Check interrupter
  //CUOPTLP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CuoptlpBackend::Solve() {
  ProblemData* problem_data = lp();

  cuopt_int_t is_mip;
  cuopt_int_t status = cuOptIsMIP(problem_data->problem, &is_mip);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Failed to check if problem is MIP, error code {}.", status));
  }
  //status = cuOptCreateSolverSettings(&problem_data->settings);
  //if (status != CUOPT_SUCCESS) {
  //  throw std::runtime_error(fmt::format("Failed to create solver settings, error code {}.", status));
  //}
  status = cuOptSolve(problem_data->problem, problem_data->settings, &problem_data->solution);
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Failed to solve problem, error code {}.", status));
  }


  WindupCUOPTLPSolve();
}

void CuoptlpBackend::WindupCUOPTLPSolve() {
}

void CuoptlpBackend::ReportResults() {
  ReportCUOPTLPResults();
  BaseBackend::ReportResults();
}

void CuoptlpBackend::ReportCUOPTLPResults() {
  SetStatus( GetSolveResult() );
  AddCUOPTLPMessages();
  if (need_multiple_solutions())
    ReportCUOPTLPPool();
}
std::vector<double> CuoptlpBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // CUOPTLP_CCALL(CUOPTLP_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double CuoptlpBackend::getPoolObjective(int i)
{
  double obj;
 // CUOPTLP_CCALL(CUOPTLP_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void CuoptlpBackend::ReportCUOPTLPPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(CUOPTLP_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void CuoptlpBackend::AddCUOPTLPMessages() {
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

std::pair<int, std::string> CuoptlpBackend::GetSolveResult() {
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


void CuoptlpBackend::FinishOptionParsing() {
  int v=-1;
 // GetSolverOption(CUOPTLP_INTPARAM_LOGGING, v);
  set_verbose_mode(v>0);

  // Nartive params
  if (storedOptions_.paramread_.size()) {
    //GRB_CALL(
    //  GRBreadparams(GRBgetenv(model()),
    //    paramfile_read().c_str()));
  }
  /// Set advanced parameters
  for (const auto& prm : storedOptions_.inlineparams_)
    this->SetSolverOption("Dummy", prm);
  // Write native params
  if (storedOptions_.paramwrite_.size()) {
    //GRB_CALL(
    //  GRBwriteparams(GRBgetenv(model()),
    //    paramfile_write().c_str()));
  }
  //lp()->SetVerbosity(storedOptions_.verbosity_);
}


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


void CuoptlpBackend::InitCustomOptions() {

  set_option_header(
      "CUOPTLP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cuoptlp_options``. For example::\n"
      "\n"
      "  ampl: option cuoptlp_options 'mipgap=1e-6';\n");

  // Use AddSolverOption() for proper solver parameters.
  // Below are examples of options stored in variables for own use.

  AddSolverOption("lim:timelimit timelimit",
    "Time limit in seconds after which the solver will stop and return the current solution",
    CUOPT_TIME_LIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:ncputhreads ncputhreads",
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

  AddSolverOption("tech:outlevel outlevel",
    "Set the output level for the solver",
    CUOPT_LOG_TO_CONSOLE,
    values_bool, "true");

  AddStoredOption("tech:logfile logfile",
    "Log file name.", storedOptions_.logFile_);

  ////////////////// CUSTOM RESULT CODES ///////////////////
  AddSolveResults( {
                     { sol::FAILURE+1, "fatal error 1" },
                     { sol::FAILURE+2, "fatal error 2" },
                     { sol::LIMIT_FEAS_NEW + 1, "AI iteration limit, feasible solution" },
                     { sol::LIMIT_NO_FEAS_NEW + 1, "AI iteration limit, no feasible solution" }
                   } );     // No replacement, make sure they are new
}


double CuoptlpBackend::MIPGap() {
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

double CuoptlpBackend::BestDualBound() {
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

double CuoptlpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}

#if 0
ArrayRef<int> CuoptlpBackend::VarStatii() {

  std::vector<int> vars(NumVars());
  /*
  if (!CUOPTLP_GetBasis(lp(), vars.data(), NULL))
    vars.clear();         // return empty if no basis
  for (auto& s : vars) {
    switch (s) {
    case CUOPTLP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case CUOPTLP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case CUOPTLP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case CUOPTLP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case CUOPTLP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Cuoptlp VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> CuoptlpBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  if (!CUOPTLP_GetBasis(lp(), NULL, cons.data()))
    cons.clear();          // return empty if no basis
  for (auto& s : cons) {
    switch (s) {
    case CUOPTLP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case CUOPTLP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case CUOPTLP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case CUOPTLP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case CUOPTLP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Cuoptlp VBasis value: {}", s));
    }
  }*/
  return cons;
}

void CuoptlpBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = CUOPTLP_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = CUOPTLP_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = CUOPTLP_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = CUOPTLP_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = CUOPTLP_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!CUOPTLP_GetColInfo(lp(), CUOPTLP_DBLINFO_LB, 1, index, &lb) &&
        !CUOPTLP_GetColInfo(lp(), CUOPTLP_DBLINFO_UB, 1, index, &ub))
      {
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
  CUOPTLP_SetBasis(lp(), stt.data(), NULL);
  */
}

void CuoptlpBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = CUOPTLP_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    //
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = CUOPTLP_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  CUOPTLP_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis CuoptlpBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());         // not for constraints, can be QCP
  }
  return { std::move(varstt), std::move(constt) };
}

void CuoptlpBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

void CuoptlpBackend::AddMIPStart(
    ArrayRef<double> x0, ArrayRef<int> sparsity) {
  //CUOPTLP_CCALL(CUOPTLP_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));
}
#endif

} // namespace mp


// AMPLs
void* AMPLSOpenCuoptlp(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::CuoptlpBackend()},
        cb);
}

void AMPLSCloseCuoptlp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetCuoptlpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CuoptlpBackend*>(AMPLSGetBackend(slv))->lp();
}

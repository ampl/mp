#include <vector>
#include <fstream>
#include <stdexcept>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "knitrompbackend.h"

extern "C" {
  #include "knitromp-ampls-c-api.h"    // Knitromp AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptKnitromp(void* prob) {
  //return KNITROMP_Interrupt((knitromp_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateKnitrompBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::KnitrompBackend()};
}


namespace mp {

/// Create Knitromp Model Manager
/// @param gc: the Knitromp common handle
/// @param e: environment
/// @param pre: presolver to be names[Solver::returned]= "";
/// need it to convert solution data
/// @return KnitrompModelMgr
std::unique_ptr<BasicModelManager>
CreateKnitrompModelMgr(KnitrompCommon&, Env&, pre::BasicValuePresolver*&);


int  callbackNewPoint(KN_context_ptr        kc,
    const double* const  x,
    const double* const  lambda,
    void* userParams)
{
    int i;
    int n;
    int nFC;
    int error;
    /** Demonstrate user-defined termination */
    /** (Uncomment to activate) */
    /*
    double dObj;
    */
    double dFeasError;
    if (x == nullptr) return 0;
    /** Get the number of variables in the model */
    error = KN_get_number_vars(kc, &n);

    printf("\n>> New point computed by Knitro: (");
    for (i = 0; i < (n - 1); i++)
        printf("%20.12e, ", x[i]);
    printf("%20.12e)\n", x[n - 1]);

    /** Query information about the current problem. */
    error = KN_get_number_FC_evals(kc, &nFC);
    printf("Number FC evals=%d, ", nFC);
    error = KN_get_abs_feas_error(kc, &dFeasError);
    printf("Current feasError=%e\n", dFeasError);

    /** Demonstrate user-defined termination */
    /** (Uncomment to activate) */
    /*
    error = KN_get_obj_value(kc, &dObj);
    if (dObj > 0.2 && dFeasError <= 1.0e-4)
    {
        return(  KN_RC_USER_TERMINATION );
    }
    */

    return(0);
}


KnitrompBackend::KnitrompBackend() {
  

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateKnitrompModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

KnitrompBackend::~KnitrompBackend() {
  CloseSolver();
}

void KnitrompBackend::OpenSolver() {
    const auto create_fn = GetCallbacks().init;
    if (create_fn)
        set_lp((KN_context_ptr)create_fn());
    else
    {
        KN_context_ptr prob;
        KNITROMP_CCALL(KN_new(&prob));
        set_lp(prob); // Assign it
    }
    
}


void KnitrompBackend::CloseSolver() {

  if ( lp() != NULL ) {
      KN_free(&lp_ref());
  }
}

void KnitrompBackend::InitOptionParsing()
{
    OpenSolver();
}


const char* KnitrompBackend::GetBackendName()
  { return "KnitrompBackend"; }

std::string KnitrompBackend::GetSolverVersion() {
  // TODO Return version from solver API
  return "0.0.0";
  //return fmt::format("{}.{}.{}", KNITROMP_VERSION_MAJOR, 
  //  KNITROMP_VERSION_MINOR, KNITROMP_VERSION_TECHNICAL);
}


bool KnitrompBackend::IsMIP() const {
    return false;
  // TODO. Use most precise information
  // (nonconvexities etc.)
  //return getIntAttr(Solver::NVARS_INT) > 0;
  //return getIntAttr(KNITROMP_INTATTR_ISMIP);
}

bool KnitrompBackend::IsQCP() const {
    return false;
  //return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_QUAD) > 0;
// return getIntAttr(KNITROMP_INTATTR_QELEMS) > 0;
}

ArrayRef<double> KnitrompBackend::PrimalSolution() {
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  int sstatus;
  double obj;
  int error = KN_get_solution(lp(), &sstatus, &obj, x.data(), NULL);
  if (error)
    x.clear();
  return x;
}

pre::ValueMapDbl KnitrompBackend::DualSolution() {
//  return { {
//    { CG_Linear, DualSolution_LP() },
//    { CG_Quadratic, DualSolution_QP() } } };
    return {};
}


double KnitrompBackend::ObjectiveValue() const {
    double obj;
    KNITROMP_CCALL(KN_get_obj_value(lp(), &obj));
    return obj;
}

double KnitrompBackend::NodeCount() const {
  return 0;
//  return getIntAttr(KNITROMP_INTATTR_NODECNT);
}

double KnitrompBackend::SimplexIterations() const {
  return 0;
//  return getIntAttr(KNITROMP_INTATTR_SIMPLEXITER);
}

int KnitrompBackend::BarrierIterations() const {
  return 0;
//  return getIntAttr(KNITROMP_INTATTR_BARRIERITER);
}


void KnitrompBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptKnitromp, lp());
}

void KnitrompBackend::Solve() {

      solstatus_ = KN_solve(lp());
      WindupKNITROMPSolve();
}

void KnitrompBackend::WindupKNITROMPSolve() { 
}

void KnitrompBackend::ReportResults() {
  ReportKNITROMPResults();
  BaseBackend::ReportResults();
}

void KnitrompBackend::ReportKNITROMPResults() {
  SetStatus( GetSolveResult() );
  AddKNITROMPMessages();
}

void KnitrompBackend::AddKNITROMPMessages() {
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

void KnitrompBackend::DoWriteProblem(const std::string& name) {

    auto* common = get_other();
    fmt::MemoryWriter w;
	common->call_format_model(w);

    std::ofstream outfile(name, std::ios::out | std::ios::trunc);
    if (!outfile.is_open()) {
        throw std::runtime_error(
            fmt::format("Cannot open file '{}' for writing", name));
    }
    outfile.write(w.data(), w.size());
    if (!outfile.good()) {
        throw std::runtime_error(
            fmt::format("Error writing to file '{}'", name));
    }

    fmt::print("Model exported to: {}\n", name);
}

// TODO Populate mp solver msgs
typedef struct { char* msg; int code, wantsol; } Sol_info;
static Sol_info solinfo[] = {
  { /* 0  */  "Locally optimal or satisfactory solution.", 0, 1 },
  { /*-100*/  "Current feasible solution estimate cannot be improved. Nearly optimal.", 100, 1 },
  { /*-101*/  "Relative change in feasible solution estimate < xtol for xtol_iters.", 101, 1 },
  { /*-102*/  "Current feasible solution estimate cannot be improved.", 102, 1 },
  { /*-103*/  "Relative change in feasible objective < ftol for ftol_iters.", 103, 1 },
  { /*-104*/  "Returning best feasible iterate.", 104, 1 },
  { /*-105*/  "Multistart: Feasible point found.", 105, 1 },
  { /*-200*/  "Convergence to an infeasible point. Problem may be locally infeasible.", 200, 1 },
  { /*-201*/  "Relative change in infeasible solution estimate < xtol for xtol_iters.", 201, 1 },
  { /*-202*/  "Current infeasible solution estimate cannot be improved.", 202, 1 },
  { /*-203*/  "Multistart: No primal feasible point found.", 203, 1 },
  { /*-204*/  "Problem determined to be infeasible with respect to constraint bounds.", 204, 1 },
  { /*-205*/  "Problem determined to be infeasible with respect to variable bounds.", 205, 1 },
  { /*-300*/  "Problem appears to be unbounded.", 300, 1 },
  { /*-301*/  "Problem is unbounded or infeasible.", 301, 1 },
  { /*-400*/  "Iteration limit reached. Current point is feasible.", 400, 1 },
  { /*-401*/  "Time limit reached. Current point is feasible.", 401, 1 },
  { /*-402*/  "Function evaluation limit reached. Current point is feasible.", 402, 1 },
  { /*-403*/  "MIP: All nodes have been explored. Integer feasible point found.", 403, 1 },
  { /*-404*/  "MIP: Integer feasible point found.", 404, 1 },
  { /*-405*/  "MIP: Subproblem solve limit reached. Integer feasible point found.", 405, 1 },
  { /*-406*/  "MIP: Node limit reached. Integer feasible point found.", 406, 1 },
  { /*-410*/  "Iteration limit reached. Current point is infeasible.", 410, 1 },
  { /*-411*/  "Time limit reached. Current point is infeasible.", 411, 1 },
  { /*-412*/  "Function evaluation limit reached. Current point is infeasible.", 412, 1 },
  { /*-413*/  "MIP: All nodes have been explored. No integer feasible point found.", 413, 1 },
  { /*-415*/  "MIP: Subproblem solve limit reached. No integer feasible point found.", 415, 1 },
  { /*-416*/  "MIP: Node limit reached. No integer feasible point found.", 416, 1 },
  { /*-501*/  "LP solver error.", 501, 1 },
  { /*-502*/  "Evaluation error.", 502, 1 },
  { /*-503*/  "Not enough memory.", 503, 1 },
  { /*-504*/  "Terminated by user.", 504, 1 },
  { /*-524*/  "Terminated after derivative check.", 505, 1 },
  { /*-505:-599*/ "Input or other API error.", 506, 0 },
  { /*-600*/  "Internal Knitro error.", 507, 0 },
  { /*    */  "Unknown termination.", 508, 0 },
  { /*    */  "Illegal objno value.", 509, 0 }
};

static int  setSolutionMessage(const int info)
{
    int infonum;

    switch (info) {
    case KN_RC_OPTIMAL_OR_SATISFACTORY: infonum = 0;    break;
    case KN_RC_NEAR_OPT:                infonum = 1;    break;
    case KN_RC_FEAS_XTOL:               infonum = 2;    break;
    case KN_RC_FEAS_NO_IMPROVE:         infonum = 3;    break;
    case KN_RC_FEAS_FTOL:               infonum = 4;    break;
    case KN_RC_FEAS_BEST:               infonum = 5;    break;
    case KN_RC_FEAS_MULTISTART:         infonum = 6;    break;
    case KN_RC_INFEASIBLE:              infonum = 7;    break;
    case KN_RC_INFEAS_XTOL:             infonum = 8;    break;
    case KN_RC_INFEAS_NO_IMPROVE:       infonum = 9;    break;
    case KN_RC_INFEAS_MULTISTART:       infonum = 10;    break;
    case KN_RC_INFEAS_CON_BOUNDS:       infonum = 11;    break;
    case KN_RC_INFEAS_VAR_BOUNDS:       infonum = 12;   break;
    case KN_RC_UNBOUNDED:               infonum = 13;   break;
    case KN_RC_UNBOUNDED_OR_INFEAS:     infonum = 14;   break;
    case KN_RC_ITER_LIMIT_FEAS:         infonum = 15;   break;
    case KN_RC_TIME_LIMIT_FEAS:         infonum = 16;   break;
    case KN_RC_FEVAL_LIMIT_FEAS:        infonum = 17;   break;
    case KN_RC_MIP_EXH_FEAS:            infonum = 18;   break;
    case KN_RC_MIP_TERM_FEAS:           infonum = 19;   break;
    case KN_RC_MIP_SOLVE_LIMIT_FEAS:    infonum = 20;   break;
    case KN_RC_MIP_NODE_LIMIT_FEAS:     infonum = 21;   break;
    case KN_RC_ITER_LIMIT_INFEAS:       infonum = 22;   break;
    case KN_RC_TIME_LIMIT_INFEAS:       infonum = 23;   break;
    case KN_RC_FEVAL_LIMIT_INFEAS:      infonum = 24;   break;
    case KN_RC_MIP_EXH_INFEAS:          infonum = 25;   break;
    case KN_RC_MIP_SOLVE_LIMIT_INFEAS:  infonum = 26;   break;
    case KN_RC_MIP_NODE_LIMIT_INFEAS:   infonum = 27;   break;
    case KN_RC_LP_SOLVER_ERR:           infonum = 28;   break;
    case KN_RC_EVAL_ERR:                infonum = 29;   break;
    case KN_RC_OUT_OF_MEMORY:           infonum = 30;   break;
    case KN_RC_USER_TERMINATION:        infonum = 31;   break;
    case KN_RC_DERIV_CHECK_TERMINATE:   infonum = 32;   break;
    case KN_RC_INTERNAL_ERROR:          infonum = 33;   break;
    default:
        infonum = info >= -599 && info <= -505 ? 33 : 35;
    }
    return infonum;
}

std::pair<int, std::string> KnitrompBackend::GetSolveResult() {
  namespace sol = mp::sol;
  
  int infonum = setSolutionMessage(solstatus_);
  int solve_result_num = solinfo[infonum].code;
  return { infonum, solinfo[infonum].msg};
}


void KnitrompBackend::FinishOptionParsing() {
  this->SetSolverOption(KN_PARAM_OUTLEV, storedOptions_.outlev);
  set_verbose_mode(storedOptions_.outlev >0);
  if(storedOptions_.outlev > 5)
     KNITROMP_CCALL(KN_set_newpt_callback(lp(), callbackNewPoint, NULL));

  useHessian = storedOptions_.hessian == 1;
  useJacobian = storedOptions_.jacobian == 1;
  // Override to enable recording of information even when we're writing out
  if (exportFileMode() > 0)
      printProblem = -exportFileMode();
  // Then override to enable actualy printing; writing is taken care of by
  // the framework, so we don't need to worry about it here.
  // So we record the model if printProblem < 0 (for writing)
  // and if ==1 (for printing)
  if (storedOptions_.printProblem == 1)
      printProblem = 1;
  copy_common_info_to_other();
}



void KnitrompBackend::AddPrimalDualStart(Solution sol0_unpres) {
    auto mv = GetValuePresolver().PresolveSolution(
        { sol0_unpres.primal, sol0_unpres.dual });
    auto& x0 = mv.GetVarValues()();
    auto& pi0 = mv.GetConValues()(CG_Linear);

    int status;
    
    status = KN_set_var_primal_init_values_all(lp(),
        x0.data());

    if (status)
        Print("warmstart: solution is not loaded "
            "because the problem is in presolved status.\n");
}


////////////////////////////// OPTIONS /////////////////////////////////



static const mp::OptionValueInfo outlev_values_[] = {
    {"0", "Nothing", 0},
    {"1", "Only final summary information", 1},
    {"2", "Information every 10 iterations is printed", 2},
    {"3", "Information at each iteration is printed", 3},
    {"4", "More verbose information at each iteration is printed", 4},
    {"5", "In addition, values of solution vector(x) are printed", 5},
    {"6", "In addition, constraints(c) and multipliers(lambda)", 6}
};


void KnitrompBackend::InitCustomOptions() {

  set_option_header(
      "KnitroMP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``knitro_options``. For example::\n"
      "\n"
      "  ampl: option knitro_options 'mipgap=1e-6';\n");

      AddSolverOption("alg:multistart multistart",
          "Use multistart (default 0)", 
          KN_PARAM_MULTISTART, KN_MULTISTART_NO, KN_MULTISTART_YES);
  
      AddStoredOption("tech:outlev outlev",
          "Specifies the verbosity of output:\n"
          "\n.. value-table::\n", storedOptions_.outlev,
          outlev_values_);

      AddStoredOption("tech:print print printproblem",
          "Print problem to screen before solving (default 0)",
          storedOptions_.printProblem);

      AddStoredOption("tech:threads threads",
          "Specifies the number of threads to use",
          storedOptions_.threads);


      AddStoredOption("alg:jacobian calculatejacobian jacobian",
          "Use calculated jacobian/gradient instead of relying on Knitro's (default 0)",
          storedOptions_.jacobian);

      AddStoredOption("alg:hessian calculatehessian hessian",
          "Use calculated hessian instead of relying on Knitro's (default 0)",
          storedOptions_.hessian);

      AddSolverOption("lim:time timelim timelimit",
          "limit on solve time (in seconds; default: no limit).",
          KN_PARAM_MAXTIME, 0.0, DBL_MAX);

  ////////////////// CUSTOM RESULT CODES ///////////////////



  AddSolveResults( {
                     //{sol::SOLVED,  "Locally optimal or satisfactory solution."},
                    // {sol::UNCERTAIN, "Current feasible solution estimate cannot be improved. Nearly optimal."},
                     {sol::UNCERTAIN+1, "Relative change in feasible solution estimate < xtol for xtol_iters."},
                     {sol::UNCERTAIN+2, "Current feasible solution estimate cannot be improved."},

                     
                   } );     // No replacement, make sure they are new
}



} // namespace mp


// AMPLs
void* AMPLSOpenKnitromp(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::KnitrompBackend()},
        cb);
}

void AMPLSCloseKnitromp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetKnitrompmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::KnitrompBackend*>(AMPLSGetBackend(slv))->lp();
}

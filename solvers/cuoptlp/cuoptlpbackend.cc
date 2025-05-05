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
  int status = 0;
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
  Solver::SolverModel* prob = Solver::CreateSolverModel();
  set_lp(prob); // Assign it
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  /* TODO Typically check call */
  /// Turn off verbosity by default
  // CUOPTLP_CCALL(CUOPTLP_SetIntParam(prob, "Logging", 0));

}

void CuoptlpBackend::CloseSolver() {
  /* TODO Cleanup: close problem and environment
  if ( lp() != NULL ) {
    CUOPTLP_CCALL(CUOPTLP_DeleteProb(&lp_) );
  }
  if ( env() != NULL ) {
    CUOPTLP_CCALL(CUOPTLP_DeleteEnv(&env_) );
  }
  */
}

const char* CuoptlpBackend::GetBackendName()
  { return "CuoptlpBackend"; }

std::string CuoptlpBackend::GetSolverVersion() {
  // TODO Return version from solver API
  int32_t major, minor, patch;
  cuOptGetSemanticVersion(&major, &minor, &patch);
  return fmt::format("{}.{}.{}", major, minor, patch);
  //return fmt::format("{}.{}.{}", CUOPTLP_VERSION_MAJOR,
  //  CUOPTLP_VERSION_MINOR, CUOPTLP_VERSION_TECHNICAL);
}


bool CuoptlpBackend::IsMIP() const {
  // TODO. Use most precise information
  // (nonconvexities etc.)
  return getIntAttr(Solver::NVARS_INT) > 0;
  //return getIntAttr(CUOPTLP_INTATTR_ISMIP);
}

bool CuoptlpBackend::IsQCP() const {
  return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_QUAD) > 0;
// return getIntAttr(CUOPTLP_INTATTR_QELEMS) > 0;
}

ArrayRef<double> CuoptlpBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error = 1;
  std::vector<double> x(num_vars);
  // We should always return a solution when available,
  // even if infeasible/suboptimal etc.
  // User decides on it using solve_result.
  /*
  if (IsMIP())
    error = CUOPTLP_GetSolution(lp(), x.data());
  else
    error = CUOPTLP_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);
    */
  if (error)
    x.clear();
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
 // int error = CUOPTLP_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 1;
  if (error)
    pi.clear();
  return pi;
}
ArrayRef<double> CuoptlpBackend::DualSolution_QP() {
  int num_cons = NumQPCons();
  std::vector<double> pi(num_cons);
  // int error = CUOPTLP_GetQpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 1;
  if (error)
    pi.clear();
  return pi;
}

double CuoptlpBackend::ObjectiveValue() const {
 /* if (IsMIP())
    return getDblAttr(CUOPTLP_DBLATTR_BESTOBJ);
  else
    return getDblAttr(CUOPTLP_DBLATTR_LPOBJVAL);
    */
  return 0;
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
  inter->SetHandler(InterruptCuoptlp, lp());
  // TODO Check interrupter
  //CUOPTLP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CuoptlpBackend::Solve() {
  //CUOPTLP_CCALL(CUOPTLP_Solve(lp()));
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

void CuoptlpBackend::printModelStats() {

  std::map<Solver::ConsType, std::string> names;
    names[Solver::CONS_LIN] = "Linear";
    names[Solver::CONS_QUAD]= "Quadratic";
    names[Solver::CONS_QUAD_CONE]= "Cone quadratic";
    names[Solver::CONS_QUAD_CONE_ROTATED]= "Cone rotated";
    names[Solver::CONS_QUAD_CONE_EXP] = "Cone exponential";
    names[Solver::CONS_INDIC]= "Indicator";
    names[Solver::CONS_SOS]= "SOS";

    names[Solver::CONS_MAX]= "Max";
    names[Solver::CONS_MIN]= "Min";
    names[Solver::CONS_ABS]= "Abs";
    names[Solver::CONS_AND]= "And";
    names[Solver::CONS_OR]= "Or";

    names[Solver::CONS_EXP]= "Exp";
    names[Solver::CONS_EXPA]= "ExpA";
    names[Solver::CONS_LOG]= "Log";
    names[Solver::CONS_LOGA]= "LogA";

    names[Solver::CONS_POW]= "Pow";
    names[Solver::CONS_SIN]= "Sin";
    names[Solver::CONS_COS]= "Cos";
    names[Solver::CONS_TAN]= "Tan";
    names[Solver::CONS_PL] = "Piecewise linear";

    AddToSolverMessage("\n\n##### Model stats #####\n");

    AddToSolverMessage(fmt::format("Variables: ({})\n", getIntAttr(Solver::NVARS)));
    AddToSolverMessage(fmt::format("  {} continuous\n", getIntAttr(Solver::NVARS_CONT)));
    AddToSolverMessage(fmt::format("  {} integer\n", getIntAttr(Solver::NVARS_INT)));
    AddToSolverMessage(fmt::format("  {} binary\n", getIntAttr(Solver::NVARS_BIN)));


    AddToSolverMessage(fmt::format("\nObjectives: ({})\n", getIntAttr(Solver::NOBJS)));
    AddToSolverMessage(fmt::format("  {} linear\n", getIntAttr(Solver::NOBJS, Solver::CONS_LIN)));
    AddToSolverMessage(fmt::format("  {} quadratic\n", getIntAttr(Solver::NOBJS, Solver::CONS_QUAD)));
    AddToSolverMessage(fmt::format("  {} non-linear\n", getIntAttr(Solver::NOBJS, Solver::CONS_NL)));


    AddToSolverMessage(fmt::format("\nConstraints: ({})\n", getIntAttr(Solver::NCONS)));
    for (const auto& i : names) {
      auto n = getIntAttr(Solver::NCONS_TYPE, i.first);
      if (n == 0) continue;
      AddToSolverMessage(fmt::format("  {} {}\n", n, i.second));
    }
}
void CuoptlpBackend::AddCUOPTLPMessages() {
  printModelStats();
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
  if (IsMIP())
  {
  //
  }
  else {
  //
  }
  return { sol::UNKNOWN, "not solved" };
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
  lp()->SetVerbosity(storedOptions_.verbosity_);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo verbosity_values_[] = {
  { "0", "Only statistics", 0},
  { "1", "All info", 1}
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

  AddStoredOption("tech:option_example opt_example example_opt",
      "Example option. "
      "Default = \"\" (don't work too hard).",
      storedOptions_.option_example_);

  AddStoredOption("tech:flag1 flag1",
      "Flag option. Use without value. Can only be set to True.",
      storedOptions_.flag_option_);

  AddListOption("tech:list_option opt_list multi_valued_option",
      "Multi-valued option when repeated.",
      storedOptions_.list_option_);

  // Native solver options handling.
  // Actual processing of these options can be done in FinishOptionParsing().
  AddListOption("tech:optionnative optionnative optnative tech:param",
      "General way to specify values of both documented and "
      "undocumented Gurobi parameters; value should be a quoted "
      "string (delimited by ' or \") containing a parameter name, a "
      "space, and the value to be assigned to the parameter.  Can "
      "appear more than once.  Cannot be used to query current "
      "parameter values.",
      storedOptions_.inlineparams_);
  AddStoredOption("tech:optionnativeread tech:param:read param:read optnative:read",
      "Name of Gurobi parameter file (surrounded by 'single' or "
      "\"double\" quotes if the name contains blanks). "
      "The suffix on a parameter file should be .prm, optionally followed "
      "by .zip, .gz, .bz2, or .7z.\n"
      "\n"
      "Lines that start with # are ignored.  Otherwise, each nonempty "
      "line should contain a name and a value, separated by a space.",
      storedOptions_.paramread_);
  AddStoredOption("tech:optionnativewrite tech:param:write param:write optnative:write",
      "Name of Gurobi parameter file (surrounded by 'single' or \"double\" quotes if the "
      "name contains blanks) to be written.",
      storedOptions_.paramwrite_);

  // Example of an option with possible values defined in a table. Look at other drivers
  // for predefined value tables (eg. values_autonoyes_, ...)
  AddStoredOption("tech:verbosity verbosity",
    "Set the verbosity of this run, deciding what to print to console:\n"
    "\n.. value-table::\n", storedOptions_.verbosity_, verbosity_values_);

  ////////////////// CUSTOM RESULT CODES ///////////////////
  AddSolveResults( {
                     { sol::FAILURE+1, "fatal error 1" },
                     { sol::FAILURE+2, "fatal error 2" },
                     { sol::LIMIT_FEAS_NEW + 1, "AI iteration limit, feasible solution" },
                     { sol::LIMIT_NO_FEAS_NEW + 1, "AI iteration limit, no feasible solution" }
                   } );     // No replacement, make sure they are new
}


double CuoptlpBackend::MIPGap() {
  return 0;
//  return getDblAttr(CUOPTLP_DBLATTR_BESTGAP);
}
double CuoptlpBackend::BestDualBound() {
  return 0;
  //return getDblAttr(CUOPTLP_DBLATTR_BESTBND);
}

double CuoptlpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


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


void CuoptlpBackend::ComputeIIS() {
  //CUOPTLP_CCALL(CUOPTLP_ComputeIIS(lp()));
  SetStatus(GetSolveResult());   // could be new information
}

IIS CuoptlpBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> CuoptlpBackend::VarsIIS() {
  return ArrayRef<int>();
//  return getIIS(lp(), NumVars(), CUOPTLP_GetColLowerIIS, CUOPTLP_GetColUpperIIS);
}
pre::ValueMapInt CuoptlpBackend::ConsIIS() {
  /*auto iis_lincon = getIIS(lp(), NumLinCons(), CUOPTLP_GetRowLowerIIS, CUOPTLP_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  CUOPTLP_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  CUOPTLP_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };
      */
  return { {{ 0, std::vector<int>()}} };
}

void CuoptlpBackend::AddMIPStart(
    ArrayRef<double> x0, ArrayRef<int> sparsity) {
  //CUOPTLP_CCALL(CUOPTLP_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));
}


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

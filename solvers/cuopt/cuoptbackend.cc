#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cuoptbackend.h"

extern "C" {
  #include "cuopt-ampls-c-api.h"    // Cuopt AMPLS C API
}
#include "mp/ampls-cpp-api.h"

#include <iostream>

namespace {


bool InterruptCuopt(void* prob) {
  //return CUOPT_Interrupt((cuopt_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCuoptBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CuoptBackend()};
}


namespace mp {

/// Create Cuopt Model Manager
/// @param gc: the Cuopt common handle
/// @param e: environment
/// @param pre: presolver to be names[Solver::returned]= "";
/// need it to convert solution data
/// @return CuoptModelMgr
std::unique_ptr<BasicModelManager>
CreateCuoptModelMgr(CuoptCommon&, Env&, pre::BasicValuePresolver*&);


CuoptBackend::CuoptBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateCuoptModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

CuoptBackend::~CuoptBackend() {
  CloseSolver();
}

void CuoptBackend::OpenSolver() {
  int status = 0;

  const std::string server_ip = "0.0.0.0";
  const int server_port = 5000;


  httplib::Client* client = new httplib::Client(server_ip, server_port);
  set_client(client);
  
  json* prob = new json;
  json* sol = new json;
  set_json_prob(prob);
  set_json_sol(sol);
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
}

void CuoptBackend::CloseSolver() {
  /* TODO Cleanup: close problem and environment
  if ( lp() != NULL ) {
    CUOPT_CCALL(CUOPT_DeleteProb(&lp_) );
  }
  if ( env() != NULL ) {
    CUOPT_CCALL(CUOPT_DeleteEnv(&env_) );
  }
  */
}

const char* CuoptBackend::GetBackendName()
{ return "CuoptBackend"; }

std::string CuoptBackend::GetSolverVersion() {
  return "24.11";
}


bool CuoptBackend::IsMIP() const {
  // TODO. Use most precise information
  // (nonconvexities etc.)
  return getIntAttr(Solver::NVARS_INT) > 0;
  //return getIntAttr(CUOPT_INTATTR_ISMIP);
}

bool CuoptBackend::IsQCP() const {
  return false;
}

ArrayRef<double> CuoptBackend::PrimalSolution() {
  json* sol = get_json_sol();
  std::vector<double> x = (*sol)["response"]["solver_response"]["solution"]["primal_solution"];
  return x;
}

pre::ValueMapDbl CuoptBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> CuoptBackend::DualSolution_LP() {
  json* sol = get_json_sol();
  std::vector<double> pi = (*sol)["response"]["solver_response"]["solution"]["dual_solution"];
  return pi;
}


double CuoptBackend::ObjectiveValue() const {
  json* sol = get_json_sol();

  return (*sol)["response"]["solver_response"]["solution"]["primal_objective"];
}

double CuoptBackend::NodeCount() const {
  return 0;
}

double CuoptBackend::SimplexIterations() const {
  return 0;
}

int CuoptBackend::BarrierIterations() const {
  return 0;
}


void CuoptBackend::SetInterrupter(mp::Interrupter *inter) {
  //inter->SetHandler(InterruptCuopt, lp());
  // TODO Check interrupter
  //CUOPT_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void CuoptBackend::Solve() {
  json* prob = get_json_prob();
  httplib::Client* client = get_client();
  httplib::Headers headers = {
    {"Content-Type", "application/json"},
    {"CLIENT-VERSION", "custom"}
  };

  // Send the POST request
  std::cout << "Sending request to cuopt server" << (*prob).dump(2) << std::endl;
  auto res = (*client).Post("/cuopt/request", headers, (*prob).dump(), "application/json");

  json response = json::parse(res->body);
  set_uuid(response["reqId"]);

  std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(10));

  auto res_sol = (*client).Get("/cuopt/request/" + get_uuid(), headers);
  json response_sol = json::parse(res_sol->body);
  json* copy = new json(response_sol);
  set_json_sol(copy);

  WindupCUOPTSolve();
}

void CuoptBackend::WindupCUOPTSolve() { 
}

void CuoptBackend::ReportResults() {
  ReportCUOPTResults();
  BaseBackend::ReportResults();
}

void CuoptBackend::ReportCUOPTResults() {
  SetStatus( GetSolveResult() );
  AddCUOPTMessages();
  if (need_multiple_solutions())
    ReportCUOPTPool();
}
std::vector<double> CuoptBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // CUOPT_CCALL(CUOPT_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double CuoptBackend::getPoolObjective(int i)
{
  double obj;
 // CUOPT_CCALL(CUOPT_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void CuoptBackend::ReportCUOPTPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(CUOPT_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}

void CuoptBackend::AddCUOPTMessages() {
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

std::pair<int, std::string> CuoptBackend::GetSolveResult() {
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
  json* solution = get_json_sol();

  if ((*solution)["response"]["solver_response"]["status"] == 1) {
    return { sol::SOLVED, "optimal solution" };
  }
}


void CuoptBackend::FinishOptionParsing() {
  int v=-1;
 // GetSolverOption(CUOPT_INTPARAM_LOGGING, v);
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



void CuoptBackend::InitCustomOptions() {

  set_option_header(
      "cuOpt Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cuopt_options``. For example::\n"
      "\n"
      "  ampl: option cuopt_options 'mipgap=1e-6';\n");

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

  AddStoredOption("tech:infeasibility_detection infeasibility_detection",
      "Detect and leave if the problem is detected as infeasible. Default = true.",
      storedOptions_.infeasibility_detection_);

  AddStoredOption("tech:solver_mode solver_mode",
      "Solver mode to set. Only possible values are 0, 1, and 2. Default = 0.",
      storedOptions_.solver_mode_);


  ////////////////// CUSTOM RESULT CODES ///////////////////
  AddSolveResults( {
                     { sol::FAILURE+1, "fatal error 1" },
                     { sol::FAILURE+2, "fatal error 2" },
                     { sol::LIMIT_FEAS_NEW + 1, "AI iteration limit, feasible solution" },
                     { sol::LIMIT_NO_FEAS_NEW + 1, "AI iteration limit, no feasible solution" }
                   } );     // No replacement, make sure they are new
}


double CuoptBackend::MIPGap() {
  return 0;
//  return getDblAttr(CUOPT_DBLATTR_BESTGAP);
}
double CuoptBackend::BestDualBound() {
  return 0;
  //return getDblAttr(CUOPT_DBLATTR_BESTBND);
}

double CuoptBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> CuoptBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  return vars;
}

ArrayRef<int> CuoptBackend::ConStatii() {
  std::vector<int> cons(NumLinCons());
  return cons;
}

void CuoptBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
}

void CuoptBackend::ConStatii(ArrayRef<int> cst) { }

SolutionBasis CuoptBackend::GetBasis() {
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

void CuoptBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

ArrayRef<int> CuoptBackend::VarsIIS() {
  return ArrayRef<int>();
}
pre::ValueMapInt CuoptBackend::ConsIIS() {
  return { {{ 0, std::vector<int>()}} };
}

void CuoptBackend::AddPrimalDualStart(Solution sol0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution(
        { sol0_unpres.primal, sol0_unpres.dual } );
  auto x0 = mv.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);
  json* prob = get_json_prob();

  (*prob)["initial_solution"]["primal"] = json::array();
  for (auto& value : x0) {
    (*prob)["initial_solution"]["primal"].push_back(value);
  }
  
  (*prob)["initial_solution"]["dual"] = json::array();
  for (auto& value : pi0) {
    (*prob)["initial_solution"]["dual"].push_back(value);
  }
}


} // namespace mp


// AMPLs
void* AMPLSOpenCuopt(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::CuoptBackend()},
        cb);
}

void AMPLSCloseCuopt(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetCuoptmodel(AMPLS_MP_Solver* slv) {
  return nullptr;
  //return dynamic_cast<mp::CuoptBackend*>(AMPLSGetBackend(slv))->lp();
}

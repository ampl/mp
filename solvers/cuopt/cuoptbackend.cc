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
  /*httplib::Client* client = get_client();
  std::string uuid = get_uuid();
  httplib::Headers headers = {
    {"Content-Type", "application/json"},
    {"CLIENT-VERSION", "custom"}
  };
  auto res_sol = (*client).Get("/cuopt/request/" + get_uuid(), headers);
  json response_sol = json::parse(res_sol->body);
  json* copy = new json(response_sol);
  set_json_sol(copy);

  auto res = (*client).Delete("/cuopt/request" + get_uuid(), headers);

  std::cout << "Interrupted!" << std::endl;
  return true;*/
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
  json *prob = get_json_prob();
  json *sol = get_json_sol();

  delete prob;
  delete sol;
}

const char* CuoptBackend::GetBackendName()
{ return "CuoptBackend"; }

std::string CuoptBackend::GetSolverVersion() {
  return "24.11";
}


bool CuoptBackend::IsMIP() const {
  return isMIP();
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

}

void CuoptBackend::Solve() {
  json* prob = get_json_prob();
  httplib::Client* client = get_client();
  httplib::Headers headers = {
    {"Content-Type", "application/json"},
    {"CLIENT-VERSION", "custom"}
  };


  //if (isMIP()) {
    (*prob)["solver_config"]["time_limit"] = storedOptions_.time_limit_;
  //}

  // Send the POST request
  //std::cout << "Sending request to cuopt server" << (*prob).dump(2) << std::endl; //debug
  auto res = (*client).Post("/cuopt/request", headers, (*prob).dump(), "application/json");

  json response = json::parse(res->body);
  set_uuid(response["reqId"]);

  json response_sol;

  if (isMIP()) {
    std::this_thread::sleep_until(std::chrono::system_clock::now() 
      + std::chrono::seconds(int(storedOptions_.time_limit_)) + std::chrono::seconds(5));
    auto res_sol = (*client).Get("/cuopt/request/" + get_uuid(), headers);
    response_sol = json::parse(res_sol->body);
  }
  else {
    do {
      std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(10));
      auto res_sol = (*client).Get("/cuopt/request/" + get_uuid(), headers);
      std::cout << "Current solution: " << response_sol["response"]["solver_response"]["solution"]["primal_solution"] << std::endl;
      response_sol = json::parse(res_sol->body);
    } while(response_sol["response"]["solver_response"]["status"] != 1);
  }
  
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
  else if ((*solution)["response"]["solver_response"]["status"] == 5) {
    return { sol::LIMIT_FEAS_TIME, "time limit, feasible solution" };
  }
}


void CuoptBackend::FinishOptionParsing() {
  json* prob = get_json_prob();
  int v=-1;
  set_verbose_mode(v>0);

  if (storedOptions_.time_limit_)
    (*prob)["solver_config"]["time_limit"] = storedOptions_.time_limit_;
  if (storedOptions_.iteration_limit_)
    (*prob)["solver_config"]["iteration_limit"] = storedOptions_.iteration_limit_;
  (*prob)["solver_config"]["infeasibility_detection"] = storedOptions_.infeasibility_detection_;
  (*prob)["solver_config"]["solver_mode"] = storedOptions_.solver_mode_;

  if (storedOptions_.optimality_)
    (*prob)["solver_config"]["tolerances"]["optimality"] = storedOptions_.optimality_;
  if (storedOptions_.absolute_primal_)
    (*prob)["solver_config"]["tolerances"]["absolute_primal"] = storedOptions_.absolute_primal_;
  if (storedOptions_.absolute_dual_)
    (*prob)["solver_config"]["tolerances"]["absolute_dual"] = storedOptions_.absolute_dual_;
  if (storedOptions_.absolute_gap_)
    (*prob)["solver_config"]["tolerances"]["absolute_gap"] = storedOptions_.absolute_gap_;
  if (storedOptions_.relative_primal_)
    (*prob)["solver_config"]["tolerances"]["relative_primal"] = storedOptions_.relative_primal_;
  if (storedOptions_.relative_dual_)
    (*prob)["solver_config"]["tolerances"]["relative_dual"] = storedOptions_.relative_dual_;
  if (storedOptions_.relative_gap_)
    (*prob)["solver_config"]["tolerances"]["relative_gap"] = storedOptions_.relative_gap_;
  if (storedOptions_.primal_infeasible_)
    (*prob)["solver_config"]["tolerances"]["primal_infeasible"] = storedOptions_.primal_infeasible_;
  if (storedOptions_.dual_infeasible_)
    (*prob)["solver_config"]["tolerances"]["dual_infeasible"] = storedOptions_.dual_infeasible_;
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

  AddStoredOption("tech:infeasibility_detection infeasibility_detection",
      "Detect and leave if the problem is detected as infeasible. Default = true.",
      storedOptions_.infeasibility_detection_);

  AddStoredOption("tech:solver_mode solver_mode",
      "Solver mode to set. Only possible values are 0, 1, and 2. Default = 0.",
      storedOptions_.solver_mode_);


  AddStoredOption("tol:optimality optimality",
      "Absolute and relative tolerance on the primal feasibility, dual feasibility, and gap. Default = 1e-4.",
      storedOptions_.optimality_);

  AddStoredOption("tol:absolute_primal absolute_primal",
      "Absolute primal tolerance. Default = 1e-4.",
      storedOptions_.absolute_primal_);

  AddStoredOption("tol:absolute_dual absolute_dual",
      "Absolute dual tolerance. Default = 1e-4.",
      storedOptions_.absolute_dual_);

  AddStoredOption("tol:absolute_gap absolute_gap",
      "Absolute gap tolerance. Default = 1e-4.",
      storedOptions_.absolute_gap_);

  AddStoredOption("tol:relative_primal relative_primal",
      "Relative primal tolerance. Default = 1e-4.",
      storedOptions_.relative_primal_);

  AddStoredOption("tol:relative_dual relative_dual",
      "Relative dual tolerance. Default = 1e-4.",
      storedOptions_.relative_dual_);

  AddStoredOption("tol:relative_gap relative_gap",
      "Relative gap tolerance. Default = 1e-4.",
      storedOptions_.relative_gap_);

  AddStoredOption("tol:primal_infeasible primal_infeasible",
      "Primal infeasible tolerance. Default = 1e-4.",
      storedOptions_.primal_infeasible_);

  AddStoredOption("tol:dual_infeasible dual_infeasible",
      "Dual infeasible tolerance. Default = 1e-4.",
      storedOptions_.dual_infeasible_);


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

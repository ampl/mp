#include <vector>
#include <climits>
#include <cfloat>
#include <iostream> // for cerr


#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "baronmpbackend.h"


#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include "ConsoleApi2.h"
#else
#include <csignal>  // for kill()
#endif


extern "C" {
  #include "baronmp-ampls-c-api.h"    // Baronmp AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptBaronmp(void* prob) {
  #ifndef WIN32
  kill(mp::BaronmpCommon::pid, SIGINT);
  #else
  // CTRL_C_EVENT = 0
  if (!GenerateConsoleCtrlEvent(0, mp::BaronmpCommon::pid)) {
        std::cerr << "GenerateConsoleCtrlEvent failed (" << GetLastError() << ").\n";
        return false;
    }
  #endif
  return true;
}

}  

std::unique_ptr<mp::BasicBackend> CreateBaronmpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::BaronmpBackend()};
}


namespace mp {

/// Create Baronmp Model Manager
/// @param gc: the Baronmp common handle
/// @param e: environment
/// @param pre: presolver to be names[Solver::returned]= "";
/// need it to convert solution data
/// @return BaronmpModelMgr
std::unique_ptr<BasicModelManager>
CreateBaronmpModelMgr(BaronmpCommon&, Env&, pre::BasicValuePresolver*&);


BaronmpBackend::BaronmpBackend() {
  OpenSolver();
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateBaronmpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

BaronmpBackend::~BaronmpBackend() {
  CloseSolver();
}

void BaronmpBackend::OpenSolver() {
  set_lp(new BaronProblemInfo());
  int status = 0;
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
}


void BaronmpBackend::CloseSolver() {
  if(!baronDir.empty())
  {
    if(baronOptions().keepsol)
      fmt::print("Keeping temporary files in {}\n", baronDir);
    else
    {
      deinitBaronFile();
      recrmdir(baronDir);

    }
  } 
}

const char* BaronmpBackend::GetBackendName()
  { return "BaronmpBackend"; }

std::string BaronmpBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", v_year, v_month, v_day);
}


bool BaronmpBackend::IsMIP() const {
  return false;
}

bool BaronmpBackend::IsQCP() const {
  return false;
}

std::map<std::string, std::variant<int, double, std::string>>
BaronmpBackend::SolutionStats() {
    std::map<std::string, std::variant<int, double, std::string>> stats;
    stats["barrier_iterations"] = timFileData_.itera;
    return stats;
}

ArrayRef<double> BaronmpBackend::PrimalSolution() {
  auto vec = resFileData_.PrimalSolution();
  int num_vars = vec.size();
  std::vector<double> x(num_vars);
  for (int i = 0; i < num_vars; i++)
    x[lp()->baronToAMPLIndices[i]] = vec[i];
  return x;
}

pre::ValueMapDbl BaronmpBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> BaronmpBackend::DualSolution_LP() {
  return resFileData_.DualSolution();
}

double BaronmpBackend::ObjectiveValue() const {
  return resFileData_.ObjectiveValue();
}

double BaronmpBackend::NodeCount() const {
  return 0;
}

double BaronmpBackend::SimplexIterations() const {
  return 0;
}

int BaronmpBackend::BarrierIterations() const {
  return timFileData_.itera;
}


void BaronmpBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptBaronmp, lp());
}

void BaronmpBackend::Solve() {
  deinitBaronFile(false);
  

  BaronGlobals g;
  g.baronDir = baronDir;
  g.initialDir = initialDir;
  g.lsolmsg = baronOptions().lsolmsg;
  if (!baronOptions().lsolver.empty())
    g.lsolver = baronOptions().lsolver;
  if (baronOptions().lpsolver == "cplex")
    g.lpsol_dll = baronOptions().cplexlibname;
  g.nlfile = nlFilePath;
  g.nvars = nVarsInteger() + nVarsContinuous()+ nVarsBinary();
  g.verbuf = version();
  g.v_keepsol = baronOptions().keepsol;
  g.serialize(filePathAMPL());
  int status =  BaronmpCommon::runBaron(filePathBar(), baronOptions().maxtime);
  if (status)
    errorLevel = -1; // Problems running process

  WindupBARONMPSolve();
}

void BaronmpBackend::WindupBARONMPSolve() {
  timFileData_ = TimFileData::fromFile(BARON_TIM);
  resFileData_ = ResFileData::fromFile(BARON_RES, timFileData_.nodeopt);

}

void BaronmpBackend::ReportResults() {
  ReportBARONMPResults();
  BaseBackend::ReportResults();
}

void BaronmpBackend::ReportBARONMPResults() {
  SetStatus( GetSolveResult() );
  AddBARONMPMessages();
  if (need_multiple_solutions())
    ReportBARONMPPool();
}

void BaronmpBackend::ReportBARONMPPool() {
  for(int i=0; i<resFileData_.NumSol();  i++)
  for(const auto &sol : resFileData_.PrimalSolutions()){
    ReportIntermediateSolution(
      { resFileData_.PrimalSolutions()[i],
        {}, {  resFileData_.Objectives()[i]}});
      }
}

void BaronmpBackend::AddBARONMPMessages() {
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} iterations\n", nbi));
}

std::pair<int, std::string> BaronmpBackend::GetSolveResult() {
  namespace sol = mp::sol;

  if(errorLevel == -1) // could not call Baron process
    return { sol::UNKNOWN, "not solved, could not call Baron process" };
  switch (timFileData_.modelstatus) {
  case 1:
    return { sol::SOLVED, "optimal solution" };
  case 2:
    return { sol::INFEASIBLE, "infeasible problem" };
  case 3:
    return { sol::UNBOUNDED_NO_FEAS, "unbounded, no feasible solution returned" };
  case 4:
    if (timFileData_.barstatus == 1)
      return { sol::UNCERTAIN, "numerical difficulties but possibly optimal" };
  default:
    ;
  }
  switch (timFileData_.barstatus)
  {
  case 2:
    return { sol::LIMIT_NO_FEAS_NODES, "node limit, without a feasible solution" };
  case 3:
    return { sol::LIMIT_NO_FEAS_ITER, "iteration limit, without a feasible solution" };
  case 4:
    return { sol::LIMIT_NO_FEAS_TIME, "time limit, without a feasible solution" };
  case 5:
    return { sol::NUMERIC, "terminated due to unrecoverable numerical issues" };
  case 6:
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted, no feasible solution" };
  case 7:
    return { sol::FAILURE+2, "too little memory" };
  case 9:
    return { sol::FAILURE+1, "BARON syntax error (should not happen)" };
  case 10:
    return { sol::FAILURE, "licensing error" };
  case 11:
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted by Control-C" };
  }
  return { sol::UNKNOWN, fmt::format("Unknown solution status: {}", timFileData_.barstatus) };
}


void BaronmpBackend::FinishOptionParsing() {
  initDirectories(stub(), baronOptions().scratch,
    baronOptions().overwrite);
  baronOptions().setProblemName(stub());
  copy_common_info_to_other();
  set_verbose_mode(baronOptions().outlev > 0);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo iismethod_values[] = {
  { "1", "Fast heuristic", 1},
  { "2", "Deletion filtering algorithm", 2},
  { "3", "Addition filtering algorithm", 3},
  { "4", "Addition-deletion filtering algorithm", 4},
  { "5", "Depth-first search algorithm", 5}
};


static const mp::OptionValueInfo iisorder_values[] = {
  { "-1", "Automatic", -1},
  { "1", "Problem order (as in reformulated model)", 1},
  { "2", "Ascending order by degree", 2},
  { "3", "Descending order by degree", 3},
  { ">=4", "Random order with seed iisorder", 4} 
}; 


void BaronmpBackend::InitCustomOptions() {

  set_option_header(
      "BARONMP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``baronmp_options``. For example::\n"
      "\n"
      "  ampl: option baronmp_options 'mipgap=1e-6';\n");



  #define ADDOPTION(var, n)\
  AddStoredOption(n, baronOptions().description(n), baronOptions().var)
  #define ADDOPTIONV(var, n, valuetable)\
  AddStoredOption(n, baronOptions().description(n), baronOptions().var, valuetable)

  #define ADDALGOPTION(n) ADDOPTION(n, "alg:"#n " "#n)
  #define ADDTECHOPTION(n) ADDOPTION(n,"tech:"#n " "#n)
  
  ADDALGOPTION(deltaa);
  ADDALGOPTION(deltaterm);
  ADDALGOPTION(deltar);
  ADDALGOPTION(deltat);
  ADDALGOPTION(epsa);
  ADDALGOPTION(epsr);
  ADDALGOPTION(firstfeas);
  ADDALGOPTION(firstloc);

  // ADDALGOPTION(iisint);
  // ADDALGOPTION(iismethod);
  // ADDALGOPTION(iisorder);

  ADDTECHOPTION(barstats);
  ADDTECHOPTION(keepsol);
  ADDTECHOPTION(lpsolver);
  ADDTECHOPTION(cplexlibname);
  ADDTECHOPTION(xpresslibname);
  ADDTECHOPTION(lsolmsg);
  ADDTECHOPTION(lsolver);

  ADDOPTION(maxiter, "lim:iter iterlimit maxiter");
  ADDOPTION(maxtime, "lim:time timelim timelimit maxtime");

  ADDOPTION(nlpsolver, "tech:nlpsolver nlpsol nlpsolver");
  
  ADDTECHOPTION(numsol);
  ADDTECHOPTION(objbound);
  ADDTECHOPTION(objno);
  ADDTECHOPTION(optfile);

  
  ADDTECHOPTION(outlev);
  ADDTECHOPTION(overwrite); 
  ADDTECHOPTION(prfreq);
  ADDOPTIONV(prloc, "tech:prloc prloc", values_01_noyes_0default_);
  ADDTECHOPTION(problem);
  ADDTECHOPTION(prtime);
  ADDTECHOPTION(scratch);
  ADDTECHOPTION(seed);
  ADDTECHOPTION(sumfile);
  ADDTECHOPTION(threads);
  ADDTECHOPTION(trace);

  AddListOption("tech:optionnative optionnative optnative tech:param",
    "General way to specify values of both documented and "
    "undocumented Baron parameters; value should be a quoted "
    "string (delimited by ' or \") containing a parameter name, a "
    "space, and the value to be assigned to the parameter.",
    baronOptions().inlineparams);

  ////////////////// CUSTOM RESULT CODES ///////////////////
  AddSolveResults( {
                    { sol::LIMIT_NO_FEAS_ITER, "iteration limit, without a feasible soluton" },
                    { sol::LIMIT_NO_FEAS_NODES, "node limit, without a feasible soluton" },
                    { sol::LIMIT_NO_FEAS_TIME, "time limit, without a feasible solution" },
                    { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted, no feasible solution" },
                    { sol::FAILURE + 2, "too little memory" },
                    { sol::FAILURE + 1, "BARON syntax error (should not happen)" },
                    { sol::FAILURE, "licensing error" },
                    { sol::UNCERTAIN, "numerical difficulties but possibly optimal" },
                    { sol::UNKNOWN, "Unknown solution status returned from Baron"}
                   } );     
}

void BaronmpBackend::AddPrimalDualStart(Solution sol0_unpres) {
  if (IsMIP())
    return;
  auto mv = GetValuePresolver().PresolveSolution(
        { sol0_unpres.primal } );
  auto x0 = mv.GetVarValues()();
  if(x0.size() > 0){
    writeBaron("\nSTARTING_POINT {\n");
    for(int i=0; i<x0.size(); i++)
    writeBaron("{}: {};\n", lp()->varNames[i], x0[i]);
    writeBaron("}\n");
  }
}

void BaronmpBackend::AddMIPStart(
    ArrayRef<double> x0_unpres, ArrayRef<int> s0_unpres) {
    if (!IsMIP())
      return;
    auto mv = GetValuePresolver().PresolveSolution( { x0_unpres } );
    auto ms = GetValuePresolver().PresolveGenericInt( { s0_unpres } );
    auto x0 = mv.GetVarValues()();
    auto s0 = ms.GetVarValues()();
    bool headerWritten=false;
    for (int i=0; i<(int)x0.size(); ++i) {
    
      if (s0[i]) {
        if(!headerWritten) {
          headerWritten=true;
          writeBaron("\nSTARTING_POINT {\n");
        }
         writeBaron("{}: {};\n", lp()->varNames[i], x0[i]);
      }
    }
    if(headerWritten) writeBaron("}\n");
}

void BaronmpBackend::ComputeIIS() {
  // TODO
  //BARONMP_CCALL(BARONMP_ComputeIIS(lp()));
  SetStatus(GetSolveResult());   // could be new information
}

IIS BaronmpBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> BaronmpBackend::VarsIIS() {
  return ArrayRef<int>();
}
pre::ValueMapInt BaronmpBackend::ConsIIS() {

  return { {{ 0, std::vector<int>()}} };
}

double BaronmpBackend::MIPGap() {
  // implementation taken from ASL driver for consistency
  return MIPGapAbs() / (1e-10 + std::abs(ObjectiveValue()));
}
double BaronmpBackend::BestDualBound() {
  if(lp()->obj_sense == mp::obj::MAX)
    return timFileData_.printub;
  else
    return timFileData_.printlb;
}

double BaronmpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}

} // namespace mp


// AMPLs
void* AMPLSOpenBaronmp(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::BaronmpBackend()},
        cb);
}

void AMPLSCloseBaronmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetBaronmpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::BaronmpBackend*>(AMPLSGetBackend(slv))->lp();
}

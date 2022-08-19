#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "xpressmpbackend.h"

extern "C" {
  #include "xpressmp-ampls-c-api.h"    // Xpressmp AMPLS C API
  #include "xprs_mse_defaulthandler.h"
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptXpressmp(void* prob) {
  //return XPRESSMP_Interrupt((xpressmp_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateXpressmpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::XpressmpBackend()};
}


namespace mp {
  int XpressmpBackend::outlev_ = 0;
/// Create Xpressmp Model Manager
/// @param gc: the Xpressmp common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return XpressmpModelMgr
std::unique_ptr<BasicModelManager>
CreateXpressmpModelMgr(XpressmpCommon&, Env&, pre::BasicValuePresolver*&);


XpressmpBackend::XpressmpBackend() : msp_(NULL) {
  pre::BasicValuePresolver* pPre;
  auto data = CreateXpressmpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

XpressmpBackend::~XpressmpBackend() {
  CloseSolver();
}

void XpressmpBackend::InitOptionParsing() {
  OpenSolver();
}
void XpressmpBackend::OpenSolver() {
  int status = 0;
  // Try the registered function first; if not available 
  // call the solver's API function directly
  const auto& create_fn = GetCallbacks().init;
  if (create_fn)
    create_fn();
  else
    XPRESSMP_CCALL(XPRSinit(NULL));
  XPRESSMP_CCALL(XPRScreateprob(lp_ref()));
  copy_common_info_to_other();
  if (status)
    throw std::runtime_error("Error while creating Xpress environment");
}

void XpressmpBackend::CloseSolver() {
  if(lp())
    XPRESSMP_CCALL(XPRSdestroyprob(lp()));
}

const char* XpressmpBackend::GetBackendName()
  { return "XpressmpBackend"; }

std::string XpressmpBackend::GetSolverVersion() {
  char v[128];
  XPRSgetversion(v);
  return std::string(v);
}


bool XpressmpBackend::IsMIP() const {
  auto types = IsVarInt();
  for (auto a : types)
    if (a)
      return true;
  return false;
}

bool XpressmpBackend::IsQCP() const {
  return numQuadCons() > 0;
}

Solution XpressmpBackend::GetSolution() {
  auto mv = GetValuePresolver().PostsolveSolution(
    { PrimalSolution(), DualSolution(), GetObjectiveValues() });
  return { mv.GetVarValues()(), mv.GetConValues()(), mv.GetObjValues()() };
}

ArrayRef<double> XpressmpBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  if (IsMIP()) 
    error = XPRSgetmipsol(lp(), x.data(), NULL);
  else
    error = XPRSgetlpsol(lp(), x.data(), NULL, NULL, NULL);
  if (error)
    x.clear();
  return x;
}

pre::ValueMapDbl XpressmpBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> XpressmpBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  if (!IsMIP()) {
    int error = XPRSgetlpsol(lp(), NULL, NULL, pi.data(), NULL);
    if (error)
      pi.clear();
  }
  return pi;
}

double XpressmpBackend::ObjectiveValue() const {
  if (IsMIP())
    return getDblAttr(XPRS_MIPOBJVAL);
  else
    return getDblAttr(XPRS_LPOBJVAL);
}

double XpressmpBackend::NodeCount() const {
  return getIntAttr(XPRS_NODES);
}

double XpressmpBackend::SimplexIterations() const {
  return getIntAttr(XPRS_SIMPLEXITER);
}

int XpressmpBackend::BarrierIterations() const {
  return getIntAttr(XPRS_BARITER);
}

void XpressmpBackend::ExportModel(const std::string &file) {
  const char* s;
  char* wpflags = NULL;
  if (s = strrchr(file.c_str(), '.')) {
    if (!strcmp(s, ".mps"))
      wpflags = "";
    else if (!strcmp(s, ".lp"))
      wpflags = "l";
  }
  if (wpflags)
    XPRESSMP_CCALL(XPRSwriteprob(lp(), file.c_str(), wpflags));
  else 
    throw std::runtime_error(fmt::format("Expected \"writeprob=...\" to specify a filename ending in \".lp\"\n"
      "or \".mps\"; got \"{}\".\n", file));
}


void XpressmpBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptXpressmp, lp());
}

void XpressmpBackend::Solve() {
  int nsols = 10;
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  if (IsMIP()) {
    if (need_multiple_solutions() && storedOptions_.nbest_ > 1)
        XPRESSMP_CCALL(XPRS_mse_opt(mse_, lp(), 
          msp_, XPRS_mse_defaulthandler, 0, &storedOptions_.nbest_));
    else
      XPRESSMP_CCALL(XPRSmipoptimize(lp(), NULL));
  }
  else
    XPRESSMP_CCALL(XPRSlpoptimize(lp(), NULL));
  WindupXPRESSMPSolve();
}

void XpressmpBackend::WindupXPRESSMPSolve() { }

void XpressmpBackend::ReportResults() {
  ReportXPRESSMPResults();
  BaseBackend::ReportResults();
}

void XpressmpBackend::ReportXPRESSMPResults() {
  SetStatus( ConvertXPRESSMPStatus() );
  AddXPRESSMPMessages();
  if (need_multiple_solutions())
    ReportXPRESSMPPool();
}
std::vector<double> XpressmpBackend::getPoolSolution(int id)
{
  std::vector<double> vars(NumVars());
  int j, nx;
  XPRS_msp_getsol(msp_, id, &j, vars.data(), 0, NumVars() - 1, &nx);
  return vars;
}
double XpressmpBackend::getPoolObjective(int id)
{
  double obj;
  int j;
  XPRS_msp_getdblattribprobsol(msp_, lp(), id, &j, XPRS_MSP_SOLPRB_OBJ, &obj);
  return obj;
}
void XpressmpBackend::ReportXPRESSMPPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nPool;  
  XPRS_msp_getintattrib(msp_, XPRS_MSP_SOLUTIONS, &nPool);
  auto sid = std::vector<int>(nPool);
  int nret, nsols;
  XPRS_msp_getsollist(msp_, 0, 0, 0, 1, nPool, sid.data(),
    &nret, &nsols);
  nPool = nret;
  int id;
  for (int i = 0; i < nPool; i++) {
    id = sid[i++];
    auto mv = GetValuePresolver().PostsolveSolution({
          getPoolSolution(id) });
    ReportIntermediateSolution(
      { mv.GetVarValues()(), mv.GetConValues()(),
        { getPoolObjective(id) } });   // not when multiobj
  }
}


void XpressmpBackend::AddXPRESSMPMessages() {
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

std::pair<int, std::string> XpressmpBackend::ConvertXPRESSMPStatus() {
  namespace sol = mp::sol;
  if (IsMIP())
  {
    auto status = getIntAttr(XPRS_MIPSTATUS);
    switch (status) {
    case XPRS_MIP_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case XPRS_MIP_NO_SOL_FOUND:
      return { sol::INFEASIBLE, "infeasible problem" };
    case XPRS_MIP_UNBOUNDED:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case XPRS_MIP_LP_NOT_OPTIMAL:
      return { sol::INTERRUPTED, "interrupted" };
    }
  }
  else {
    auto status = getIntAttr(XPRS_LPSTATUS);
    switch (status) {
    case XPRS_LP_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case XPRS_LP_INFEAS:
      return { sol::INFEASIBLE, "infeasible problem" };
    case XPRS_LP_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case XPRS_LP_UNFINISHED:
      return { sol::INTERRUPTED, "interrupted" };
    case XPRS_LP_UNSTARTED:
      return { sol::UNKNOWN, "unstarted" };
    default:
      return { sol::UNKNOWN, "unfinished" };
    }
  }
  return { sol::UNKNOWN, "not solved" };
}

void XpressmpBackend::CreateSolutionPoolEnvironment() {
  // Solution pool must be created in all cases we specify a stub
  XPRESSMP_CCALL(XPRS_msp_create(&msp_));
  XPRESSMP_CCALL(XPRS_msp_probattach(msp_, lp()));

  // Create solutions enumerator only if explicitly needed
  if (storedOptions_.nbest_ > 1) {
    XPRESSMP_CCALL(XPRS_mse_create(&mse_));
    if(outlev_>0)
      XPRESSMP_CCALL(XPRS_mse_addcbmsghandler(mse_, xp_mse_display, NULL, 0));
  }

  SetSolverOption(XPRS_HEURSTRATEGY, 0);
  SetSolverOption(XPRS_MIPDUALREDUCTIONS, 2);

  if (storedOptions_.pooldualred_ != 2 || storedOptions_.pooldupcol_ != 2)
  {
    // Adapting presolve operations to solution pool specific controls
    int i, j = 0, k = 0;
    GetSolverOption(XPRS_PRESOLVEOPS, i);
    if (storedOptions_.pooldualred_ != 2) {
      j = 1 << 3;
      if (storedOptions_.pooldualred_ == 1)
        k = j;
    }
    if (storedOptions_.pooldupcol_ != 2) {
      j |= 1 << 5;
      if (storedOptions_.pooldupcol_ == 1)
        k |= 1 << 5;
    }
    i = (i & ~j) | k;
    SetSolverOption(XPRS_PRESOLVEOPS, i);
  }
}
void XpressmpBackend::FinishOptionParsing() {
  set_verbose_mode(outlev_ >0);
  if (outlev_ > 0)
    XPRSaddcbmessage(lp(), xpdisplay, NULL, 0);
  
  if (need_multiple_solutions())
    CreateSolutionPoolEnvironment();
}


////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo pool_values_[] = {
  { "0", "Yes (default, can be expensive)", 0},
  { "1", "No", 1},
  { "2", "Honor presolveops bit 3 (2^3 = 8)", 2}
};
static const mp::OptionValueInfo presolve_values_[] = {
  { "0", "No", 0},
  { "1", "Yes, removing redundant bounds (default)", 1},
  { "2", "Yes, retaining redundant bounds", 2}
};

static const mp::OptionValueInfo presolveops_values_[] = {
  { "1 = 2^0", "Remove singleton columns", 1},
  { "2 = 2^1", "Remove singleton constraints (rows)", 2},
  { "4 = 2^2", "Forcing row removal", 4},
  { "8 = 2^3", "Dual reductions", 8},
  { "16 = 2^4", "Redundant row removal", 16},
  { "32 = 2^5", "Duplicate column removal", 32},
  { "64 = 2^6", "Duplicate row removal", 64},
  { "128 = 2^7", "Strong dual reductions", 4},
  { "256 = 2^8", "Variable eliminations", 4},
  { "512 = 2^9", "No IP reductions", 4},
  { "1024 = 2^10", "No semi-continuous variable detection", 10},
  { "2048 = 2^11", "No advanced IP reductions", 11},
  { "4096 = 2^12", "No eliminations on integers", 12},
  { "16384 = 2^14", "Linearly dependant row removal", 14},
  { "32768 = 2^15", "Linearly dependant row removal", 15},
};

static const mp::OptionValueInfo pooldups_values_[] = {
  { "0", "Retain all duplicates", 0},
  { "1", "Discard exact matches", 1},
  { "2", "discard exact matches of continuous variables "
        "and matches of rounded values of discrete variables", 2},
  { "3", "Discard matches of rounded values of discrete "
  " variables (default)", 3}
};


void XpressmpBackend::InitCustomOptions() {

  set_option_header(
      "XPRESSMP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``xpressmp_options``. For example::\n"
      "\n"
      "  ampl: option xpressmp_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  AddStoredOption("tech:outlev outlev",
    "0*/1: Whether to write xpress log lines (chatter) to stdout and to file.",
    outlev_);

  // Solution pool params
  AddStoredOption("sol:pooldualred pooldualred",
    "Whether to suppress removal of dominated solutions(via "
      "\"dual reductions\") when poolstub is specified:\n"
    "\n.. value-table::\n",
    storedOptions_.pooldualred_, pool_values_);

  AddStoredOption("sol:pooldupcol pooldupcol",
    "Whether to suppress duplicate variable removal when "
		"poolstub is specified:\n"
    "\n.. value-table::\n",
    storedOptions_.pooldupcol_, pool_values_);

  AddSolverOption("sol:pooldups pooldups",
    "How poolstub should handle duplicate solutions:\n"
    "\n.. value-table::\nRounding of discrete variables is affected by"
    "poolmiptol and poolfeastol",
    XPRS_MSP_DUPLICATESOLUTIONSPOLICY, pooldups_values_, 3);

  AddSolverOption("sol:poolfeastol poolfeastol",
    "Zero tolerance for discrete variables in the solution "
		"pool (default 1e-6)",
    XPRS_MSP_SOL_FEASTOL, 0, 1);

  AddSolverOption("sol:poolmiptol poolmiptol",
    "Error (nonintegrality) allowed in discrete variables "
		"in the solution pool (default 5e-6)",
    XPRS_MSP_SOL_MIPTOL, 0, 1);

  AddStoredOption("sol:poolnbest poolnbest",
    "Whether the solution pool (see poolstub) should contain "
    "inferior solutions.  When poolnbest = n > 1, the "
    "solution pool is allowed to keep the n best solutions.",
    storedOptions_.nbest_);

  AddSolverOption("pre:solve presolve",
    "Whether to use Xpress' presolve:\n"
    "\n.. value-table::\n",
    XPRS_PRESOLVE, presolve_values_, 1);

  AddSolverOption("pre:ops presolveops", 
    "Reductions to use in XPRESS's presolve, sum of:\n"
    "\n.. value-table::\n(default 511 = bits 0-8 set)",
    XPRS_PRESOLVEOPS, presolveops_values_, 511);


}


double XpressmpBackend::MIPGap() {
  // implementation taken from ASL driver for consistency
  return MIPGapAbs() / (1e-10 + std::abs(ObjectiveValue()));
}
double XpressmpBackend::BestDualBound() {
  return getDblAttr(XPRS_BESTBOUND);
}

double XpressmpBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> XpressmpBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  XPRESSMP_CCALL(XPRSgetbasis(lp(), vars.data(), NULL));
  for (auto& s : vars) {
    switch (s) {
    case 1:
      s = (int)BasicStatus::bas;
      break;
    case 0:
      s = (int)BasicStatus::low;
      break;
    case 2:
      s = (int)BasicStatus::upp;
      break;
    case 3:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Xpressmp VBasis value: {}", s));
    }
  }
  return vars;
}

ArrayRef<int> XpressmpBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  XPRESSMP_CCALL(XPRSgetbasis(lp(), NULL, cons.data()));
  for (auto& s : cons) {
    switch (s) {
    case 1:
      s = (int)BasicStatus::bas;
      break;
    case 0:
      s = (int)BasicStatus::low;
      break;
    case 2:
      s = (int)BasicStatus::upp;
      break;
    case 3:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Xpressmp VBasis value: {}", s));
    }
  }
  return cons;
}

void XpressmpBackend::VarStatii(ArrayRef<int> vst) {
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  std::vector<double> lb, ub;
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = 1;
      break;
    case BasicStatus::low:
      s = 0;
      break;
    case BasicStatus::equ:
      s = 1;
      break;
    case BasicStatus::upp:
      s = 2;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = 3;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      if (lb.size() == 0)
      {
        lb.resize(vst.size());
        ub.resize(vst.size());
        XPRSgetlb(lp(), lb.data(), 0, vst.size());
        XPRSgetub(lp(), ub.data(), 0, vst.size());
      }
      if (lb[j] >= -1e-6)
        s = 0;
      else if (ub[j] <= 1e-6)
        s = 2;
      else
        s = 3;  
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }
  XPRESSMP_CCALL(XPRSloadbasis(lp(), NULL, stt.data()));
}

void XpressmpBackend::ConStatii(ArrayRef<int> cst) {
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = 1;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = 3;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  XPRESSMP_CCALL(XPRSloadbasis(lp(), stt.data(), NULL));
}

SolutionBasis XpressmpBackend::GetBasis() {
  // TODO: The following crashes 
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { varstt,
        {{{ CG_Linear, constt }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
    assert(constt.size());
  }
  
  
  //std::vector<int> varstt;
  //std::vector<int> constt;
  return { varstt,constt };

}

void XpressmpBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void XpressmpBackend::ComputeIIS() {
  int status;
  XPRESSMP_CCALL(XPRSiisfirst(lp(), 1, &status));
}

IIS XpressmpBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}
static IISStatus IIS_VarToAMPL(char c) {
  switch (c)
  {
  case 'L':
    return IISStatus::low;
  case 'U':
    return IISStatus::upp;
  case 'F':
    return IISStatus::fix;
  }
  return IISStatus::mem;
}
static IISStatus IIS_ConsToAMPL(char c) {
  switch (c)
  {
  case 'G':
    return IISStatus::low;
  case 'L':
    return IISStatus::upp;
  case '1':
  case '2':
  case 'E':
    return IISStatus::pmem;
  }
  return IISStatus::mem;
}

ArrayRef<int> XpressmpBackend::VarsIIS() {
  int nconsiis, nvarsiis;
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, 0, 0, 0, 0, 0, 0, 0, 0));
  std::vector<int> vars(nvarsiis);
  std::vector<char> bounds(nvarsiis), isolvars(nvarsiis);
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, 0,
    vars.data(), 0, bounds.data(), 0, 0, 0, isolvars.data()));
  std::vector<int> iis(NumVars(), 0);
  for (int i = 0; i < nvarsiis; i++) 
    iis[i] = (int)IIS_VarToAMPL(bounds[i]);
  return iis;
}
pre::ValueMapInt XpressmpBackend::ConsIIS() {
  int nconsiis, nvarsiis;
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, 0, 0, 0, 0, 0, 0, 0, 0));
  std::vector<int> cons(nconsiis);
  std::vector<char> contype(nconsiis), isolrows(nconsiis);
  XPRESSMP_CCALL(XPRSgetiisdata(lp(), 1, &nconsiis, &nvarsiis, cons.data(),
    0, contype.data(), 0, 0, 0, isolrows.data(), 0));
  std::vector<int> iis_lincon(NumLinCons(), 0);
  for (int i = 0; i < nvarsiis; i++)
    iis_lincon[i] = (int)IIS_VarToAMPL(contype[i]);
  return { {{ CG_Linear, iis_lincon }} }; // TODO other constraint types
}

void XpressmpBackend::AddMIPStart(ArrayRef<double> x0) {
  int status;
  XPRSloadmipsol(lp(), x0.data(), &status);
}

void XpressmpBackend::xpdisplay(XPRSprob prob, void* data, const char* ch, int n, int msglvl)
{
  /*
   msglvl gives the message level as follows:
   * 1 dialogue
   * 2 information
   * 3 warnings
   * 4 errors
   * a negative value indicates the XPRESS is about to finish and
   * buffers should be flushed.
   */
  if (msglvl < 0)
    fflush(NULL);
  else if (msglvl >= outlev_ && (msglvl != 4 || strncmp(ch, "?899 ", 5)))
    fmt::print("{}\n", ch);
}

int XpressmpBackend::xp_mse_display(XPRSobject o, void* context, void* thread, 
  const char* ch, int msglvl, int msgnumber)
{
  if (msglvl < 0)
    fflush(NULL);
  else if (msglvl >= outlev_ && (msglvl != 4))
    fmt::print("{}\n", ch);
  return 0;
}
} // namespace mp

// AMPLs
AMPLS_MP_Solver* AMPLSOpenXpressmp(
  const char* slv_opt) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::XpressmpBackend()},
    slv_opt);
}

void AMPLSCloseXpressmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

XPRSprob GetXpressmpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::XpressmpBackend*>(AMPLSGetBackend(slv))->lp();
}

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
  int XpressmpBackend::outlev_ = 3;
/// Create Xpressmp Model Manager
/// @param gc: the Xpressmp common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return XpressmpModelMgr
std::unique_ptr<BasicModelManager>
CreateXpressmpModelMgr(XpressmpCommon&, Env&, pre::BasicValuePresolver*&);


XpressmpBackend::XpressmpBackend() : msp_(NULL), mse_(NULL) {
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

bool XpressmpBackend::IsQCP() const {
  return numQuadCons() > 0;
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
  char const* wpflags = NULL;
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

void XpressmpBackend::DoXPRESSTune() {
  SetSolverOption(XPRS_TUNEROUTPUTPATH, tunebase().data());
  if (tunename().size())
    SetSolverOption(XPRS_TUNERSESSIONNAME, tunename().data());
  XPRESSMP_CCALL(XPRStune(lp(), ""));
}

void XpressmpBackend::Solve() {
  int nsols = 10;
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }

  if (tunebase().size())
    DoXPRESSTune();

  if (IsMIP()) {
    if (need_multiple_solutions() || storedOptions_.nbest_ > 1) {
      if (storedOptions_.nbest_ == 0) {
        storedOptions_.nbest_ = 20;
      }
      XPRESSMP_CCALL(XPRS_mse_opt(mse_, lp(),
        msp_, XPRS_mse_defaulthandler, 0,
        &storedOptions_.nbest_));
    }
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
      id = sid[i];
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
      case XPRS_MIP_INFEAS:
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
    if (storedOptions_.nbest_ >= 0) {
      XPRESSMP_CCALL(XPRS_mse_create(&mse_));
      if (outlev_ > 0)
        XPRESSMP_CCALL(XPRS_mse_addcbmsghandler(mse_, xp_mse_display, NULL, 0));
    }

    SetSolverOption(XPRS_HEURSTRATEGY, 0);
    SetSolverOption(XPRS_MIPDUALREDUCTIONS, 2);

    if (storedOptions_.pooldualred_ != 2 || storedOptions_.pooldupcol_ != 2)
    {
      // Adapt presolve operations to solution pool specific controls
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
    bool doLog = outlev_ > 0 && outlev_ < 5;
    set_verbose_mode(doLog);
    if (doLog)
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
    { "1 = 2^0", "Remove singleton columns", XPRS_PRESOLVEOPS_SINGLETONCOLUMNREMOVAL},
    { "2 = 2^1", "Remove singleton constraints (rows)", XPRS_PRESOLVEOPS_SINGLETONROWREMOVAL},
    { "4 = 2^2", "Forcing row removal", XPRS_PRESOLVEOPS_FORCINGROWREMOVAL},
    { "8 = 2^3", "Dual reductions", XPRS_PRESOLVEOPS_DUALREDUCTIONS},
    { "16 = 2^4", "Redundant row removal", XPRS_PRESOLVEOPS_REDUNDANTROWREMOVAL},
    { "32 = 2^5", "Duplicate column removal", XPRS_PRESOLVEOPS_DUPLICATECOLUMNREMOVAL},
    { "64 = 2^6", "Duplicate row removal", XPRS_PRESOLVEOPS_DUPLICATEROWREMOVAL},
    { "128 = 2^7", "Strong dual reductions", XPRS_PRESOLVEOPS_STRONGDUALREDUCTIONS},
    { "256 = 2^8", "Variable eliminations", XPRS_PRESOLVEOPS_VARIABLEELIMINATIONS},
    { "512 = 2^9", "No IP reductions", XPRS_PRESOLVEOPS_NOIPREDUCTIONS},
    { "1024 = 2^10", "No semi-continuous variable detection", XPRS_PRESOLVEOPS_NOGLOBALDOMAINCHANGE},
    { "2048 = 2^11", "No advanced IP reductions", XPRS_PRESOLVEOPS_NOADVANCEDIPREDUCTIONS},
    { "4096 = 2^12", "No eliminations on integers", 12},
    { "16384 = 2^14", "Linearly dependant row removal", XPRS_PRESOLVEOPS_LINEARLYDEPENDANTROWREMOVAL},
    { "32768 = 2^15", "No integer variable and SOS detection", XPRS_PRESOLVEOPS_NOINTEGERVARIABLEANDSOSDETECTION},
    { "536870912 = 2^29", "No dual reduction on globals", XPRS_PRESOLVEOPS_NODUALREDONGLOBALS},
  };

  static const mp::OptionValueInfo pooldups_values_[] = {
    { "0", "Retain all duplicates", 0},
    { "1", "Discard exact matches", 1},
    { "2", "Discard exact matches of continuous variables "
          "and matches of rounded values of discrete variables", 2},
    { "3", "Discard matches of rounded values of discrete "
    " variables (default)", 3}
  };

  static const mp::OptionValueInfo values_method[] = {
    { "1", "Automatic choice (default)", 1},
    { "2", "Dual simplex", 2},
    { "3", "Primal simplex", 3},
    { "4", "Netwon Barrier", 4}
  };

  static const mp::OptionValueInfo values_baralg[] = {
    { "-1", "Automatic choice (default)", -1},
    { "1", "Infeasible-start barrier algorithm", 1},
    { "2", "Homogeneous self-dual barrier algorithm", 2},
    { "3", "Start with 2 and maybe switch to 1 while solving", 3}
  };

  static const mp::OptionValueInfo values_barstart[] = {
    { "-1", "Use incoming solution for warm start", -1},
    { "0",  "Automatic choice (default)", 0},
    { "1", "Heuristics based on magnitudes of matrix entries", 1},
    { "2", "Use pseudoinverse of constraint matrix", 2},
    { "3", "Unit starting point for homogeneous self - dual "
          "barrier algorithm.", 3}
  };

  static const mp::OptionValueInfo values_barcholeskyalg[] = {
    {"-1", "Automatic choice (default)", -1},
    {"1", "manual matrix blocking", 1},
    {"2", "manual blocking: single-pass (multi-pass if not set)", 2},
    {"4", "nonseparable QP relaxation", 4},
    {"8", "manual corrector weight (honor \"16\" bit)", 8},
    {"16", "manual corrector weight \"on\"", 16},
    {"32", "manual refinement", 32},
    {"64", "use preconditioned conjugate gradients", 64},
    {"128", "refine with QMR (quasi-minimal residual)", 128},
    {"256", "perform refinement on the augmented system", 256},
    {"512", "force highest accuracy in refinement", 512}

  };
  static const mp::OptionValueInfo values_barcrash[] = {
    { "0", "No crash", 0},
    { "1-6", "Available strategies", 2}
  };

  static const mp::OptionValueInfo values_barkernel[] = {
    { ">= +1.0", "More emphasis on centrality (default 1.0)", 0},
    { "<= -1.0", "Each iteration, adaptively select a value "
    "from [+1, -barkernel]", 2}
  };


  static const mp::OptionValueInfo values_barobjperturb[] = {
    { "n > 0", "automatic decison, scale n", 1},
    { "n = 0", "turn off perturbation", 0},
    { "n < 0", "force perturbation by abs(n)", -1}
  };

  static const mp::OptionValueInfo values_barobjscale[] = {
    {"-1", "Automatic choice (default)", -1},
    {"0", "Scale by the geometric mean of the objective coefficients", 0},
    {"> 0", "Scale so the argest objective coefficient in absolute "
            "value is <= barobjscale.", 1}
  };

  static const mp::OptionValueInfo values_barorder[] = {
    { "0", "automatic choice(default)", 0},
    { "1", "minimum degree", 1},
    { "2", "minimum local fill", 2},
    { "3", "nested dissection",3}
  };

  static const mp::OptionValueInfo values_baroutput[] = {
    { "0", "no output", 0},
    { "1", "each iteration (default)", 1} 
  };

  static const mp::OptionValueInfo barpresolve_desc[] = {
    { "0", "use standard presolve (default)", 0},
    { "1", "use more effort", 1},
    {"2", "do full matrix eliminations for size reduction", 2}
  };

  static const mp::OptionValueInfo values_barrefiter[] = {
    { "0", "default", 0},
    { "n > 0", "perform n refinement iterations", 1}
  };

  static const mp::OptionValueInfo values_barregularize[] = {
    {"1", "use \"standard\" regularization", 1},
    {"2", "use \"reduced\" regularization: less perturbation than \"standard\" regularization", 2},
    {"4", "keep dependent rows in the KKT system",4},
    {"8", "keep degenerate rows in the KKT system", 8}
 };
  static const mp::OptionValueInfo values_bigmmethod[] = {
    { "0", "phase I / II", 0},
    { "1", "bigM method (default)", 1}
  };

  static const mp::OptionValueInfo values_branchchoice[] = {
    {"0", "explore branch with min.estimate first(default)", 0},
    {"1", "explore branch with max.estimate first", 1},
    {"2", "if an incumbent solution exists, first explore "
          "the branch satisfied by the incumbent; "
          "otherwise use choice 0 (min.est.first).", 2},
    {"3", "(default) explore the first branch that moves the "
          "branching variable away from its value at the "
          "root node; if the branching entity is not a "
          "simple variable, assume branchchoice=0.", 3}
  };

  static const mp::OptionValueInfo values_branchdisj[] = {
    {"-1", "automatic choice (default)", 0},
    {"0", "disabled", 0},
    {"1", "cautious strategy : create branches only for "
        "general integers with a wide range", 1},
    {"2", "moderate strategy", 2},
    {"3", "aggressive strategy : create disjunctive branches "
    "for both binaryand integer variables", 3}
  };

  static const mp::OptionValueInfo values_clamping[] = {
    {"-1", "Determined automatically", -1},
    {"0", "Adjust primal solution to always be within primal bounds (default)", 0},
    {"1", "Adjust primal slack values to always be within constraint bounds", 1},
    {"2", "Adjust dual solution to always be within the dual bounds implied by the slacks", 2},
    {"3", "Adjust reduced costs to always be within dual bounds implied by the primal solution", 3}
  };

  static const mp::OptionValueInfo values_cpuplatform[] = {
  {"-2", "Highest supported [Generic, SSE2, AVX or AVX2]", -2},
  {"-1", "Highest supported solve path consistent code [Generic, SSE2 or AVX] (default)", -1},
  {"0", "Use generic code compatible with all CPUs", 0},
  {"1", "Use SSE2", 1},
  {"2", "Use AVX", 2},
  {"3", "Use AVX2", 3}
  };

  static const mp::OptionValueInfo values_cputime[] = {
    {"-1", "disable the timer", -1},
    {"0", "use elapsed time (default)", 0},
    {"1", "use process time", 1}
  };

  static const mp::OptionValueInfo values_crash[] = {
    {"0", "none", 0},
    {"1", "one-pass search for singletons", 1},
    {"2", "multi-pass search for singletons", 2},
    {"3", "multi-pass search including slacks", 3},
    {"4", "at most 10 passes, only considering slacks at the end", 4},
    {"n>10", " like 4, but at most n-10 passes", 5},
    {"0 (dual)", "perform standard crash.", 0},
    {"1 (dual)", "perform additional numerical checks during crash", 1},
    {"2 (dual)", "extend the set of column candidates for crash", 2},
    {"3 (dual)", "extend the set of row candidates for crash", 3},
    {"4 (dual)", "force crash", 4}
  };

  static const mp::OptionValueInfo values_barcrossover[] = {
  {"-1", "automatic choice (default)", -1},
  { "0", "none: return an interior solution", 0},
  { "1", "primal crossover first", 1},
  { "2", "dual crossover first", 2}
  };

  static const mp::OptionValueInfo values_barcrossoverops[] = {
 { "1", "return the barrier solution (rather than the last intermediate solution) when crossover stop early", 1},
 { "2", "skip the second crossover stage", 2},
 { "4", "skip pivots that are \"less numerically reliable\"", 4},
 { "8", "do a slower but more numerically stable crossover", 8},
  };

  static const mp::OptionValueInfo values_cutselect[] = {
    {"32", "clique cuts", 32},
    {"64", "mixed - integer founding(MIR) cuts", 64},
    {"128", "lifted cover cuts", 128},
    {"2048", "flow path cuts", 2048},
    {"4096", "implication cuts", 4096},
    {"8192", "automatic lift - and -project strategy", 8192},
    {"16384", "disable cutting from cut rows", 16384},
    {"32768", "lifted GUB cover cuts", 32768},
    {"65536", "zero - half cuts", 65536},
    {"131072", "indicator - constraint cuts", 131072},
    {"-1", "all available cuts(default)", -1}
  };
  static const mp::OptionValueInfo values_cutstrategy[] = {
    {"-1", "automatic (default)", -1},
    {"0", "no cuts", 0},
    {"1", "conservative strategy", 1},
    {"2", "moderate strategy", 2},
    {"3", "aggressive strategy", 3}
  };

  static const mp::OptionValueInfo values_deterministic[] = {
    {"0", "no", 0},
    {"1", "yes", 1},
    {"2", "yes, with opportunistic root LP solve", 1},
  };

  static const mp::OptionValueInfo values_lpdualgradient[] = {
   {"-1", "automatic (default)", -1},
   {"0", "devex", 0},
   {"1", "steepest edge", 1},
   {"2", "direct steepest edge", 2},
   {"3", "sparse devex", 3}
  };

  static const mp::OptionValueInfo values_lpdualize[] = {
 {"-1", "automatic (default)", -1},
 {"0", "solve the primal problem", 0},
 {"1", "solve the dual problem", 1}
  };
  static const mp::OptionValueInfo values_lpdualstrategy[] = {
    {"1", "switch to primal when dual infeasible", 1},
    {"2", "stop the solve instead of switching to primal", 2},
    {"4", "use aggressive cut-off in MIP search", 4},
    {"8", "use dual simplex to remove cost perturbations", 8},
    {"16", "aggressive dual pivoting", 16},
    {"32", "keep using dual simplex even when numerically unstable", 32}
  };
#define xstr(a) stringify(a)
#define stringify(a) #a

  static const mp::OptionValueInfo values_feasibilitypump[] = {
    {xstr(XPRS_FEASIBILITYPUMP_AUTOMATIC), "automatic (default)", XPRS_FEASIBILITYPUMP_AUTOMATIC},
    {xstr(XPRS_FEASIBILITYPUMP_NEVER), "turned off", },
    {xstr(XPRS_FEASIBILITYPUMP_ALWAYS), "always run", 2},
    {xstr(XPRS_FEASIBILITYPUMP_LASTRESORT), "run if other heuristics have failed to find an integer solution", 4}
  };
  static const mp::OptionValueInfo values_heurdivesoftrounding[] = {
    {"-1", "automatic (default)", -1},
    {"0", "do not use soft rounding", 0},
    {"1", "cautious use", 1},
    {"2", "aggressing use", 2}
  };

  static const mp::OptionValueInfo values_heurdivespeed[] = {
    {"-2", "automatic bias toward quality", -2},
    {"-1", "automatic bias toward speed (default)", -1},
    {"0", "emphasize quality", 0},
    {"1-3", "intermediate emphasis", 1},
    {"4", "emphasize speed", 4}
  };

  static const mp::OptionValueInfo values_heurdivestrategy[] = {
   {"-1", "automatic selection (default)", -1},
   {"0", "disable heuristics", 0},
   {"1-18", "available pre-set strategies for rounding infeasible global entities", 1},
  };

    static const mp::OptionValueInfo values_heuremphasis[] = {
            {"-1", "default strategy (default)", -1},
            {"0", "disable heuristics", 0},
            {"1", "focus on reducing the gap early", 1},
            {"2", "extremely aggressive heuristics", 2},
    };

    static const mp::OptionValueInfo values_heursearchfreq[] = {
            {"-1", "automatic (default)", -1},
            {"0", "disabled in the tree", 0},
            {"n>0", "number of nodes between each run", 1}
    };

    static const mp::OptionValueInfo values_heursearchrootcutfreq[] = {
            {"-1", "automatic (default)", -1},
            {"0", "disabled during cutting", 0},
            {"n>0", "number cutting rounds between each run", 1}
    };

    static const mp::OptionValueInfo values_heursearchrootcutselect[] = {
            {"1", "local search with a large neighborhood. Potentially slow but is good for finding solutions that differs significantly from the incumbent", 1},
            {"2", "local search with a small neighborhood centered around a node LP solution", 2},
            {"4", "local search with a small neighborhood centered around an integer solution. This heuristic will often provide smaller, incremental improvements to an incumbent solution", 4},
            {"8", "local search with a neighborhood set up through the combination of multiple integer solutions.", 8},
            {"32", "local search without an objective function", 32},
            {"64", "local search with an auxiliary objective function", 64}
    };

    static const mp::OptionValueInfo values_heurthreads[] = {
            {"-1", "determined from \"threads\" keyword", -1},
            {"0", "no separate threads (default)", 0},
            {"n>0", "use n threads", 1}
    };

    static const mp::OptionValueInfo values_historycosts[] = {
            {"-1", "automatic (default)", -1},
            {"0", "no update", 0},
            {"1","initialize using only regular branches from the root to the current node", 1},
            {"2", "same as 1, but initialize with strong branching results as well", 2},
            {"3", "initialize using any regular branching or strong branching information from all nodes solves before the current node", 3}
    };

    static const mp::OptionValueInfo values_keepbasis[] = {
            {"0", "ignore previous basis", 0},
            {"1", "use previous basis (default)", 1},
            {"2", "use previous basis only if the number of basic variables == number of constraints", 2}
    };

    static const mp::OptionValueInfo values_keepnrows[] = {
            {"-1", "delete N type rows from the matrix (default)", -1},
            {"0", "delete elements from N type rows leaving empty N type rows in the matrix", 0},
            {"1","keep N type rows", 1}
    };

    static const mp::OptionValueInfo values_localchoice[] = {
            {"1", "never backtrack from the first child unless it is dropped "
                  "(i.e., is infeasible or cut off) (default)", 1},
            {"2", "always solve both child nodes before deciding which child to continue with", 2},
            {"3", "automatically determined", 2},
    };

    static const mp::OptionValueInfo values_mipconcurrentnodes[] = {
           {"-1", "automatic (default)", -1},
           {"n > 0", "number of nodes to complete", 1}
    };

    static const mp::OptionValueInfo values_mipconcurrentsolves[] = {
       {"-1", "enabled, the number of concurrent solves "
				  "depends on mipthreads", -1},
       {"0", "disabled (default)", 0},
       {"1", "disabled (default)", 1},
      {"n > 1", "number of concurrent solves = n", 2},
    };

    static const mp::OptionValueInfo values_mipdualreductions[] = {
      {"0", "none", 0},
      {"1", "all (default)", 1},
      {"2", "restrict dual reductions to continuous variables", 2},
    };

  static const mp::OptionValueInfo values_mipkappafreq[] = {
  {"0", "never (default)", 0},
  {"1", "every node", 1},
  {"n > 1", "once per node at level n of the branch-and-bound tree", 2} };

  static const mp::OptionValueInfo values_qsimplexops[] = {
    { "1", "traditional primal first phase (default)", 1},
    { "2", "force Big M primal first phase", 2},
    { "4", "force traditional dual first phase", 4},
    { "8", "force BigM dual first phase", 8},
    { "16", "always use artificial bounds in dual", 16},
    { "32", "use original problem basis only when warmstarting the KKT", 32},
    { "64", "skip the primal bound flips for ranged primals (might cause "
    "more trouble than good if the bounds are very large)", 64},
    { "128", "also do the single pivot crash", 128},
    { "256", "do not apply aggressive perturbation in dual", 256}
};

  static const mp::OptionValueInfo values_mippresolve[] = {
  { "1", "reduced-cost fixing at each node", 1},
  { "2", "primal reductions will be performed at each node", 2},
  { "8", "allow changing continuous-variable bounds", 8},
  { "16", "allow dual reductions", 16},
  { "32", "allow global tightening of the problem", 32},
  { "64", "use objective function", 64},
  { "128", "allow restarting", 128},
  { "256", "allow use of symmetry", 256}
  };

  static const mp::OptionValueInfo values_miprampup[] = {
   {"-1", "automatic choice (default)", -1},
   {"0", "no: use as many tasks as possible", 0},
   {"1", "yes, until finished with initial dives", 1}
  };
  static const mp::OptionValueInfo values_miprestart[] = {
   {"-1", "automatic choice (default)", -1},
   {"0", "disable in-tree restarts", 0},
   {"1", "normal aggressiveness", 1},
   {"2", "higher aggressiveness", 2},
  };

  static const mp::OptionValueInfo values_miqcpalg[]{
    {"-1", "automatic (default)", -1},
    {"0", "barrier", 0},
    {"1", "outer approximations", 1}
  };
  static const mp::OptionValueInfo values_netstalllimit[]{
  {"-1", "automatic (default)", -1},
  {"0", "no limit", 0},
  {"n > 0", "limit to n network simplex iterations", 1}
  };

  static const mp::OptionValueInfo values_nodeselection[]{
    {"1", "local first: choose between descendant and sibling nodes if available; choose from all outstanding nodes otherwise", 1},
    {"2", "best first: choose from all outstanding nodes", 2},
    {"3", "local depth first: choose between descendant and sibling nodes if available; choose from the deepest nodes otherwise", 3},
    {"4", "best first, then local first: best first is used for the first BREADTHFIRST nodes, after which local first is used", 4},
    {"5", "pure depth first: choose from the deepest outstanding nodes", 5}
  };
  static const mp::OptionValueInfo values_numericalemphasis[]{
        {"-1", "automatic (default)", -1},
        {"0", "epmhasize speed", 0},
        {"1", "mild emphasis on numerical stability", 1},
        {"2", "medium emphasis on numerical stability", 2},
        {"3", "strong emphasis on numerical stability", 3},
  };

  static const mp::OptionValueInfo values_outlev[]{
  {"0","none", 0},
  {"1", "all", 1},
  {"2", "information", 2},
  {"3", "warnings & errors only (default)", 3},
  {"4", "errors", 4},
  {"5","none", 5}
  };
  static const mp::OptionValueInfo values_precoefelim[]{
  {"0","disabled", 0},
  {"1", "remove as many coefficients as possible", 1},
  {"2", "cautious eliminations", 2}
};
  static const mp::OptionValueInfo values_preconvertseparable[]{
 {"-1","automatic (default)", -1},
    {"0", "disable", 0},
 {"1", "enable reformulation to diagonal quadratic constraints.", 1},
 {"2", "1, plus reduction to second-order cones", 2},
 {"3", "2, plus the objective function is converted to a constraint and treated as a quadratic constraint", 3}
 };

  static const mp::OptionValueInfo values_predomcol[]{
{"-1","automatic (default)", -1},
{"0", "disable", 0},
{"1", "cautious", 1},
{"2", "aggressive: all candidate will be checked", 2}
  };

  static const mp::OptionValueInfo values_predomrow[] = {
 {"-1", "automatic choice (default)", -1},
 {"0", "disabled", 0},
 {"1", "cautious", 1},
 {"2", "medium", 2},
 {"3", "aggressive", 3}
  };

  static const mp::OptionValueInfo values_preduprow[]{
    {"-1","automatic (default)", -1},
    {"0", "disable", 0},
    {"1", "eliminate only rows that are identical in all variables", 1},
    {"2", "1 plus eliminate duplicate rows with simple penalty variable expressions", 2},
    {"3", "2 plus eliminate duplicate rows with more complex penalty variable expressions", 3}
  };
  
  static const mp::OptionValueInfo values_prepermute[]{
    {"1", "permute rows", 1},
    {"2", "permute columns", 2},
    {"4", "permute global information (for MIP)", 4}
  };

  static const mp::OptionValueInfo values_pricingalg[]{
   {"-1", "partial pricing", 1},
   {"0", "automatic choice (default)", 0},
   {"1", "devex pricing", 1},
   {"2", "steepest edge", 2},
    {"3", "steepest edge with initial weights", 3}
  };

  const mp::OptionValueInfo values_yesnoinverted_defaultyes[] = {
  {     "0", "yes (default)", 0 },
  {     "1", "no", 1}
  };
  const mp::OptionValueInfo values_yesnoinverted_defaultno[] = {
  {     "0", "yes", 0 },
  {     "1", "no (default)", 1}
  };

  const mp::OptionValueInfo values_pwlnonconvextransformation[] = {
    {"-1", "automatic (default)", -1},
    {     "0", "use a formulation based on SOS2-constraints", 0 },
    {     "1", "use a formulation based on binary variables", 1}
  };

  const mp::OptionValueInfo values_qcrootalg[] = {
  {"-1", "automatic (default)", -1},
  {     "0", "use barrier", 0 },
  {     "1", "use dual simplex on outer approximation", 1}
  };

  const mp::OptionValueInfo values_refineops[] = {
 { "1", "refine optimal LP solutions", 1},
    { "2", "refine MIP solutions", 2},
    { "8", "refine each node of the search tree", 8},
    { "16", "refine non-global solutions", 16},
    { "32", "apply the iterative refiner to refine the solution", 32},
    { "64", "use higher precision in the iterative refinement", 64},
    { "128", "iterative refiner will use the primal simplex algorithm", 128},
    { "256", "iterative refiner will use the dual simplex algorithm", 256},

    { "512", "refine MIP solutions such that rounding them keeps the problem "
    "feasible when reoptimized", 512},
    { "1024", "ttempt to refine MIP solutions such that rounding them keeps the "
    "problem feasible when reoptimized, but accept integers solutions even if "
    "refinement fails", 1024},
  };

  const mp::OptionValueInfo values_sbbest[] = {
{"-1", "automatic (default)", -1},
{"0", "disable strong branching", 0 },
{"n > 1", "perform strong branching on up to n entities at each node", 1}
  };
  const mp::OptionValueInfo values_sbestimate[] = {
{"-1", "automatic (default)", -1},
{"1-6", "different variants of local pseudo costs.", 1}
  };

  const mp::OptionValueInfo values_sbselect[] = {
    {"-2", "automatic - low effort (default)", -2},
    {"-1", "automatic - high effort", -1},
    {"n > 0", "include max(n, sbbest) candidates", 1},
  };
  const mp::OptionValueInfo values_scaling[] = {
    {"1", "row scaling", 1},
    {"2", "column scaling", 2},
    {"4", "row scaling again", 4},
    {"8", "maximum", 8},
    {"16", "Curtis-Red", 16},
    {"32", "0->geometric mean, 1->maximum element", 32},
    {"64", "no special handling for BigM rows", 64},
    {"128", "scale objective function for the simplex method", 128},
    {"256", "exclude the quadratic part of constraints when calculating scaling factors", 256},
    {"512", "scale before presolve", 512},
    {"1024", "do not scale constraints up", 1024},
    {"2048", "do not scale variables down", 2048},
    {"4096", "do not apply automatic global objective scaling", 4096},
    {"8192", "RHS scaling", 8192},
    {"16384", "disable aggressive quadratic scaling", 16384},
    {"32768", "explicit linear slack scaling", 32768}
  };

  const mp::OptionValueInfo values_siftpresolveops[] = {
    {"-1", "use the  \"presolveops\" setting specified for the original problem", -1},
    {">=0", "use the value (see \"presolveops\" for its semantic)", 0}
  };

  const mp::OptionValueInfo values_siftswitch[] = {
    {"-1", "dual simplex", -1},
    {"0", "barrier", 0},
    {">0", "use the barrier algorithm while the number of dual infeasibilities is larger than this value, otherwise use dual simplex", 1}
  };

  const mp::OptionValueInfo values_sleeponthreadwait[] = {
   {"-1", "automatically determined", -1},
   {"0", "no (busy-wait)", 0},
   {">0", "yes (sleep, might add overhead)", 1}
  };

  const mp::OptionValueInfo values_symmetry[] = {
    {"0", "no simmetry detection", 0},
    {"1", "conservative effort", 1},
    {"2", "intensive effort", 2}
  };
  const mp::OptionValueInfo values_symselect[] = {
    {"0", "search the whole matrix (otherwise the 0, 1 and -1 coefficients only)", 0},
    {"1", "search all entities(otherwise binaries only)", 1} };

  const mp::OptionValueInfo values_tunerhistory[] = {
    {"0", "Discard any previous result", 0},
    {"1", "Append new results but do not reuse them", 1},
    {"2", "Reuse and append new results", 2}
  };

  const mp::OptionValueInfo values_tunermethod[] = {
      {"- 1", "automatic choice(default)", -1},
      {"0", "default LP tuner", 0},
      {"1", "default MIP tuner", 1},
      {"2", "more elaborate MIP tuner", 2},
      {"3", "root - focused MIP tuner", 3},
      {"4", "tree - focused MIP tuner", 4},
      {"5", "simple MIP tuner", 5},
      {"6", "default SLP tuner", 6},
      {"7", "default MISLP tuner", 7},
      {"8", "MIP tuner using primal heuristics", 8}
  };
  const mp::OptionValueInfo values_tunertarget[] = {
    {"- 1", "automatic choice(default)",-1},
    {"0", "solution time, then integrality gap", 0},
      {"1", "solution time, then best bound",1},
       {"2", "solution time, then best integer solution",2},
       {"3", "the \"primal dual integral\", whatever that is",3},
       {"4", "just solution time (default for LPs)",4},
       {"5", "just objective value" ,5},
       {"6", "validation number (probably not relevant)",6},
       {"7", "gap only" , 7},
       {"8", "best bound only" , 8},
       {"9", "best integer solution only" , 9},
       {"10", "best primal integral - only for individual instances" , 10}
  }; 
  
  const mp::OptionValueInfo values_varselection[] = {
    {"- 1", "automatic choice(default)", -1},
    {"1", "minimum of the 'up' and 'down' pseudo - costs", 1},
    {"2", "'up' pseudo - cost + 'down' pseudo - cost", 2},
    {"3", "maximum of the 'up' and 'down' pseudo - costs plus twice their minimum", 3},
    {"4", "maximum of the 'up' and 'down' pseudo - costs", 4},
    {"5", "the 'down' pseudo - cost", 5},
    {"6", "the 'up' pseudo - cost", 6},
    {"7", "weighted combination of the 'up' and 'down' pseudo costs", 7},
    {"8", "product of 'up' and 'down' pseudo costs", 8}
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

  // ****************************
  // General
  // ****************************
  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

  AddStoredOption("tech:outlev outlev",
    "Whether to write xpress log lines (chatter) to stdout and to file:\n"
    "\n.. value-table::\n",
    outlev_, values_outlev);

  AddSolverOption("tech:cputime cputime",
    "How time should be measured when timings are reported in the log and when checking against time limits :\n"
    "\n.. value-table::\n",
    XPRS_CPUTIME, values_cputime, 0);

  AddSolverOption("tech:threads threads",
    "The default number of threads used during optimization.;"
    "default - 1 ==> automatic choice.",
    XPRS_THREADS, -1, INT_MAX);
  
  AddSolverOption("lim:time timelim timelimit",
    "Limit on solve time (in seconds; default: no limit). If n>0, stop MIP search only"
    "after a solution has been found, if n<0, stop after -n seconds nevertheless",
    XPRS_MAXTIME, -INT_MAX, INT_MAX);

  AddSolverOption("tech:sleeponthreadwait sleeponthreadwait",
    "Whether threads should sleep while awaiting work:\n"
    "\n.. value-table::\n",
    XPRS_SLEEPONTHREADWAIT, values_sleeponthreadwait,  -1);

    AddSolverOption("lim:lpiterlimit lpiterlimit",
        "The maximum number of iterations that will be performed "
        "by primal simplex or dual simplex before the optimization "
        "process terminates. For MIP problems, this is the maximum "
        "total number of iterations over all nodes.",
        XPRS_LPITERLIMIT, 2147483645, INT_MAX);

    AddSolverOption("lim:lprefineiterlimit lprefineiterlimit",
                    "This specifies the simplex iteration limit the solution "
                    "refiner can spend in attempting to increase the accuracy "
                    "of an LP solution; default=-1 (automatic).",
                    XPRS_LPREFINEITERLIMIT, -1, INT_MAX);

    AddSolverOption("lim:maxcuttime maxcuttime",
                    "The maximum amount of time allowed for generation of cutting planes and reoptimization;"
                    "default=0 (no time limit)",
                    XPRS_MAXCUTTIME, 0, INT_MAX);

    AddSolverOption("tech:globalfilemax globalfilemax",
    "Maximum megabytes for temporary files storing the global search "
    "tree: a new file is started if globalfilemax megabytes would be exceeded.",
    XPRS_MAXGLOBALFILESIZE, 0, INT_MAX);

  AddSolverOption("tech:globalfileloginterval globalfileloginterval",
    "Seconds between additions to the logfile about, additions "
		"to the \"global file\", a temporary file written during a "
		"global search. Default = 60.",
    XPRS_GLOBALFILELOGINTERVAL, 1, INT_MAX);

  AddSolverOption("lim:memlimit memlimit maxmemoryhard", 
                  "Hard limit (integer number of MB) on memory allocated, "
                  "causing early termination if exceeded; default = 0 (no limit)",
                  XPRS_MAXMEMORYHARD, 0, INT_MAX);

    AddSolverOption("lim:softmemlimit softmemlimit maxmemorysoft",
                    "Soft limit (integer number of MB) on memory allocated; "
                    "default = 0 (no limit)",
                    XPRS_MAXMEMORYSOFT, 0, INT_MAX);
  // ****************************
  // Solution pool params
  // ****************************
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

  AddSolverOption("sol:pooldups poold/ups",
    "How poolstub should handle duplicate solutions:\n"
    "\n.. value-table::\nRounding of discrete variables is affected by"
    "poolmiptol and poolfeastol",
    XPRS_MSP_DUPLICATESOLUTIONSPOLICY, pooldups_values_, 3);

  AddOptionSynonyms_Inline_Front("ams_stub", "sol:stub");

  AddSolverOption("sol:poolfeastol poolfeastol",
    "Zero tolerance for discrete variables in the solution "
    "pool (default 1e-6)",
    XPRS_MSP_SOL_FEASTOL, 0, 1);

  AddSolverOption("sol:poolmiptol poolmiptol",
    "Error (nonintegrality) allowed in discrete variables "
    "in the solution pool (default 5e-6)",
    XPRS_MSP_SOL_MIPTOL, 0, 1);

  AddStoredOption("sol:poolnbest poolnbest poollimit",
    "Whether the solution pool (see poolstub) should contain "
    "inferior solutions.  When poolnbest = n > 1, the "
    "solution pool is allowed to keep the n best solutions.",
    storedOptions_.nbest_);

    AddSolverOption("lim:maxmipsol maxmipsol",
                    "Limit on the number of MIP solutions to be found (default no limit).",
                    XPRS_MAXMIPSOL, 0, INT_MAX);

    AddSolverOption("lim:nodes nodelim nodelimit maxnode",
                    "Maximum MIP nodes to explore (default: 2147483647).",
                    XPRS_MAXNODE, 2147483647, INT_MAX);

    AddSolverOption("lim:maxstalltime maxstalltime",
                    "Maximum time in seconds that the MIP Optimizer will continue to search "
                    "for improving solution after finding a new incumbent, default=0 (no limit)",
                    XPRS_MAXSTALLTIME, 0, INT_MAX);
  // ****************************
  // Generic algorithm controls
  // ****************************
  AddSolverOption("alg:method method lpmethod defaultalg",
    "Which algorithm to use for non-MIP problems or for the root node of MIP problems:\n"
    "\n.. value-table::\n", XPRS_DEFAULTALG, values_method, -1);


  AddSolverOption("alg:clamping clamping",
    "Control adjustements of the returned solution values "
     "such that they are always within bounds:\n"
    "\n.. value-table::\n", XPRS_CLAMPING, values_clamping,0);

  AddSolverOption("alg:feastol feastol",
    "Primal feasibility tolerance (default 1e-6).",
    XPRS_FEASTOL, 0.0, DBL_MAX);

  AddSolverOption("alg:feastolperturb feastolperturb",
    "How much a feasible primal basic solution is allowed to "
    "be perturbed when performing basis changes.  The tolerance "
    "specified by \"alg:feastol\" is always considered as an upper "
    "limit for the perturbations; default = 1.0E-06",
    XPRS_FEASTOLPERTURB, 0.0, DBL_MAX);

  AddSolverOption("alg:feastoltarget feastoltarget",
    "Specifies the target feasibility tolerance for the solution refiner. "
    "Default = 0 (use the value of \"alg:feastol\")",
    XPRS_FEASTOLTARGET, 0.0, DBL_MAX);

    AddSolverOption("alg:indlinbigm indlinbigm",
                    "Largest \"big M\" value to use in converting indicator "
                    "constraints to regular constraints, default = 1e5",
                    XPRS_INDLINBIGM, 0.0, DBL_MAX);

    AddSolverOption("alg:lpfolding lpfolding",
                    "Simplex and barrier: whether to fold an LP problem before solving it:"
                    "\n.. value-table::\n",
                    XPRS_LPFOLDING, values_autonoyes_, -1);

    AddSolverOption("alg:maxiis maxiis",
                    "Maximum number of IIS to find; default=-1 (no limit)",
                    XPRS_MAXIIS, -1, INT_MAX);

    AddSolverOption("alg:zerotol matrixtol",
                    "The zero tolerance on matrix elements. If the value of a matrix element is less "
                    "than or equal to this in absolute value, it is treated as zero, default=1e-9.",
                    XPRS_MATRIXTOL, 0.0, DBL_MAX);

AddSolverOption("lp:pivtol pivtol markowitztol",
                "Markowitz pivot tolerance (default = 0.01)",
                XPRS_MARKOWITZTOL, 0.0, DBL_MAX);

AddSolverOption("alg:cutoff cutoff",
  "If the optimal objective value is worse than cutoff, "
  "report \"objective cutoff\" and do not return a solution. "
  "Default: 1.0E+40 for minimizing, -1.0E+40 for maximizing.",
  XPRS_MIPABSCUTOFF, MinusInfinity(), Infinity());

AddSolverOption("alg:addcutoff addcutoff mipaddcutoff",
  "Amount to add to the objective function of the best integer\n\
		solution found to give the new MIP cutoff; default -1e-5.",
  XPRS_MIPADDCUTOFF, -1e-10, DBL_MAX);

AddSolverOption("alg:relcutoff relcutoff miprelcutoff",
  "If the optimal objective value is (relatively) worse than relcutoff, "
  "report \"objective cutoff\" and do not return a solution. "
  "Default: 1.0E-4.",
  XPRS_MIPRELCUTOFF, 1e-4, Infinity());

AddSolverOption("alg:randomseed randomseed",
  "Sets the initial seed to use for the pseudo-random number generator in the "
  "Optimizer; default=1",
  XPRS_RANDOMSEED, -INT_MAX,  INT_MAX);

AddSolverOption("alg:refactor refactor",
  "Whether the optimization should restart using the current representation of the "
  "factorization in memory:\n"
  "\n.. value-table::\n",
  XPRS_REFACTOR, values_autonoyes_, -1);

AddSolverOption("alg:refineops refineops",
  "Bit vector: specifies wmhen the solution refiner should be executed to "
  "reduce solution infeasibilities. "
  "The refiner will attempt to satisfy the target tolerances for all original linear "
  "constraints before presolve or scaling has been applied:\n"
  "\n.. value-table::\n",
  XPRS_REFINEOPS, values_refineops, 19);

AddSolverOption("alg:resourcestrategy resourcestrategy",
  "Wether to allow nondeterministic decisions to cope with "
		"low memory (affected by maxmemory and maxmemoryhard):\n"
  "\n.. value-table::\n",
  XPRS_RESOURCESTRATEGY, values_01_noyes_0default_, 0);

//endalg
    // ****************************
  // PRESOLVER
  // ****************************
  AddSolverOption("pre:solve presolve",
    "Whether to use Xpress' presolve:\n"
    "\n.. value-table::\n",
    XPRS_PRESOLVE, presolve_values_, 1);

  AddSolverOption("pre:ops presolveops",
    "Reductions to use in XPRESS's presolve, sum of:\n"
    "\n.. value-table::\n(default 511 = bits 0-8 set)",
    XPRS_PRESOLVEOPS, presolveops_values_, 511);

  AddSolverOption("pre:maxscalefactor maxscalefactor",
  "Maximum log2 factor that can be applied during scaling, must be >=0 and <=64; default=64.",
  XPRS_MAXSCALEFACTOR, 0, 64);

  AddSolverOption("pre:elimfillin elimfillin",
    "Maximum fillins allowed for a presolve elimination; default = 10",
    XPRS_ELIMFILLIN, 0, INT_MAX);

  AddSolverOption("pre:elimtol elimtol",
    "The Markowitz tolerance for the elimination phase of the presolve; default=0.001",
    XPRS_ELIMTOL, 0.0, DBL_MAX);

  AddSolverOption("pre:genconsdualreductions genconsdualreductions",
    "Whether dual reductions should be applied to reduce the number "
    "of columns and rows added when transforming general constraints to MIP structs:\n"
    "\n.. value-table::\n",
    XPRS_GENCONSDUALREDUCTIONS, values_01_noyes_1default_, 1);

    AddSolverOption("pre:indlinbigm indprelinbigm",
            "Largest \"big M\" value to use in converting indicator "
            "constraints to regular constraints during XPRESS "
            "presolve; default = 100.0",
            XPRS_INDPRELINBIGM, 0.0, DBL_MAX);

    AddSolverOption("pre:maximpliedbound maximpliedbound",
                "When preprocessing MIP problems, only use computed bounds "
                "at most maximpliedbound (default 1e8) in absolute value",
                XPRS_MAXIMPLIEDBOUND, 0.0, DBL_MAX);


    AddSolverOption("pre:objscalefactor objscalefactor",
      "Power of 2 (default 0) by which the objective is scaled. "
      "Nonzero objscalfactor values override automatic global "
      "objective scaling",
      XPRS_OBJSCALEFACTOR, 0, INT_MAX);

    AddSolverOption("pre:basisred prebasisred",
      "Determines if a lattice basis reduction algorithm should be "
      "attempted as part of presolve:\n"
      "\n.. value-table::\n",
      XPRS_PREBASISRED, values_autonoyes_, -1);

    AddSolverOption("pre:bndredcone prebndredcone",
      "Determines if second order cone constraints should be used for "
      "inferring bound reductions on variables when solving a MIP:\n"
      "\n.. value-table::\n", 
      XPRS_PREBNDREDCONE, values_autonoyes_, -1);

    AddSolverOption("pre:bndredquad prebndredquad",
      "Determines if convex quadratic contraints should be used for "
      "inferring bound reductions on variables when solving a MIP",
      XPRS_PREBNDREDQUAD, values_autonoyes_, -1);

    AddSolverOption("pre:cliquestrategy precliquestrategy",
      "Determines how much effort to spend on clique covers "
      "in presolve; default=-1.",
      XPRS_PRECLIQUESTRATEGY, -1, INT_MAX);

    AddSolverOption("pre:coefelim precoefelim",
      "Specifies whether the optimizer should attempt to "
      "recombine constraints:\n"
      "\n.. value-table::\n",
      XPRS_PRECOEFELIM, values_precoefelim, 2);
    
    AddSolverOption("pre:components precomponents",
      "Determines whether small independent components should "
      "be detected and solved as individual subproblems during root "
      "node processing:\n"
      "\n.. value-table::\n",
      XPRS_PRECOMPONENTS, values_autonoyes_, -1);

    AddSolverOption("pre:componentseffort precomponentseffort",
      "adjusts the overall effort for the independent component "
      "presolver; default = 1.0.",
      XPRS_PRECOMPONENTSEFFORT, 0.0, DBL_MAX);

    AddSolverOption("pre:convertseparable preconvertseparable",
      "Reformulate problem with non-diagonal quadratic objective "
      "and/or constraints as diagonal quadratic or second-order conic "
      "constraints:\n"
      "\n.. value-table::\n",
      XPRS_PRECONVERTSEPARABLE, values_preconvertseparable, -1);

    AddSolverOption("pre:domcol predomcol",
      "Whether presolve should remove variables when solving MIP problems:\n"
      "\n.. value-table::\n",
      XPRS_PREDOMCOL, values_predomcol, -1);

    AddSolverOption("pre:domrow predomrow",
      "Whether presolve should remove constraints when solving MIP problems:\n"
      "\n.. value-table::\n",
      XPRS_PREDOMROW, values_predomrow, -1);

    AddSolverOption("pre:duprow preduprow",
      "How presolve should deal with duplicate rows in MIP problems:\n"
      "\n.. value-table::\n",
      XPRS_PREDUPROW, values_preduprow, -1);

    AddSolverOption("pre:elimquad preelimquad",
      "Allows for elimination of quadratic variables via doubleton rows:\n"
      "\n.. value-table::\n",
      XPRS_PREELIMQUAD, values_autonoyes_, -1);

    AddSolverOption("pre:folding prefolding",
      "Determines if a folding procedure should be used to aggregate "
      "continuous columns in an equitable partition:\n"
      "\n.. value-table::\n",
      XPRS_PREFOLDING, values_autonoyes_, -1);

    AddSolverOption("pre:implications preimplications",
      "Determines whether to use implication structures to remove redundant rows:\n"
      "\n.. value-table::\n",
      XPRS_PREIMPLICATIONS, values_autonoyes_, -1);

    AddSolverOption("pre:lindep prelindep",
      "Determines whether to check for and remove linearly dependent equality constraints when presolving a problem:\n"
      "\n.. value-table::\n",
      XPRS_PRELINDEP, values_autonoyes_, -1);

    AddSolverOption("pre:objcutdetect preobjcutdetect",
      "MIP: Determines whether to check for constraints that are "
      "parallel or near parallel to a linear objective function, "
      "and which can safely be removed:\n"
      "\n.. value-table::\n",
      XPRS_PREOBJCUTDETECT, values_01_noyes_1default_, 1);

    AddSolverOption("pre:permute prepermute",
      "Bit vector: specifies whether to randomly permute rows, columns "
      "and global information when starting the presolve (default 0):\n"
      "\n.. value-table::\n",
      XPRS_PREPERMUTE, values_prepermute, 0);

    AddSolverOption("pre:permuteseed prepermuteseed",
      "Sets the seed for the pseudo-random number generator for permuting; "
      "default=0",
      XPRS_PREPERMUTESEED, -INT_MAX, INT_MAX);

    AddSolverOption("pre:probing preprobing",
      "Amount of probing to perform on binary variables during presolve. "
      "This is done by fixing a binary to each of its values in turn and analyzing "
      "the implications:\n"
      "\n.. value-table::\n",
      XPRS_PREPROBING, values_predomrow, -1);

    AddSolverOption("pre:protectdual preprotectdual",
      "Specifies whether the presolver should protect a given dual solution "
      "by maintaining the same level of dual feasibility:"
      "\n.. value-table::\n",
      XPRS_PREPROTECTDUAL, values_01_noyes_0default_, 1);

    AddSolverOption("pre:maxgrow presolvemaxgrow",
      "Limit on how much the number of non-zero coefficients "
      "is allowed to grow during presolve, specified as a ratio of the "
      "number of non-zero coefficients in the original problem; default=0.1",
      XPRS_PRESOLVEMAXGROW, 0.0, DBL_MAX);

    AddSolverOption("pre:passes presolvepasses",
      "Number of reduction rounds to be performed in presolve; "
      "default=1.",
      XPRS_PRESOLVEPASSES, 0, INT_MAX);

    AddSolverOption("pre:pwldualreductions pwldualreductions",
      "Whether dual reductions should be applied to reduce the number "
      "of columns, rows and SOS-constraints added when transforming piecewise "
      "linear objectives and constraints to MIP structs:\n"
      "\n.. value-table::\n",
      XPRS_PWLDUALREDUCTIONS, values_01_noyes_1default_, 1);

    AddSolverOption("pre:pwlnonconvextransformation pwlnonconvextransformation",
      "Reformulation method for piecewise linear constraints at the "
      "beginning of the search:\n"
      "\n.. value-table::\n",
      XPRS_PWLNONCONVEXTRANSFORMATION, values_pwlnonconvextransformation, -1);

    AddSolverOption("pre:rootpresolve rootpresolve",
      "Whether to presolve after root cutting and heuristics:\n"
      "\n.. value-table::\n",
      XPRS_ROOTPRESOLVE, values_autonoyes_, -1);

    AddSolverOption("pre:scaling scaling",
      "Bit vector determining how to scale the constraint matrix before optimizing:\n"
      "\n.. value-table::\n",
      XPRS_SCALING, values_scaling, -1);

    AddSolverOption("pre:trace trace",
      "Display the infeasibility diagnosis during presolve:\n"
      "\n.. value-table::\n",
      XPRS_TRACE, values_01_noyes_0default_, 0);

    AddSolverOption("pre:sosreftol sosreftol",
      "Minimum relative gap between the ordering values of elements in a special "
      "ordered set; default=1e-6.",
      XPRS_SOSREFTOL, 1e-6, Infinity());

    //endpre
  // ****************************
  // BARRIER ALGORITHM CONTROLS
  // ****************************
  AddSolverOption("bar:alg baralg", "Which barrier algorithm to use ",
    XPRS_BARALG, values_baralg, -1);

  AddSolverOption("bar:cachesize cachesize",
    "Newton Barrier: L2 or L3 (see notes) cache size in kB (kilobytes) of the CPU (default -1). "
    "On Intel (or compatible) platforms a value of -1 may be used to determine the cache size automatically.",
    XPRS_CACHESIZE, -1, INT_MAX);

    AddSolverOption("bar:l1cache l1cache",
                    "Newton barrier: L1 cache size in kB (kilo bytes) of the CPU. On Intel (or compatible) "
                    "platforms a value of -1 may be used to determine the cache size automatically.",
                    XPRS_L1CACHE, -1, INT_MAX);

  AddSolverOption("bar:corespercpu corespercpu",
    "Newton Barrier: number of cores to assume per cpu. Barrier cache = cachesize/corespercpu. "
    " Default -1 = automatic. ",
    XPRS_CORESPERCPU, -1, INT_MAX);

  AddSolverOption("bar:cores barcores",
    "If positive, number of CPU cores to assume present when "
    "using the barrier algorithm.  Default = -1, which means "
    "automatic choice", XPRS_BARCORES, -1, INT_MAX);

  AddSolverOption("bar:choleskyalg choleskyalg",
    "Type of Cholesky factorization used for barrier, sum of:\n:",
    XPRS_CHOLESKYALG, values_barcholeskyalg, -1);

  AddSolverOption("bar:choleskytol choleskytol",
    "Zero tolerance for Cholesky pivots in the\n\
		Newton Barrier algorithm; default = 1e-15",
    XPRS_CHOLESKYTOL, 1e-15, DBL_MAX);

  AddSolverOption("bar:cpuplatform cpuplatform",
    "Type of Cholesky factorization used for barrier, sum of:\n:",
    XPRS_CPUPLATFORM, values_cpuplatform, -1);

  AddSolverOption("bar:crash barcrash",
    "Choice of crash procedure for crossover, higher number "
    "means more aggressive procedure:",
    XPRS_BARCRASH, values_barcrash, 4);

  AddSolverOption("bar:crossover crossover",
    "How to transform a barrier solution to a basic one:\n"
    "\n.. value-table::\n", XPRS_CROSSOVER, values_barcrossover, -1);

  AddSolverOption("bar:crossoverops crossoverops",
    "Bit vector affecting crossover after the barrier "
    "algorithm; sum of:\n"
    "\n.. value-table::\n", XPRS_CROSSOVEROPS, values_barcrossoverops, -1);

  AddSolverOption("bar:crossoverthreads crossoverthreads",
    "Limit on threads used during crossover; default -1 (determined by "
    "bar:threads).", XPRS_CROSSOVERTHREADS, -1, INT_MAX);

  AddSolverOption("lim:crossoveriterlim bar:crossoveriterlim crossoveriterlim crossoveritlim",
    "Limit on crossover iterations after the barrier "
    "algorithm; default = 2147483645", XPRS_CROSSOVERITERLIMIT, 1, INT_MAX);

  AddSolverOption("bar:crossovertol crossovertol crossoveraccuracytol",
    "Tolerance (default 1e-6) for deciding whether to adjust the "
    "relative pivot tolerance during crossover when a new basis "
    "factorization is necessary.  Errors in the recalculated "
    "basic solution above this tolerance cause the pivot "
    "tolerance to be adjusted.",
    XPRS_CROSSOVERACCURACYTOL, 0.0, DBL_MAX);

  AddSolverOption("bar:densecollimit densecollimit",
    "Number of nonzeros above which a column is treated as dense "
    "in the barrier algorithm's Cholesky factorization. Default=0 (automatic).",
    XPRS_DENSECOLLIMIT, 0, INT_MAX);

  AddSolverOption("bar:dualstop bardualstop",
    "Barrier method convergence tolerance on "
    "dual infeasibilities; default = 0 (automatic choice)",
    XPRS_BARDUALSTOP, 0.0, DBL_MAX);

  AddSolverOption("bar:gapstop bargapstop",
    "Barrier method convergence tolerance on the relative"
    " duality gap; default = 0", XPRS_BARGAPSTOP, 0.0, DBL_MAX);

  AddSolverOption("bar:gap bargaptarget", // todo DEFAULT?
    "Barrier algorithm target tolerance for the relative duality"
    " gap.If not satisfied and no further progress is possible"
    " but barstopgap is satisfied, then the current solution is"
    " considered optimal.", XPRS_BARGAPTARGET, 0.0, DBL_MAX);

  AddSolverOption("bar:indeflimit barindeflimit",
    "Maximum indefinite factorizations to allow in the barrier "
    "algorithm for solving a QP: stop when the limit is hit "
    "default = 15", XPRS_BARINDEFLIMIT, 0, INT_MAX);

  AddSolverOption("lim:bariterlim bar:iterlim bariterlim",
    "Limit on the number of barrier iterations (default 500).",
    XPRS_BARITERLIMIT, 1, INT_MAX);

  AddSolverOption("bar:start barstart",
    "Choice of starting point for barrier method:\n"
    "\n.. value-table::\n",
    XPRS_BARSTART, values_barstart, 0);

  AddSolverOption("bar:kernel barkernel",
    "How the barrier algorithm weights centrality:\n"
    "\n.. value-table::\n",
    XPRS_BARKERNEL, values_barkernel, 1.0);

  AddSolverOption("bar:objperturb barobjperturb",
    "Defines how the barrier perturbs the objective (default "
    "1e-6); values > 0 let the optimizer decide if to perturb the "
    "objective, values < 0 force the perturbation:\n"
    "\n.. value-table::\n",
    XPRS_BAROBJPERTURB, values_barobjperturb, 1e-6);

  AddSolverOption("bar:objscale barobjscale",
    "How the barrier algorithm scales the objective; when the objective "
    "is quadratic, the quadratic diagonal is used in determining the scale:\n"
    "\n.. value-table::\n",
    XPRS_BAROBJSCALE, values_barobjscale, -1);

  AddSolverOption("bar:order barorder",
    "Cholesky factorization pivot order for barrier algorithm:\n"
    "\n.. value-table::\n",
    XPRS_BARORDER, values_barorder, 0);

  AddSolverOption("bar:orderthreads barorderthreads",
    "Number of threads to use when choosing a pivot order for "
    "Cholesky factorization; default 0 (automatic choice).",
    XPRS_BARORDERTHREADS, 0, INT_MAX);

  AddSolverOption("bar:output baroutput",
    "Amount of output for the barrier method:\n"
    "\n.. value-table::\n",
    XPRS_BAROUTPUT, values_baroutput, 0);

  AddSolverOption("bar:presolve barpresolve",
    "Level of barrier-specific presolve effort:\n"
    "\n.. value-table::\n",
    XPRS_BARPRESOLVEOPS, values_baroutput, 0);

  AddSolverOption("bar:primalstop barprimalstop",
    "Barrier method convergence tolerance on "
    "primal infeasibilities; default = 0 (automatic choice)",
    XPRS_BARPRIMALSTOP, 0.0, DBL_MAX);

  AddSolverOption("bar:refiter barrefiter",
    "Maximum number of refinement iterations, helpful when the "
    "the solution is near to the optimum using barrier or crossover:\n"
    "\n.. value-table::\n",
    XPRS_BARREFITER, values_barrefiter, 0);

  AddSolverOption("bar:regularize barreg barrregularize",
    "Regularization to use with \"barrier. Default=-1 (automatic choice), "
    "else sum of:\n"
    "\n.. value-table::\n",
    XPRS_BARREGULARIZE, values_barregularize, -1);

  AddSolverOption("bar:stepstop barstepstop",
    "Barrier method convergence tolerance: stop when "
		"step size <= barstepstop; default = 1e-10", XPRS_BARSTEPSTOP, 1e-10, DBL_MAX);

  AddSolverOption("bar:threads threads",
    "number of threads used in the Newton Barrier algorithm;\n\
		default = -1 (determined by \"threads\")", XPRS_BARTHREADS, -1, INT_MAX);
  //endbarrier
  // ****************************
  // SIMPLEX RELATED
  // ****************************
  AddSolverOption("lp:bigmmethod bigmmethod",
    "Simplex: This specifies whether to use the \"Big M\" method, or the standard phase I "
    "(achieving feasibility) and phase II (achieving optimality). "
    "the \"Big M\" method, the objective coefficients of the variables are considered during "
    "the feasibility phase, possibly leading to an initial feasible basis which is closer to "
    "optimal. The side-effects involve possible round-off errors.\n"
    "\n.. value-table::\n",
    XPRS_BIGMMETHOD, values_bigmmethod, 1);

  AddSolverOption("lp:bigm bigm bigmpenalty",
    "Infeasibility penalty to be used if \"BigM\" method is used; default = 1024",
    XPRS_BIGM, 1024.0, DBL_MAX);

  AddSolverOption("lp:crash crash",
    "Simplex: This determines the type of crash used when the algorithm begins."
    "For primal simplex, the choices are listed below; for dual simplex the choices "
    "follow and are interpreted a bit-vector:\n"
    "\n.. value-table::\n",
    XPRS_CRASH, values_crash, 2);

  AddSolverOption("lp:dualgradient dualgradient",
    "dual simplex pricing strategy:\n"
    "\n.. value-table::\n",
    XPRS_CRASH, values_lpdualgradient, -1);

  AddSolverOption("lp:dualize dualize",
    "Whether to convert the primal problem to its dual and solve "
    "the converted problem:\n"
    "\n.. value-table::\n",
    XPRS_DUALIZE, values_lpdualize, -1);

  AddSolverOption("lp:dualstrategy dualstrategy",
    "Bit vector controlling the dual simplex strategy (default 1):\n"
    "\n.. value-table::\n",
    XPRS_DUALSTRATEGY, values_lpdualstrategy, 1);

  AddSolverOption("lp:dualthreads dualthreads ",
    "Limit on threads used by parallel dual simplex; default -1 (determined by "
    "tech:threads).", XPRS_DUALTHREADS, -1, INT_MAX);

  AddSolverOption("lp:dualizeops dualizeops",
    "When solving the dual problem after deriving it from the "
    "primal, whether to use primal simplex if dual simplex was "
    "specified and vice versa:\n"
    "\n.. value-table::\n",
    XPRS_DUALIZEOPS, values_01_noyes_1default_, 1);

  AddSolverOption("lp:dualperturb dualperturb",
    "Factor by which the problem will be perturbed prior "
    "to optimization by dual simplex. Default -1 (automatic); note that "
    "a value of 0 implies no perturbation",
    XPRS_DUALPERTURB, -1.0, DBL_MAX);

  AddSolverOption("lp:etatol etatol",
    "Zero tolerance on eta elements, elements of eta vectors whose "
    "absolute value is smaller than etatol are taken to be zero.",
    XPRS_ETATOL, 0.0, DBL_MAX);

  AddSolverOption("lp:dualforceparallel forceparalleldual dualforceparallel",
    "Specifies whether the dual simplex solver should always use the "
    "parallel simplex algorithm",
    XPRS_FORCEPARALLELDUAL, values_01_noyes_0default_, 0);



    AddSolverOption("lp:invertfreq invertfreq",
                    "Maximum simplex iterations before refactoring the basis; "
                    "default -1 (automatic)",
                    XPRS_INVERTFREQ, -1, INT_MAX);

    AddSolverOption("lp:invertmin invertmin",
                    "Minimum simplex iterations before refactoring the basis; "
                    "default = 3",
                    XPRS_INVERTMIN, 3, INT_MAX);

    AddSolverOption("lp:keepbasis keepbasis",
                    "Basis choice for the next LP iteration:\n"
                    "\n.. value-table::\n",
                    XPRS_KEEPBASIS, values_keepbasis, 1);

    AddSolverOption("lp:keepnrows keepnrows",
                    "Status for nonbinding rows:\n"
                    "\n.. value-table::\n",
                    XPRS_KEEPNROWS, values_keepnrows, 1);

    AddSolverOption("lp:log lplog",
                    "Frequency of printing simplex iteration log; default = 100."
                    "Values n < 0 display detailed outputs every -n iterations.",
                    XPRS_LPLOG, -INT_MAX, INT_MAX);

    AddSolverOption("lp:netstalllimit netstalllimit",
      "Limit the number of degenerate pivots of the network "
      "simplex algorithmm before switching to primal or dual:\n"
      "\n.. value-table::\n",
      XPRS_NETSTALLLIMIT, values_netstalllimit, -1);

   
    AddSolverOption("lp:optimalitytol optimalitytol",
      "This is the zero tolerance for reduced costs. On each "
      "iteration, the simplex method searches for a variable to enter "
      "the basis which has a negative reduced cost. The candidates are "
      "only those variables which have reduced costs less than the "
      "negative value of optimalitytol; default=1e-6",
      XPRS_OPTIMALITYTOL, 0.0, DBL_MAX);


    AddSolverOption("lp:optimalitytoltarget optimalitytoltarget",
      "Target optimality tolerance for the solution refiner; default=0 "
      "(use the value specified by lp:optimalitytol)",
      XPRS_OPTIMALITYTOLTARGET, 0.0, DBL_MAX);

    AddSolverOption("lp:penalty penalty",
      "Minimum absolute penalty variable coefficient; "
	    "default = automatic choice",
      XPRS_PENALTY, 0.0, DBL_MAX);

    AddSolverOption("lp:pricingalg pricingalg",
      "Primal simplex pricing method:\n"
      "\n.. value-table::\n",
      XPRS_PRICINGALG, values_pricingalg, 0);

    AddSolverOption("lp:primalunshift primalunshift",
      "Whether the primal alg. calls the dual to unshift:\n"
      "\n.. value-table::\n",
      XPRS_PRICINGALG, values_01_noyes_1default_, 0);

    AddSolverOption("lp:relpivottol relpivottol",
      "Relative pivot tolerance; default = 1e-6",
      XPRS_RELPIVOTTOL, 0.0, Infinity());

    AddSolverOption("lp:sifting sifting",
      "When using dual simplex, whether to enable sifting, "
		"which can speed up the solve when there are many more "
		"variables than constraints:\n"
      "\n.. value-table::\n",
      XPRS_SIFTING, values_autonoyes_, -1);

    AddSolverOption("lp:siftpasses siftpasses",
      "Determines how quickly we allow to grow the worker problems "
      "during the sifting algorithm; default 4.",
      XPRS_SIFTPASSES, 1, INT_MAX);

    AddSolverOption("lp:siftpresolveops siftpresolveops",
      "Presolve operations for solving the subproblems during sifting:\n"
      "\n.. value-table::\n",
      XPRS_SIFTPRESOLVEOPS, values_siftpresolveops, -1);

    AddSolverOption("lp:siftswitch siftswitch",
      "Determines which algorithm to use for solving the subproblems during sifting:\n"
      "\n.. value-table::\n",
      XPRS_SIFTSWITCH, values_siftswitch, -1);


    // endlp
  // MIP
  // ****************************
  AddSolverOption("mip:branchchoice branchchoice",
    "Control the choice of branching when solving a MIP problem:\n"
    "\n.. value-table::\n",
    XPRS_BRANCHCHOICE, values_branchchoice, 3);

  AddSolverOption("mip:branchdisj branchdisj",
    "Whether to branch on general split disjunctions while solving MIPs:\n"
    "\n.. value-table::\n",
    XPRS_BRANCHDISJ, values_branchdisj, -1);

  AddSolverOption("mip:branchstructural branchstructural branchstruct",
    "Whether to search for special structure during branch and bound:\n"
    "\n.. value-table::\n",
    XPRS_BRANCHSTRUCTURAL, values_autonoyes_, -1);

  AddSolverOption("mip:breadthfirst breadthfirst",
    "Number of MIP nodes included in best-first search before switching to local-first search; "
    "default=11.", XPRS_BREADTHFIRST, 11, INT_MAX);
  
  AddSolverOption("mip:deterministic deterministic",
    "Whether a MIP search should be deterministic:\n", 
    XPRS_DETERMINISTIC, values_deterministic, 1);

  AddSolverOption("mip:feasibilitypump feasibilitypump",
    "Decides whether to run the Feasibility Pump heuristic at the top "
    "node during branch-and-bound:\n"
    "\n.. value-table::\n", XPRS_FEASIBILITYPUMP,
    values_feasibilitypump, -1);

    AddSolverOption("mip:heurbeforelp heurbeforelp",
    "Whether primal heuristics should be run before the initial LP relaxation has been solved:\n"
    "\n.. value-table::\n",
    XPRS_HEURBEFORELP, values_autonoyes_, -1);

    AddSolverOption("lim:heurdiveiterlimit heurdepth mip:heurdiveiterlimit",
      "Simplex iteration limit for reoptimizing during the diving heuristic; "
      "default = -1 (automatic selection); a value of 0 implies no iteration limit",
      XPRS_HEURDIVEITERLIMIT, -INT_MAX, INT_MAX);

    AddSolverOption("mip:heurdiverandomize hdive_rand heurdiverandomize",
      "The level of randomization to apply in the diving heuristic; values "
      "range from 0.0=none to 1.0=full.",
      XPRS_HEURDIVERANDOMIZE, 0.0, 1.0);

    AddSolverOption("mip:heurdivesoftrounding hdive_rounding heurdivesoftrounding",
      "Whether to use soft rounding in the MIP diving heuristic "
      "(to push variables to their bounds via the objective rather "
      "than fixing them):\n"
      "\n.. value-table::\n",
      XPRS_HEURDIVESOFTROUNDING, values_heurdivesoftrounding, -1);

    AddSolverOption("mip:heurdivespeedup hdive_speed heurdivespeedup",
      "Controls tradeoff between speed and solution quality in the diving heuristic:"
      "\n.. value-table::\n",
      XPRS_HEURDIVESPEEDUP, values_heurdivespeed, -1);

    AddSolverOption("mip:heurdivestrategy hdive_strategy heurdivestrategy",
      "Chooses the strategy for the diving heuristic:\n"
      "\n.. value-table::\n",
      XPRS_HEURDIVESTRATEGY, values_heurdivestrategy, -1);

    AddSolverOption("mip:heuremphasis heuremphasis",
                    "Chooses the strategy for the diving heuristic:\n"
                    "\n.. value-table::\n",
                    XPRS_HEUREMPHASIS, values_heuremphasis, -1);

    AddSolverOption("mip:heurforcespecialobj heurforcespecobj heurforcespecialobj" ,
                    "Whether to use special objective heuristics on large problems and even if an incumbant exists:\n"
                    "\n.. value-table::\n",XPRS_HEURFORCESPECIALOBJ, values_01_noyes_0default_, 0);

    AddSolverOption("mip:heurfreq heurfreq",
                    "During branch and bound, heuristics are applied at nodes whose depth "
                    "from the root is zero modulo \"heurfreq\"; default -1 (automatic).",
                    XPRS_HEURFREQ, -1, INT_MAX);

    AddSolverOption("mip:heursearcheffort heursearcheffort",
                    "Adjusts the overall level of the local search heuristics; default 1.0 (normal level).",
                    XPRS_HEURSEARCHEFFORT, 1.0, DBL_MAX);

    AddSolverOption("mip:heursearchfreq heurfreq heursearchfreq",
                    "Specifies how often the local search heuristic should be run in the tree:\n"
                    "\n.. value-table::\n",
                    XPRS_HEURSEARCHFREQ, values_heursearchfreq,-1);

    AddSolverOption("mip:heursearchrootcutfreq heurrootcutfreq heursearchrootcutfreq",
                    "How often to run the local search heuristic while cutting at the root node:\n"
                    "\n.. value-table::\n",
                    XPRS_HEURSEARCHROOTCUTFREQ, values_heursearchrootcutfreq,-1);


    AddSolverOption("mip:heursearchrootselect  heursearchrootselect",
                    "A bit vector control for selecting which local search heuristics to apply on the root node of a global solve; "
                    "default 117:\n"
                    "\n.. value-table::\n",
                    XPRS_HEURSEARCHROOTSELECT, values_heursearchrootcutselect,-1);

    AddSolverOption("mip:heursearchtreeselect  heursearchtreeselect",
                    "A bit vector control for selecting which local search heuristics to apply during the tree search of a global solve, "
                    "default 17:\n"
                    "\n.. value-table::\n",
                    XPRS_HEURSEARCHTREESELECT, values_heursearchrootcutselect,-1);

    AddSolverOption("mip:heurthreads heurtreads",
                    "Number of threads to dedicate to running heuristics on the root node:\n",
                    XPRS_HEURTHREADS, values_heurthreads, 0);

    AddSolverOption("mip:historycosts historycosts",
                    "How to update the pseudo cost for a global entity when a strong branch or a regular branch is applied:\n",
                    XPRS_HISTORYCOSTS, values_historycosts, -1);

    AddSolverOption("mip:localchoice localchoice",
                    "when to backtrack between two child nodes during a \"dive\":\n"
                    "\n.. value-table::\n",
                    XPRS_LOCALCHOICE, values_localchoice, 1);

AddSolverOption("mip:maxlocalbacktrack maxlocalbacktrack maxlocalbt",
                "Max height above current node to look for a local backtrack "
                "candidate node; default=-1(automatic)",
                XPRS_MAXLOCALBACKTRACK, -1, INT_MAX);

AddSolverOption("mip:maxtasks maxmiptasks",
                "Maximum tasks to run in parallel during a MIP solve; default = -1 "
                "(use mip:threads)."
                "For mip:maxtasks > 0, branch-and-bound nodes are solved in a "
		        "deterministic way, but the barrier algorithm (if used) may "
                "cause a nondeterministic MIP solve unless bar:threads = 1.",
                XPRS_MAXMIPTASKS, -1, INT_MAX);

AddSolverOption("mip:gap mipgap",
  "Max. relative MIP optimality gap (default 1e-4).",
  XPRS_MIPRELSTOP, 1e-4, DBL_MAX);

AddSolverOption("mip:gapabs mipgapabs",
  "Max. absolute MIP optimality gap (default 0).",
  XPRS_MIPABSSTOP, 0.0, DBL_MAX);

AddSolverOption("mip:components mipcomponents",
  "Determines whether disconnected components in a MIP should\n\
		be solved as separate MIPs:\n"
  "\n.. value-table::\n",
  XPRS_MIPCOMPONENTS, values_autonoyes_, -1);

AddSolverOption("mip:concurrentnodes mipconcurrentnodes",
  "Node limit to choose the winning solve when concurrent\n\
		solves are enabled:\n"
  "\n.. value-table::\n",
  XPRS_MIPCONCURRENTNODES, values_mipconcurrentnodes, -1);

AddSolverOption("mip:concurrentsolves mipconcurrentsolves",
  "Select the number of concurrent solves to start for a MIP:\n"
  "\n.. value-table::\n",
  XPRS_MIPCONCURRENTSOLVES, values_mipconcurrentsolves, 0);

AddSolverOption("mip:dualreductions mipdualreductions",
  "Kinds of dual reductions allowed during branch and bound:\n"
  "\n.. value-table::\n"
  "\nIf poolnbest > 1 is specified, specifying "
  "mipdualreductions = 2 might be prudent.",
  XPRS_MIPDUALREDUCTIONS, values_mipdualreductions, 1);


AddSolverOption("mip:kappafreq mipkappafreq",
  "During branch-and-bound, how often to compute "
  "basis condition numbers:\n"
  "\n.. value-table::\n"
  "\nWhen mipkappafreq > 0, a final summary shows the number of "
  "sampled nodes that are:\n\
			\"stable\": kappa < 10^7\n\
			\"suspicious\": 10^7 <= kappa < 10^10\n\
			\"unstable\": 10^10 <= kappa < 10^13\n\
			\"ill-posed\": 10^13 <= kappa.\n"
    "A \"Kappa attention level\" between 0 and 1 is also reported. "
    "Condition numbers use the Frobenius norms of the basis "
    "and its inverse.", 
  XPRS_MIPKAPPAFREQ, values_mipkappafreq,0);

AddSolverOption("mip:log miplog",
  "Frequency of printing MIP iteration log; default = -100."
  "Values n < 0 display detailed outputs every -n iterations.",
  XPRS_LPLOG, -INT_MAX, INT_MAX);

AddSolverOption("mip:presolve mippresolve",
  "Type of integer processing to be performed. "
  "If set to 0, no processing will be performed (default automatic):\n"
  "\n.. value-table::\n",
  XPRS_MIPPRESOLVE, values_mippresolve, 1);

AddSolverOption("mip:rampup miprampup",
  "Whether to limit the number of parallel tasks\n\
		during the ramp-up phase of the parallel MIP algorithm:\n"
  "\n.. value-table::\n",
  XPRS_MIPRAMPUP, values_miprampup, -1);

AddSolverOption("mip:miprefineiterlimit miprefiterlim miprefineiterlimit",
  "Max. simplex iterations per reoptimization in MIP refiner "
  "when refineops is 2 or 3; default -1 (automatic).",
  XPRS_MIPREFINEITERLIMIT, -1, INT_MAX);

AddSolverOption("mip:restart miprestart",
  "Control strategy for in-tree restarts:\n"
  "\n.. value-table::\n",
  XPRS_MIPRESTART, values_miprestart, -1);

AddSolverOption("mip:restartgapthreshold miprestartgapthreshold",
  "Initial gap threshold to delay in-tree restart; "
    "the restart is delayed if the relative gap is below the "
    "threshold (default 0.02)",
  XPRS_MIPRESTARTGAPTHRESHOLD, 0.0, DBL_MAX);

AddSolverOption("mip:restartfactor miprestartfactor",
  "Fine tune initial conditions to trigger an in-tree "
  "restart; values > 1 increase the aggressiveness, < 1 "
  "decrease it (default 1.0)",
  XPRS_MIPRESTARTFACTOR, 0.0, DBL_MAX);

AddSolverOption("mip:threads mipthreads",
  "Determines the number of threads implemented to run the parallel "
  "MIP code; default -1: alg:threads will determine the number of threads.",
  XPRS_MIPTHREADS, -1, INT_MAX);

AddSolverOption("mip:intfeastol intfeastol",
  "Feasibility tolerance for integer variables (default 5e-06).",
  XPRS_MIPTOL, 0.0, DBL_MAX);

AddSolverOption("mip:toltarget miptoltarget",
  "Value of miptol used for refining equalities on MIP "
		"problems when refineops is 2 or 3; default = 0.",
  XPRS_MIPTOLTARGET, 0.0, DBL_MAX);

AddSolverOption("mip:nodeprobingeffort nodeprobingeffort",
  "Multiplier on the default amount of work node probing "
  "should do. Setting the control to zero disables node probing.",
  XPRS_NODEPROBINGEFFORT, 0.0, DBL_MAX);

AddSolverOption("mip:nodeselection nodeselection",
  "Determines which nodes will be considered for solution "
  "once the current node has been solved:\n"
  "\n.. value-table::\n",
  XPRS_NODESELECTION, values_nodeselection, 0);

AddSolverOption("mip:pseudocost pseudocost",
  "Default pseudo-cost assumed for forcing an integer variable\n\
		to an integer value; default = 0.01",
  XPRS_PSEUDOCOST, 0.0, DBL_MAX);



AddSolverOption("mip:qcrootalg qcrootalg",
  "when using miqcpalg = 1 to solve a mixed - integer problem that "
  "has quadratic constraints or second - order cone constraints, "
  "the algorithm for solving the root node:\n"
  "\n.. value-table::\n",
  XPRS_QCROOTALG, values_qcrootalg, -1);

AddSolverOption("mip:relaxtreememorylimit relaxtreemem relaxtreememorylimit",
  "Fraction of memory limit by which to relax \"treememlimit\" "
    "when too much structural data appears; default 0.1. Set to 0 to never relax "
  "the memory limit in this way.",
  XPRS_RELAXTREEMEMORYLIMIT, 0.0, DBL_MAX);

AddSolverOption("mip:sbbest sbbest",
  "Number of infeasible global entities to initialize pseudo costs for on each node:\n"
  "\n.. value-table::\n",
  XPRS_SBBEST, values_sbbest, -1);

AddSolverOption("mip:sbeffort sbeffort",
  "Adjusts the overall amount of effort when using strong branching to select an "
  "infeasible global entity to branch on; default = 1.",
  XPRS_SBEFFORT, 0.0, DBL_MAX);

AddSolverOption("mip:sbestimate sbestimate",
  "How to compute pseudo costs from the local node "
  "when selecting an infeasible entity to branch on:\n"
  "\n.. value-table::\n",
  XPRS_SBESTIMATE, values_sbestimate, -1);

AddSolverOption("mip:sbiterlimit sbiterlimit",
  "Number of dual iterations to perform the strong branching; "
		"0=none, default = -1 (automatic choice)",
  XPRS_SBITERLIMIT, -1, INT_MAX);

AddSolverOption("mip:sbselect sbselect",
  "size of candidate list for strong branching:\n"
  "\n.. value-table::\n",
  XPRS_SBSELECT, values_sbselect, -2);

AddSolverOption("mip:symmetry symmetry",
  "Amount of effort to detect symmetry in MIP problems:\n"
  "\n.. value-table::\n",
  XPRS_SYMMETRY, values_symmetry, 1);

AddSolverOption("mip:symselect symselect",
  "Adjusts the overall amount of effort for symmetry detection:\n"
  "\n.. value-table::\n",
  XPRS_SYMSELECT, values_symselect, -1);

AddSolverOption("mip:varselection varselection",
  "How to score the integer variables at a MIP node, for "
    "branching on a variable with minimum score:\n"
    "\n.. value-table::\n",
  XPRS_VARSELECTION, values_varselection, -1);


// endmip
    // ****************************
  // Cuts
  // ****************************
  AddSolverOption("cut:cover covercuts",
    "The number of rounds of lifted cover inequalities at the top node."
    "Default=-1, automatic.",
    XPRS_COVERCUTS, -1, INT_MAX);

  AddSolverOption("cut:gomory gomcuts",
    "The number of rounds of Gomory or lift-and-project cuts at the top node."
    "Default=-1, automatic.",
    XPRS_GOMCUTS, -1, INT_MAX);

  AddSolverOption("cut:qccuts qccuts",
    "when using miqcpalg=1 to solve a mixed-integer problem that "
    "has quadratic constraints or second-order cone constraints, "
    "the number of rounds of outer approximation cuts at the top "
    "node; default = -1 (automatic choice).",
    XPRS_QCCUTS, -1, INT_MAX);

    AddSolverOption("cut:lnpbest lnpbest",
                    "Number of infeasible global entities to create lift-and-project cuts "
                    "for during each round of Gomory cuts at the top node",
                    XPRS_LNPBEST, 1, INT_MAX);

    AddSolverOption("cut:lnpiterlimit lnpiterlimit",
                    "Number of iterations to perform in improving each lift-and-project cut; "
                    "default=-1 (automatic)",
                    XPRS_LNPITERLIMIT, -1, INT_MAX);

  AddSolverOption("cut:treecover treecovercuts",
    "The number of rounds of lifted cover inequalities at MIP nodes "
    "other than the top node. Default=-1 (automatic).",
    XPRS_TREECOVERCUTS, -1, INT_MAX);

  AddSolverOption("cut:treegomory treegomcuts",
    "The number of rounds of Gomory or lift-and-project cuts at MIP nodes "
    "other than the top node. Default=-1 (automatic).",
    XPRS_TREEGOMCUTS, -1, INT_MAX);

  AddSolverOption("cut:treeqccuts treeqccuts",
    "when using miqcpalg=1 to solve a MIP that "
    "has quadratic constraints or second-order cone constraints, "
    "the number of rounds of outer approximation cuts during the "
    "tree search; default = -1 (automatic choice).",
    XPRS_TREEQCCUTS, -1, INT_MAX);

  AddSolverOption("cut:depth cutdepth",
    "Maximum MIP tree depth at which to generate cuts. Default "
    "-1 (automatic); a value of 0 will disable cuts generation.",
    XPRS_CUTDEPTH, -1, INT_MAX);


  AddSolverOption("cut:factor cutfactor",
    "Limit on number of cuts and cut coefficients added "
    "while solving MIPs. Default=-1 (automatic); a value of 0 "
    "will disable cuts generation.",
    XPRS_CUTFACTOR, -1, INT_MAX);

  AddSolverOption("cut:freq cutfreq",
    "Cuts are only generated at tree depths that are integer\n\
		multiples of cutfreq. Default=-1 (automatic choice).",
    XPRS_CUTFREQ, -1, INT_MAX);

  AddSolverOption("cut:select cutselect",
    "Detailed control of cuts at MIP root node; sum of:\n"
    "\n.. value-table::\n",
    XPRS_CUTSELECT, values_cutselect, -1);

  AddSolverOption("cut:treeselect treecutselect",
    "Detailed control of cuts created during the tree search; sum of:\n"
    "\n.. value-table::\n",
    XPRS_TREECUTSELECT, values_cutselect, -1);

  AddSolverOption("cut:strategy cutstrategy",
    "How aggressively to generate MIP cuts; more ==> fewer nodes "
    "but more time per node:\n"
    "\n.. value-table::\n",
    XPRS_CUTSTRATEGY, values_cutstrategy, -1);

  // ****************************
  // Q(C)P
  // ****************************
  AddSolverOption("qp:nonconvex nonconvex",
    "Determines if the convexity of the problem is checked before optimization:\n"
    "\n.. value-table::\n",
    XPRS_IFCHECKCONVEXITY, values_01_noyes_1default_, 1);

  AddSolverOption("qp:eigenvaluetol eigenvaluetol",
    "Regard the matrix in a quadratic form as indefinite if its "
    "smallest eigvenalue is < -eigevnaltol; default = 1e-6",
    XPRS_EIGENVALUETOL, 0.0, DBL_MAX);

  AddSolverOption("qp:simplexops qsimplexops",
    "Bit vector, controls the behavior of the quadratic simplex solvers:\n"
    "\n.. value-table::\n",
    XPRS_QSIMPLEXOPS, values_qsimplexops, 1);

  AddSolverOption("qp:miqcpalg miqcpalg",
    "Which algorithm is to be used to solve mixed integer quadratic constrained and mixed integer second order cone problems:\n"
    "\n.. value-table::\n",
    XPRS_MIQCPALG, values_miqcpalg, 1);

  AddSolverOption("qp:unshift quadunshift quadraticunshift",
    "whether quadratic simplex should do an extra "
	  "purification after finding a solution:\n"
    "\n.. value-table::\n",
    XPRS_QUADRATICUNSHIFT, values_autonoyes_, -1);

  AddSolverOption("qp:repairindefiniteq repairindefq repairindefiniteq",
    "whether to repair indefinite quadratic forms:\n"
    "\n.. value-table::\n",
    XPRS_REPAIRINDEFINITEQ, values_yesnoinverted_defaultno, -1);

  //endqp



  AddStoredOption("tech:tunebase tunerdir tunebase",
    "Base name for results of running XPRESS's search for best "
    "parameter settings. The search is run only when tunebase "
    "is specified.  This control only defines the root path for the tuner "
    "output. For each problem, the tuner result will be output to a subfolder "
    "underneath this path. For example, by default, the tuner result for a "
    "problem called prob will be located at tuneroutput/prob/",
    storedOptions_.tunebase_);


  AddStoredOption("tech:tunename tunesessionname",
    "Set problem name within the tuner \"tunebase\" is specified.",
    storedOptions_.tunename_);


  AddSolverOption("tech:tuneoutput tuneroutput tuneoutput",
    "Output tuner results and logs to the file system when \"tunebase\" is specified:\n"
    "\n.. value-table::\n",
    XPRS_TUNEROUTPUT, values_01_noyes_1default_, 2);

  AddSolverOption("tech:tunerhistory tunerhistory",
    "Reuse and append to previous tuner results of the same problem:\n"
    "\n.. value-table::\n",
    XPRS_TUNERHISTORY, values_tunerhistory, 2);

  AddSolverOption("tech:tunetimelim tunermaxtime tunetimelim lim:tunetime",
    "Time limit (in seconds) on tuning when \"tunebase\" "
    "is specified; default 0 (no time limit).",
    XPRS_TUNERMAXTIME, 0, INT_MAX);

  AddSolverOption("tech:tunerthreads tunerthreads",
    "Number of tuner threads to run in parallel; "    
    "default=-1 (automatic)",
    XPRS_TUNERTHREADS, -1, INT_MAX);

  AddSolverOption("tech:tunermethod tunermethod",
    "Method for tuning when \"tunebase\" is specified:\n"
    "\n.. value-table::\n",
    XPRS_TUNERMETHOD, values_tunermethod, -1);

  AddSolverOption("tech:tunertarget tunertarget",
    "What to measure to compare two problem solutions "
		"when running the XPRESS tuner:\n"
    "\n.. value-table::\n",
    XPRS_TUNERTARGET, values_tunertarget, -1);

  AddSolverOption("tech:tunerverbose tunerverbose",
    "whether the tuner should prints detailed information for each run:\n"
    "\n.. value-table::\n",
    XPRS_TUNERVERBOSE, values_01_noyes_1default_, 1);

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
  int status = XPRSgetbasis(lp(), NULL, vars.data());
  if (status)
    vars.clear();
  else 
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
  int status = XPRSgetbasis(lp(), cons.data(), NULL);
  if (status)
    cons.clear();
  else
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

std::vector<int> XpressmpBackend::VarStatii(ArrayRef<int> vst) {
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
        XPRESSMP_CCALL(XPRSgetlb(lp(), lb.data(), 0, vst.size()-1));
        XPRESSMP_CCALL(XPRSgetub(lp(), ub.data(), 0, vst.size()-1));
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
  return stt;
}

std::vector<int> XpressmpBackend::ConStatii(ArrayRef<int> cst) {
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
  return stt;
}

SolutionBasis XpressmpBackend::GetBasis() {
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
  return { varstt,constt};

}

void XpressmpBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  auto convertedVarBasis = VarStatii(varstt);
  auto convertedConBasis =ConStatii(constt);
  XPRESSMP_CCALL(XPRSloadbasis(lp(), convertedConBasis.data(), convertedVarBasis.data()));

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

void XpressmpBackend::AddPrimalDualStart(Solution sol0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution(
    { sol0_unpres.primal, sol0_unpres.dual });
  auto x0 = mv.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);

  int status;
  XPRESSMP_CCALL(XPRSloadlpsol(lp(), x0.data(), NULL,
    pi0.data(), NULL, &status));
  if (status)
    fmt::print("warmstart: solution is not loaded because the problem is in presolved status.\n");
}

void XpressmpBackend::AddMIPStart(ArrayRef<double> x0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution({ x0_unpres });
  auto x0 = mv.GetVarValues()();
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
  if (outlev_ == 0)
    return;
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
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::XpressmpBackend()},
    slv_opt, cb);
}

void AMPLSCloseXpressmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

XPRSprob GetXpressmpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::XpressmpBackend*>(AMPLSGetBackend(slv))->lp();
}

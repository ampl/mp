#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mosekbackend.h"

extern "C" {
  #include "mosek-ampls-c-api.h"    // Mosek AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptMosek(void* prob) {
  // TODO
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateMosekBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::MosekBackend()};
}


namespace mp {

/// Create Mosek Model Manager
/// @param gc: the Mosek common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return MosekModelMgr
std::unique_ptr<BasicModelManager>
CreateMosekModelMgr(MosekCommon&, Env&, pre::BasicValuePresolver*&);


MosekBackend::MosekBackend() {
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateMosekModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

MosekBackend::~MosekBackend() {
  CloseSolver();
}

static void MSKAPI printstr(void* handle, const char* str)
{
  fmt::print("{}", str);
  fflush(stdout);
}

void MosekBackend::OpenSolver() {
  int status = MSK_RES_OK;
  const auto& initialize = GetCallbacks().init;

  MSKenv_t env = NULL;
  MSKtask_t task;
  if (initialize) { // If an initialization callback is provided,
    // use it to create the environment
    MOSEK_CCALL(MSK_makeenv(&env, NULL));
    env = (MSKenv_t)initialize();
    set_env(env);
  }
  else {
    MSKrescodee res = MSK_makeenv(&env, NULL);
    set_env(env);
    // Check license here
    if (MSK_checkoutlicense(env, MSK_FEATURE_PTS))
    {
      const auto diag = GetCallbacks().diagnostics;
      if (diag) {
        // If a diagnostic function is provided, do not print more
        // information
        diag();
        exit(res);
      }
    }
  }
  status = MSK_maketask(env, 0, 0, &task);
  if (status)
    throw std::runtime_error(fmt::format(
      "Failed to create task, error code {}.", status));
  set_lp(task); // Assign it

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();

  /// Turn off verbosity by default
  MOSEK_CCALL(MSK_putintparam(task, MSK_IPAR_LOG, 0));
  
  // Register callback for console logging (controlled by the outlev param
  // in all AMPL solver drivers)
  MSK_linkfunctotaskstream(lp(), MSK_STREAM_LOG, NULL, printstr);

}

void MosekBackend::CloseSolver() {
  if (lp()) MSK_deletetask(&lp_ref());
  if (env()) MSK_deleteenv(&env_ref());
}

const char* MosekBackend::GetBackendName()
  { return "MosekBackend"; }

std::string MosekBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", MSK_VERSION_MAJOR,
    MSK_VERSION_MINOR, MSK_VERSION_REVISION);
}

void MosekBackend::InitOptionParsing() {
  OpenSolver();
}

bool MosekBackend::IsMIP() const {
  int is_int_sol = 0;
  MOSEK_CCALL(MSK_solutiondef(lp(), MSK_SOL_ITG, &is_int_sol));
  return is_int_sol;
}

bool MosekBackend::IsQCP() const {
  // TODO
// return getIntAttr(MOSEK_INTATTR_QELEMS) > 0;
  return false;
}

ArrayRef<double> MosekBackend::PrimalSolution() {
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  // TODO get appropriate solution
  MSK_getxx(lp(), solToFetch_, x.data());
  return x;
}

pre::ValueMapDbl MosekBackend::DualSolution() {
	return {{
			{ CG_Algebraic, DualSolution_LP() }
		}};
}

ArrayRef<double> MosekBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  // TODO get appropriate solution
  MSKrescodee error = MSK_gety(lp(), solToFetch_, pi.data());
  if (error != MSK_RES_OK)
    pi.clear();
  return pi;
}

ArrayRef<double> MosekBackend::DualSolution_Cones() {
	// TODO can we return anything sensible?
	// Variable suffixes?
	return {};
}

double MosekBackend::ObjectiveValue() const {
  double v;
  // TODO get appropriate solution
  MOSEK_CCALL(MSK_getprimalobj(lp(), solToFetch_, &v));
  return v;
}

double MosekBackend::NodeCount() const {
  return getIntAttr(MSK_IINF_MIO_NUM_ACTIVE_NODES); // TODO check
}

double MosekBackend::SimplexIterations() const {
  return getIntAttr(MSK_IINF_SIM_PRIMAL_ITER);
}

int MosekBackend::BarrierIterations() const {
  return getIntAttr(MSK_IINF_INTPNT_ITER);
}

void MosekBackend::DoWriteProblem(const std::string &file) {
  MOSEK_CCALL(MSK_writedata(lp(), file.data()));
}


void MosekBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptMosek, lp());
  // TODO Check interrupter
}

MSKsoltypee MosekBackend::GetSolutionTypeToFetch() {
  int b = 0;

  // Handle MIP
  MOSEK_CCALL(MSK_solutiondef(lp(), MSK_SOL_ITG, &b));
  if (b) return MSK_SOL_ITG;

  // Handle LP
  MSKproblemtypee type;
  MOSEK_CCALL(MSK_getprobtype(lp(), &type));
  if (type == MSK_PROBTYPE_LO)
  {
    MOSEK_CCALL(MSK_solutiondef(lp(), MSK_SOL_BAS, &b));
    if (b) return MSK_SOL_BAS;
  }

  // All other case
  return MSK_SOL_ITR;
}

void MosekBackend::Solve() {
  MOSEK_CCALL(MSK_optimizetrm(lp(), &termCode_));
  solToFetch_ = GetSolutionTypeToFetch();
  MOSEK_CCALL(MSK_getsolsta(lp(), solToFetch_, &solSta_));
  WindupMOSEKSolve();
}

void MosekBackend::WindupMOSEKSolve() { }

void MosekBackend::ReportResults() {
  ReportMOSEKResults();
  BaseBackend::ReportResults();
}

void MosekBackend::ReportMOSEKResults() {
  SetStatus( ConvertMOSEKStatus() );
  AddMOSEKMessages();
  if (need_multiple_solutions())
    ReportMOSEKPool();
}

std::vector<double> MosekBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
  // TODO get solutions in the pool 
  return vars;
}

double MosekBackend::getPoolObjective(int i)
{
	double obj=0.0;
  // TODO get objective value of solution i
  return obj;
}

void MosekBackend::ReportMOSEKPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  // TODO
  /*
  while (++iPoolSolution < getIntAttr(MOSEK_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void MosekBackend::AddMOSEKMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> MosekBackend::ConvertMOSEKStatus() {
  namespace sol = mp::sol;
  std::string term_info = ConvertMOSEKTermStatus();
  switch (solSta_) {
  case MSK_SOL_STA_OPTIMAL:
  case MSK_SOL_STA_INTEGER_OPTIMAL:
    return { sol::SOLVED, "optimal" + term_info };
  case MSK_SOL_STA_PRIM_FEAS:
    return { sol::UNCERTAIN, "feasible primal" + term_info };
  case MSK_SOL_STA_DUAL_FEAS:
    return { sol::UNCERTAIN, "feasible dual" + term_info };
  case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
    return { sol::UNCERTAIN, "feasible solution" + term_info };
  case MSK_SOL_STA_PRIM_INFEAS_CER:
    return { sol::INFEASIBLE, "primal infeasible solution" + term_info };
  case MSK_SOL_STA_DUAL_INFEAS_CER:
    return { sol::INF_OR_UNB, "dual infeasible solution" + term_info };
  case MSK_SOL_STA_UNKNOWN:
  case MSK_SOL_STA_PRIM_ILLPOSED_CER:
  case MSK_SOL_STA_DUAL_ILLPOSED_CER:
  default:
    return { sol::UNKNOWN, "unknown" + term_info };
  }
}

std::string MosekBackend::ConvertMOSEKTermStatus() {
  switch (termCode_) {
  case MSK_RES_OK:
    break;
  case MSK_RES_TRM_MAX_ITERATIONS:
    return { ", max number of iterations reached" };
  case MSK_RES_TRM_MAX_TIME:
    return { ", maximum allowed time reached" };
  case MSK_RES_TRM_OBJECTIVE_RANGE:
    return { ", objective range" };
  case MSK_RES_TRM_STALL:
    return { ", stalling" };
  case MSK_RES_TRM_NUMERICAL_PROBLEM:
    return { ", numerical issues" };
  case MSK_RES_TRM_MIO_NUM_RELAXS:
  case MSK_RES_TRM_MIO_NUM_BRANCHES:
  case MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS:
    return { ", limit hit" };
  default:
    return { "termination code "
          + std::to_string(termCode_) };
  }
  return {};
}

void MosekBackend::FinishOptionParsing() {
  int v=-1;
  GetSolverOption(MSK_IPAR_LOG, v);
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo alg_values_method[] = {
  { "0", "Optimizer for conic constraints", 0},
  { "1", "Dual simplex", 1},
  { "2", "Automatic (default)", 2},
  { "3", "Free simplex", 3},
  { "4", "Interior-point method", 4},
  { "5", "Mixed-integer optimizer", 5},
  { "6", "Primal simplex", 6}
};

static const mp::OptionValueInfo alg_values_mip_presolve_use[] = {
  { "0", "Do not use presolve", 0},
  { "1", "Use presolve", 1},
  { "2", "Automatic (default)", 2}
};

void MosekBackend::InitCustomOptions() {

  set_option_header(
      "MOSEK Optimizer Options for AMPL\n"
      "--------------------------------\n\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``mosek_options``. For example::\n"
      "\n"
      "  ampl: option mosek_options 'threads=3';\n");

  AddSolverOption("alg:method method lpmethod simplex",
    "Which algorithm to use for non-MIP problems or for the root node of MIP problems:\n"
    "\n.. value-table::\n", MSK_IPAR_OPTIMIZER, alg_values_method, 2);

  AddSolverOption("lim:time timelim timelimit",
      "Limit on solve time (in seconds; default: no limit).",
      MSK_DPAR_OPTIMIZER_MAX_TIME, 0.0, DBL_MAX);

  AddStoredOption("mip:constructsol mipconstructsol",
      "Sets MSK_IPAR_MIO_CONSTRUCT_SOL. If set to MSK_ON and all integer variables "
      "have been given a value for which a feasible mixed integer solution exists, "
      "then MOSEK generates an initial solution to the mixed integer problem by "
      "fixing all integer values and solving the remaining problem."
      "Default = OFF",
      storedOptions_.MIPConstructSol_);

  AddSolverOption("mip:presolve presolve",
    "MIP presolve:\n"
                  "\n.. value-table::\n",
    MSK_IPAR_PRESOLVE_USE, alg_values_mip_presolve_use, 2);

  AddSolverOption("mip:inttol inttol",
    "MIP integrality tolerance.",
    MSK_DPAR_MIO_TOL_ABS_RELAX_INT, 1e-15, Infinity());

  AddSolverOption("tech:threads threads",
    "Controls the number of threads employed by the optimizer. "
    "Default 0 ==> number of threads used will be equal to the number "
    "of cores detected on the machine.",
    MSK_IPAR_NUM_THREADS, 0, INT_MAX);

  AddSolverOption("mip:relgapconst miorelgapconst",
    "This value is used to compute the relative gap for the solution "
    "to an integer optimization problem."
    "Default = 1.0e-10",
    MSK_DPAR_MIO_REL_GAP_CONST, 0.0, DBL_MAX);

  AddSolverOption("tech:outlev outlev",
    "0*/1: Whether to write mosek log lines to stdout.",
    MSK_IPAR_LOG, 0, 1);
}

double MosekBackend::MIPGap() {
  return getDblAttr(MSK_DINF_MIO_OBJ_REL_GAP);
}

double MosekBackend::BestDualBound() {
  return getDblAttr(MSK_DINF_MIO_OBJ_BOUND);
}

double MosekBackend::MIPGapAbs() {
  return getDblAttr(MSK_DINF_MIO_OBJ_ABS_GAP);
}

SensRangesPresolved MosekBackend::GetSensRangesPresolved()
{
  // Set sensitivity range parameters for constraints
  int lencon = NumLinCons();
  std::vector<MSKint32t> cindex(lencon);
  std::vector<MSKmarke> cmarklb(lencon);
  std::vector<MSKmarke> cmarkub(lencon);
  for (int i=0; i<lencon; i++)
  {
    cindex[i] = i;
    cmarklb[i] = MSK_MARK_LO;
    cmarkub[i] = MSK_MARK_UP;
  }
  std::vector<MSKrealt> cpricelblo(lencon);
  std::vector<MSKrealt> cpricelbhi(lencon);
  std::vector<MSKrealt> cpriceublo(lencon);
  std::vector<MSKrealt> cpriceubhi(lencon);
  std::vector<MSKrealt> crangelblo(lencon);
  std::vector<MSKrealt> crangelbhi(lencon);
  std::vector<MSKrealt> crangeublo(lencon);
  std::vector<MSKrealt> crangeubhi(lencon);

  // Set sensitivity range parameters for variables
  int lenvar = NumVars();
  std::vector<MSKint32t> vindex(lenvar);
  std::vector<MSKmarke> vmarklb(lenvar);
  std::vector<MSKmarke> vmarkub(lenvar);
  for (int i=0; i<lenvar; i++)
  {
    vindex[i] = i;
    vmarklb[i] = MSK_MARK_LO;
    vmarkub[i] = MSK_MARK_UP;
  }
  std::vector<MSKrealt> vpricelblo(lenvar);
  std::vector<MSKrealt> vpricelbhi(lenvar);
  std::vector<MSKrealt> vpriceublo(lenvar);
  std::vector<MSKrealt> vpriceubhi(lenvar);
  std::vector<MSKrealt> vrangelblo(lenvar);
  std::vector<MSKrealt> vrangelbhi(lenvar);
  std::vector<MSKrealt> vrangeublo(lenvar);
  std::vector<MSKrealt> vrangeubhi(lenvar);

  std::vector<MSKrealt> opricelo(lenvar);
  std::vector<MSKrealt> opricehi(lenvar);
  std::vector<MSKrealt> orangelo(lenvar);
  std::vector<MSKrealt> orangehi(lenvar);

  // Compute sensitivity ranges for lower bounds
  MOSEK_CCALL(
    MSK_primalsensitivity(
        lp(),
        lencon,
        cindex.data(),
        cmarklb.data(),
        lenvar,
        vindex.data(),
        vmarklb.data(),
        cpricelblo.data(),
        cpricelbhi.data(),
        crangelblo.data(),
        crangelbhi.data(),
        vpricelblo.data(),
        vpricelbhi.data(),
        vrangelblo.data(),
        vrangelbhi.data()
    )
  );

  // Compute sensitivity ranges for upper bounds
  MOSEK_CCALL(
    MSK_primalsensitivity(
        lp(),
        lencon,
        cindex.data(),
        cmarkub.data(),
        lenvar,
        vindex.data(),
        vmarkub.data(),
        cpriceublo.data(),
        cpriceubhi.data(),
        crangeublo.data(),
        crangeubhi.data(),
        vpriceublo.data(),
        vpriceubhi.data(),
        vrangeublo.data(),
        vrangeubhi.data()
    )
  );

  // Compute sensitivity ranges for objective coeffs
  MOSEK_CCALL(
    MSK_dualsensitivity(
        lp(),
        lenvar,
        vindex.data(),
        opricelo.data(),
        opricehi.data(),
        orangelo.data(),
        orangehi.data()
    )
  );

  // Get bounds and add them to sensitivity ranges
  std::vector<MSKboundkeye> vbk(lenvar);
  std::vector<MSKrealt> vbl(lenvar), vbu(lenvar);
  MOSEK_CCALL(MSK_getvarboundslice(lp(), (MSKint32t)0, (MSKint32t)lenvar, vbk.data(), vbl.data(), vbu.data()));
  std::transform(vrangelblo.begin(), vrangelblo.end(), vbl.begin(), vrangelblo.begin(), std::plus<MSKrealt>());
  std::transform(vrangelbhi.begin(), vrangelbhi.end(), vbl.begin(), vrangelbhi.begin(), std::plus<MSKrealt>());
  std::transform(vrangeublo.begin(), vrangeublo.end(), vbu.begin(), vrangeublo.begin(), std::plus<MSKrealt>());
  std::transform(vrangeubhi.begin(), vrangeubhi.end(), vbu.begin(), vrangeubhi.begin(), std::plus<MSKrealt>());
  for (int i=0; i<lenvar; i++)
  {
    if (vbk[i] == MSK_BK_UP || vbk[i] == MSK_BK_FR)
    {
      vrangelblo[i] = -MSK_INFINITY;
      vrangelbhi[i] = MSK_INFINITY;
    }
    if (vbk[i] == MSK_BK_LO || vbk[i] == MSK_BK_FR)
    {
      vrangeublo[i] = -MSK_INFINITY;
      vrangeubhi[i] = MSK_INFINITY;
    }
  }

  std::vector<MSKboundkeye> cbk(lencon);
  std::vector<MSKrealt> cbl(lencon), cbu(lencon);
  MOSEK_CCALL(MSK_getconboundslice(lp(), (MSKint32t)0, (MSKint32t)lencon, cbk.data(), cbl.data(), cbu.data()));
  std::transform(crangelblo.begin(), crangelblo.end(), cbl.begin(), crangelblo.begin(), std::plus<MSKrealt>());
  std::transform(crangelbhi.begin(), crangelbhi.end(), cbl.begin(), crangelbhi.begin(), std::plus<MSKrealt>());
  std::transform(crangeublo.begin(), crangeublo.end(), cbu.begin(), crangeublo.begin(), std::plus<MSKrealt>());
  std::transform(crangeubhi.begin(), crangeubhi.end(), cbu.begin(), crangeubhi.begin(), std::plus<MSKrealt>());
  for (int i=0; i<lencon; i++)
  {
    if (cbk[i] == MSK_BK_UP || cbk[i] == MSK_BK_FR)
    {
      crangelblo[i] = -MSK_INFINITY;
      crangelbhi[i] = MSK_INFINITY;
    }
    if (cbk[i] == MSK_BK_LO || cbk[i] == MSK_BK_FR)
    {
      crangeublo[i] = -MSK_INFINITY;
      crangeubhi[i] = MSK_INFINITY;
    }
  }

  std::vector<MSKrealt> c(lenvar);
  MOSEK_CCALL(MSK_getc(lp(), c.data()));
  std::transform(orangelo.begin(), orangelo.end(), c.begin(), orangelo.begin(), std::plus<MSKrealt>());
  std::transform(orangehi.begin(), orangehi.end(), c.begin(), orangehi.begin(), std::plus<MSKrealt>());

  /// Return sensitivity ranges
  SensRangesPresolved sensr;

  sensr.varlbhi = { {vrangelbhi} };
  sensr.varlblo = { {vrangelblo} };
  sensr.varubhi = { {vrangeubhi} };
  sensr.varublo = { {vrangeublo} };
  sensr.varobjhi = { {orangehi} };
  sensr.varobjlo = { {orangelo} };
	sensr.conlbhi = { {}, {{{CG_Algebraic, crangelbhi}}} };
	sensr.conlblo = { {}, {{{CG_Algebraic, crangelblo}}} };
	sensr.conubhi = { {}, {{{CG_Algebraic, crangeubhi}}} };
	sensr.conublo = { {}, {{{CG_Algebraic, crangeublo}}} };
  std::vector<MSKrealt> rhs(lencon, 0.0);
	sensr.conrhshi = { {}, {{{CG_Algebraic, rhs}}} };
	sensr.conrhslo = { {}, {{{CG_Algebraic, rhs}}} };

  return sensr;
}

ArrayRef<int> MosekBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  MOSEK_CCALL(MSK_getskx(lp(), solToFetch_, (MSKstakeye *)vars.data()));

  for (auto& s : vars)
  {
    switch (s)
    {
      case MSK_SK_UNK:
        s = (int)BasicStatus::none;
        break;
      case MSK_SK_BAS:
        s = (int)BasicStatus::bas;
        break;
      case MSK_SK_LOW:
        s = (int)BasicStatus::low;
        break;
      case MSK_SK_UPR:
        s = (int)BasicStatus::upp;
        break;
      case MSK_SK_SUPBAS:
        s = (int)BasicStatus::sup;
        break;
      case MSK_SK_FIX:
        s = (int)BasicStatus::equ;
        break;
      default:
        MP_RAISE(fmt::format("Unknown Mosek VBasis value: {}", s));
    }
  }

  return vars;
}

ArrayRef<int> MosekBackend::ConStatii() {
  std::vector<int> cons(NumLinCons());
  MOSEK_CCALL(MSK_getskc(lp(), solToFetch_, (MSKstakeye *)cons.data()));

  for (auto& s : cons)
  {
    switch (s)
    {
      case MSK_SK_UNK:
        s = (int)BasicStatus::none;
        break;
      case MSK_SK_BAS:
        s = (int)BasicStatus::bas;
        break;
      case MSK_SK_LOW:
        s = (int)BasicStatus::low;
        break;
      case MSK_SK_UPR:
        s = (int)BasicStatus::upp;
        break;
      case MSK_SK_SUPBAS:
        s = (int)BasicStatus::sup;
        break;
      case MSK_SK_FIX:
        s = (int)BasicStatus::equ;
        break;
      default:
        MP_RAISE(fmt::format("Unknown Mosek VBasis value: {}", s));
    }
  }
  return cons;
}

void MosekBackend::VarStatii(ArrayRef<int> vst) {
  // BasicStatus enum values: from 0 to 6: none, bas, sup, low, upp, equ, btw
  std::vector<int> skx(vst.data(), vst.data() + vst.size());

  // TODO: remove this if we are sure that it is the same.
  MSKint32t numvar;
  MOSEK_CCALL(MSK_getnumvar(lp(), &numvar));
  assert(numvar == vst.size());

  for (auto j = 0; j < skx.size(); j++)
  {
    switch ((BasicStatus)skx[j])
    {
      case BasicStatus::bas:
        skx[j] = MSK_SK_BAS;
        break;
      case BasicStatus::low:
        skx[j] = MSK_SK_LOW;
        break;
      case BasicStatus::equ:
        skx[j] = MSK_SK_FIX;
        break;
      case BasicStatus::upp:
        skx[j] = MSK_SK_UPR;
        break;
      case BasicStatus::sup:
      case BasicStatus::btw:
        skx[j] = MSK_SK_SUPBAS;
        break;
      case BasicStatus::none:
        /// 'none' is assigned to new variables (added after solving).
        /// Compute low/upp/fix/sup, depending on where initial value 0.0 is between the bounds
        MSKboundkeye bk;
        MSKrealt bl, bu, tolx;
        MOSEK_CCALL(MSK_getvarbound(lp(), (MSKint32t)j, &bk, &bl, &bu));
        MOSEK_CCALL(MSK_getdouparam(lp(), MSK_DPAR_DATA_TOL_X, &tolx));

        if (0.0 < bl + tolx) skx[j] = MSK_SK_LOW;
        else if (bu - tolx < 0.0) skx[j] = MSK_SK_UPR;
        else if (bu - bl < tolx) skx[j] = MSK_SK_FIX;
        else skx[j] = MSK_SK_SUPBAS;

        break;
      default:
        MP_RAISE(fmt::format("Unknown AMPL var status value: {}", skx[j]));
    }
  }

  // Input the variable status keys
  MOSEK_CCALL(MSK_putskx(lp(), solToFetch_, (MSKstakeye *)skx.data()));
}

void MosekBackend::ConStatii(ArrayRef<int> cst) {
  MSKint32t numcon;
  MOSEK_CCALL(MSK_getnumcon(lp(), &numcon));
	std::vector<int> skc(cst.data(), cst.data() + numcon);

  // The status key of a constraint is the status key of the logical (slack) variable assigned to it.
  // The slack is defined as: l <= a'x <= u rewritten as a'x - s = 0, l <= s <= u.
	for (size_t j = 0; j < skc.size(); j++)
  {
    skc[j] = MSK_SK_UNK;
    switch ((BasicStatus)skc[j])
    {
      case BasicStatus::bas:
        skc[j] = MSK_SK_BAS;
        break;
      case BasicStatus::sup:
      case BasicStatus::btw:
        skc[j] = MSK_SK_SUPBAS;
        break;
      case BasicStatus::low:
        skc[j] = MSK_SK_LOW;
        break;
      case BasicStatus::upp:
        skc[j] = MSK_SK_UPR;
        break;
      case BasicStatus::equ:
        skc[j] = MSK_SK_FIX;
        break;
      case BasicStatus::none:
        /// 'none' is assigned to new constraints (added after solving).
        /// Compute low/upp/fix/sup, depending on where initial value 0.0 is between the bounds
        MSKboundkeye bk;
        MSKrealt bl, bu, tolx;
        MOSEK_CCALL(MSK_getconbound(lp(), (MSKint32t)j, &bk, &bl, &bu));
        // MSK_DPAR_DATA_TOL_X is both for variables and constraints
        MOSEK_CCALL(MSK_getdouparam(lp(), MSK_DPAR_DATA_TOL_X, &tolx));

        if (0.0 < bl + tolx) skc[j] = MSK_SK_LOW;
        else if (bu - tolx < 0.0) skc[j] = MSK_SK_UPR;
        else if (bu - bl < tolx) skc[j] = MSK_SK_FIX;
        else skc[j] = MSK_SK_SUPBAS;

        break;
      default:
        MP_RAISE(fmt::format("Unknown AMPL var status value: {}", skc[j]));
    }
  }

  // Input the variable status keys
  MOSEK_CCALL(MSK_putskc(lp(), solToFetch_, (MSKstakeye *)skc.data()));
}

SolutionBasis MosekBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
				{{{ CG_Algebraic, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
    assert(constt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void MosekBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
	auto constt = mv.GetConValues()(CG_Algebraic);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

void MosekBackend::AddPrimalDualStart(Solution sol)
{
  solToFetch_ = GetSolutionTypeToFetch();
  auto mv = GetValuePresolver().PresolveSolution(
        { sol.primal, sol.dual } );
  auto x0 = mv.GetVarValues()();
	auto pi0 = mv.GetConValues()(CG_Algebraic);
  MOSEK_CCALL(MSK_putxx(lp(), solToFetch_, (MSKrealt *)x0.data()));
  MOSEK_CCALL(MSK_puty(lp(), solToFetch_, (MSKrealt *)pi0.data()));
}

void MosekBackend::AddMIPStart(ArrayRef<double> x0_unpres)
{
  solToFetch_ = GetSolutionTypeToFetch();
  auto mv = GetValuePresolver().PresolveSolution( { x0_unpres } );
  auto x0 = mv.GetVarValues()();

  if (Mosek_mip_construct_sol())
  {
    // Use integer part of given solution to construct initial solution by solving for continuous part.
    MOSEK_CCALL(MSK_putintparam(lp(), MSK_IPAR_MIO_CONSTRUCT_SOL, MSK_ON));
  }
  MOSEK_CCALL(MSK_putxx(lp(), solToFetch_, (MSKrealt *)x0.data()));
}

ArrayRef<double> MosekBackend::Ray()
{
  // Problem is checked to be dual infeasible at this point, so the ray is the primal solution variable.
  // (Primal can be unbounded or infeasible.)
  std::vector<double> xx(NumVars());
  MOSEK_CCALL(MSK_getxx(lp(), solToFetch_, (MSKrealt *)xx.data()));
  // Argument is a ModelValues<ValueMap>, which can be constructed from a ValueMap for only the variables.
  auto mv = GetValuePresolver().PostsolveSolution({xx});
  // We get the ValueMap for variables, and by calling that, the value itself.
  // (Similar to .MoveOut(0) for single key maps)
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> MosekBackend::DRay()
{
  // Problem is checked to be (primal) infeasible at this point, so the ray is the dual solution variable.
  // (Dual can be unbounded or infeasible.)
  std::vector<double> y(NumLinCons());
  MOSEK_CCALL(MSK_gety(lp(), solToFetch_, (MSKrealt *)y.data()));
  // Argument is a ModelValues<ValueMap>, which is constructed now from two ValueMaps, an empty for variables, and
  // a second one for constraints. The constraints ValueMap is given as a std::map for only linear constraints
	auto mv = GetValuePresolver().
			PostsolveSolution({
													{},
													{ { {CG_Algebraic, std::move(y)} } }
												});
  // We get the ValueMap for constraints, and by calling MoveOut(), the value itself.
  // (Similar to .MoveOut(0) for single key maps)
  return mv.GetConValues().MoveOut();
}

} // namespace mp


// AMPLs

AMPLS_MP_Solver* AMPLSOpenMosek(const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::MosekBackend()},
        slv_opt, cb);
}

void AMPLSCloseMosek(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

MSKtask_t GetMosekmodel(AMPLS_MP_Solver* slv) {
  return dynamic_cast<mp::MosekBackend*>(AMPLSGetBackend(slv))->lp();
}

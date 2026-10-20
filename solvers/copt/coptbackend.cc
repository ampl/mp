#include <vector>
#include <algorithm>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "coptbackend.h"

extern "C" {
  #include "copt-ampls-c-api.h"    // Copt AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptCopt(void* prob) {
  return COPT_Interrupt((copt_prob*)prob);
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCoptBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CoptBackend()};
}


namespace mp {

/// Create Copt Model Manager
/// @param gc: the Copt common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return CoptModelMgr
std::unique_ptr<BasicModelManager>
CreateCoptModelMgr(CoptCommon&, Env&, pre::BasicValuePresolver*&);


void CoptBackend::OpenSolver() {
  int status = 0;
  const auto create_fn = GetCallbacks().init;
  if (create_fn)
    set_env((copt_env*)create_fn());
  else
    COPT_CCALL(COPT_CreateEnv(&env_ref()));
  if (env() == NULL) {
    throw std::runtime_error(
      fmt::format("Could not open COPT environment.\n{}", status));
  }
  
  /* Create an empty model */
  status = COPT_CreateProb(env(), &lp_ref());
  if (status)
    throw std::runtime_error(fmt::format(
      "Failed to create problem, error code {}.", status));
  COPT_CCALL(COPT_SetIntParam(lp(), "Logging", 0));
  /* Copy handlers to ModelAPI */
  copy_common_info_to_other();
  SetSolverOption(COPT_INTPARAM_LOGGING, 0);
}

void CoptBackend::CloseSolver() {
  if (lp() != NULL) {
    COPT_CCALL(COPT_DeleteProb(&lp_ref()));
  }
  if (env() != NULL) {
    COPT_CCALL(COPT_DeleteEnv(&env_ref()));
  }
}

CoptBackend::CoptBackend() {
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateCoptModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

CoptBackend::~CoptBackend() {
  CloseSolver();
}

const char* CoptBackend::GetBackendName()
  { return "CoptBackend"; }

std::string CoptBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", COPT_VERSION_MAJOR, 
    COPT_VERSION_MINOR, COPT_VERSION_TECHNICAL);
}


bool CoptBackend::IsMIP() const {
  // return getIntAttr(COPT_INTATTR_ISMIP);  -- fails in 8.0.3 with nonconvex=2
    return COPT_MIPSTATUS_UNSTARTED != getIntAttr(COPT_INTATTR_MIPSTATUS);
  // @todo Is this needed before solving?
}

bool CoptBackend::IsQCP() const {
  return getIntAttr(COPT_INTATTR_QCONSTRS) > 0;
}
void CoptBackend::InitOptionParsing() {
  OpenSolver();
}


void CoptBackend::InputExtras() {
  MIPBackend::InputExtras();
  InputCOPTExtras();
}

void CoptBackend::InputCOPTExtras() {
  if (need_ray_primal() || need_ray_dual())
    SetSolverOption(COPT_INTPARAM_REQFARKASRAY, 1);
  if (sensitivity())
    SetSolverOption(COPT_INTPARAM_REQSENSITIVITY, 1);
}

void CoptBackend::MarkLazyOrUserCuts(ArrayRef<int> lazyVals) {
  auto vm = GetValuePresolver().PresolveLazyUserCutFlags({ {}, lazyVals });
  std::vector<int> lv_lin = vm.GetConValues()(CG_Linear);
  int m = NumLinCons();
  if ((int)lv_lin.size() > m)
    lv_lin.resize(m);
  auto rowlb = getConInfo(COPT_DBLINFO_LB);
  auto rowub = getConInfo(COPT_DBLINFO_UB);
  int n_lazy = 0, n_user = 0;
  std::vector<int> idx;
  std::vector<double> elem;
  for (int i=0; i<(int)lv_lin.size(); ++i) {
    if (lv_lin[i] > 0 && lazy_cuts()) {
      int req = 0, beg = 0, cnt = 0;
      COPT_CCALL(COPT_GetRows(lp(), 1, &i, NULL, NULL, NULL, NULL, 0, &req));
      idx.resize(req);
      elem.resize(req);
      COPT_CCALL(COPT_GetRows(lp(), 1, &i, &beg, &cnt,
                              idx.data(), elem.data(), req, &req));
      char name[1024] = "";
      int namelen = 0;
      COPT_GetRowName(lp(), i, name, (int)sizeof(name), &namelen);
      char sense = COPT_FREE;
      double bound = 0.0, upper = 0.0;
      bool lbfin = rowlb[i] > -COPT_INFINITY, ubfin = rowub[i] < COPT_INFINITY;
      if (lbfin && ubfin && rowlb[i] == rowub[i]) {
        sense = COPT_EQUAL; bound = rowlb[i];
      } else if (lbfin && ubfin) {
        sense = COPT_RANGE; bound = rowlb[i]; upper = rowub[i];
      } else if (ubfin) {
        sense = COPT_LESS_EQUAL; bound = rowub[i];
      } else if (lbfin) {
        sense = COPT_GREATER_EQUAL; bound = rowlb[i];
      }
      COPT_CCALL(COPT_AddLazyConstr(lp(), cnt, idx.data(), elem.data(),
                                    sense, bound, upper,
                                    namelen > 0 ? name : NULL));
      double ninf = -COPT_INFINITY, pinf = COPT_INFINITY;    // free the original row
      COPT_CCALL(COPT_SetRowLower(lp(), 1, &i, &ninf));
      COPT_CCALL(COPT_SetRowUpper(lp(), 1, &i, &pinf));
      lazy_rows_.push_back({i, rowlb[i], rowub[i]});
      ++n_lazy;
    } else if (lv_lin[i] < 0 && user_cuts())
      ++n_user;
  }
  if (n_user)
    AddWarning("USER_CUTS",
               "COPT has no user cuts; constraints with .lazy < 0 "
               "stay ordinary constraints");
  // Testing
  if (debug_mode())
    ReportProblemSuffix("test_lin_constr_lazy", n_lazy);
}

ArrayRef<double> CoptBackend::GetObjectiveValues() {
  if (HasFixedLP()) return std::vector<double>{FixedLP().objective};
  if (feasrelax_solved_)
    return std::vector<double>{ ObjectiveValue() };
  int nobjs = NumObjs();
  if (nobjs > 1) {
    std::vector<double> v(nobjs);
    for (int i = 0; i < v.size(); ++i)
    {
      if (IsMIP())
        COPT_MultiObjGetDblAttr(lp(), i, COPT_DBLATTR_BESTOBJ, &v[i]);
      else
        COPT_MultiObjGetDblAttr(lp(), i, COPT_DBLATTR_LPOBJVAL, &v[i]);
    }
    return v;
  } else {
    std::vector<double> v(1, 0.0);
    if (IsMIP())
      COPT_GetDblAttr(lp(), COPT_DBLATTR_BESTOBJ, &v[0]);
    else
      COPT_GetDblAttr(lp(), COPT_DBLATTR_LPOBJVAL, &v[0]);
    return v;
  }
}

void CoptBackend::ObjPriorities(ArrayRef<int> values) {
  // SetSolverOption(COPT_INTPARAM_MULTIOBJPARAMMODE, 1);
  for (int i = 0; i < values.size(); ++i) {
    COPT_CCALL(COPT_MultiObjSetObjParam(
        lp(), i, COPT_MULTIOBJ_PRIORITY, static_cast<double>(values[i])));
  }
}

void CoptBackend::ObjWeights(ArrayRef<double> values) {
  // SetSolverOption(COPT_INTPARAM_MULTIOBJPARAMMODE, 1);
  for (int i = 0; i < values.size(); ++i)
    COPT_CCALL(COPT_MultiObjSetObjParam(
        lp(), i, COPT_MULTIOBJ_WEIGHT, static_cast<double>(values[i])));
}

void CoptBackend::ObjAbsTol(ArrayRef<double> values) {
  // SetSolverOption(COPT_INTPARAM_MULTIOBJPARAMMODE, 1);
  for (int i = 0; i < values.size(); ++i) {
    COPT_CCALL(COPT_MultiObjSetObjParam(
        lp(), i, COPT_MULTIOBJ_ABSTOL, static_cast<double>(values[i])));
  }
}

void CoptBackend::ObjRelTol(ArrayRef<double> values) {
  // SetSolverOption(COPT_INTPARAM_MULTIOBJPARAMMODE, 1);
  for (int i = 0; i < values.size(); ++i) {
    COPT_CCALL(COPT_MultiObjSetObjParam(
        lp(), i, COPT_MULTIOBJ_RELTOL, static_cast<double>(values[i])));
  }
}

void CoptBackend::SetMultiobjOptions(BasicObjOptionSetter* setter) {
    auto passes = setter->GetPassesWithOptions();
    if (!passes.empty())   // COPT applies per-objective parameters only in this mode
      SetSolverOption(COPT_INTPARAM_MULTIOBJPARAMMODE, 1);
    for (auto pass : passes) {
        set_current_objective_options(pass);
        setter->SetOptionsForMultiobjPass(pass);
    }
}

bool CoptBackend::IsMIPModel() const {
  return getIntAttr(COPT_INTATTR_INTS) + getIntAttr(COPT_INTATTR_BINS) > 0;
}

void CoptBackend::DoCOPTFeasRelax() {
  /// AMPL modes 1..6 -> COPT FeasRelaxMode:
  /// 0/1 sum, 2/3 count, 4/5 squares; odd = optimize original objective
  int map[] = { 0, 4, 2, 1, 5, 3 };
  feasrelax_copt_mode_ = map[feasrelax()-1];
  if (feasrelax_copt_mode_ >= 4 && IsMIPModel())
    MP_RAISE("COPT does not support the squared-violation feasibility "
             "relaxation (alg:feasrelax=2 or 5) for models with integer "
             "variables; use alg:feasrelax=1, 3, 4 or 6.");
  if (feasrelax_copt_mode_ % 2)
    feasrelax().flag_orig_obj_available();
  SetSolverOption(COPT_INTPARAM_FEASRELAXMODE, feasrelax_copt_mode_);
  /// only relax linear constraints
  auto mv = GetValuePresolver().PresolveSolution({
                                              {},
                                              feasrelax().rhspen()
    });
  frx_rhspen_ = mv.GetConValues()(CG_Linear);
  /// Account for new variables (TODO: proper presolve)
  frx_lbpen_ = feasrelax().lbpen();
  if (frx_lbpen_.size() && frx_lbpen_.size() < (size_t)NumVars())
    frx_lbpen_.resize(NumVars());
  frx_ubpen_ = feasrelax().ubpen();
  if (frx_ubpen_.size() && frx_ubpen_.size() < (size_t)NumVars())
    frx_ubpen_.resize(NumVars());

  /// NULL penalties mean "no relaxation of these bounds",
  /// COPT_INFINITY entries block individual bounds, as in MP.
  COPT_CCALL(COPT_FeasRelax(lp(),
    (double*)data_or_null(frx_lbpen_),
    (double*)data_or_null(frx_ubpen_),
    (double*)data_or_null(frx_rhspen_),
    NULL));
  feasrelax_solved_ = true;
  if (feasrelax_copt_mode_ % 2
      && getIntAttr(COPT_INTATTR_HASFEASRELAXSOL))
    feasrelax().orig_obj_value() = FeasRelaxViolation();
}

double CoptBackend::FeasRelaxViolation() const {
  auto self = const_cast<CoptBackend*>(this);
  auto clb = self->getVarInfo(COPT_DBLINFO_RELAXLB);
  auto cub = self->getVarInfo(COPT_DBLINFO_RELAXUB);
  auto rlb = self->getConInfo(COPT_DBLINFO_RELAXLB);
  auto rub = self->getConInfo(COPT_DBLINFO_RELAXUB);
  int kind = feasrelax_copt_mode_ / 2;         // 0 sum, 1 count, 2 squares
  double total = 0.0;
  auto acc = [&](const std::vector<double>& v, const std::vector<double>& pen) {
    for (size_t i=0; i<v.size(); ++i) {
      double d = std::fabs(v[i]);
      if (d <= 1e-9)
        continue;
      double w = i < pen.size() ? pen[i] : 1.0;
      if (w >= COPT_INFINITY)
        continue;
      total += (0==kind) ? w*d : (1==kind) ? 1.0 : w*d*d;
    }
  };
  acc(clb, frx_lbpen_);
  acc(cub, frx_ubpen_);
  acc(rlb, frx_rhspen_);
  acc(rub, frx_rhspen_);
  return total;
}

ArrayRef<double> CoptBackend::PrimalSolution() {
  if (HasFixedLP()) return FixedLP().primal;
  if (feasrelax_solved_) {
    if (!getIntAttr(COPT_INTATTR_HASFEASRELAXSOL))
      return {};
    return getVarInfo(COPT_DBLINFO_RELAXVALUE);
  }
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  if (IsMIP()) 
    error = COPT_GetSolution(lp(), x.data());
  else
    error = COPT_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);

  if (error)
    x.clear();
  return x;
}

pre::ValueMapDbl CoptBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> CoptBackend::DualSolution_LP() {
  if (IsMIP()) return HasFixedLP() ? ArrayRef<double>(FixedLP().dual) : ArrayRef<double>{};
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  if (!IsMIP()) {
    int error = COPT_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
    if (error)
      pi.clear();
  } else
    pi.clear();
  return pi;
}

double CoptBackend::ObjectiveValue() const {
  if (feasrelax_solved_)
    return getDblAttr(COPT_DBLATTR_FEASRELAXOBJ);
  if (IsMIP())
    return getDblAttr(COPT_DBLATTR_BESTOBJ);
  else
    return getDblAttr(COPT_DBLATTR_LPOBJVAL);
}



int CoptBackend::NodeCount() const {
  return getIntAttr(COPT_INTATTR_NODECNT);
}

int CoptBackend::SimplexIterations() const {
  return getIntAttr(COPT_INTATTR_SIMPLEXITER);
}

int CoptBackend::BarrierIterations() const {
  return getIntAttr(COPT_INTATTR_BARRIERITER);
}

std::map<std::string, std::variant<int, double, std::string>>
CoptBackend::SolutionStats() {
    std::map<std::string, std::variant<int, double, std::string>> stats;
    stats["simplex_iterations"] = SimplexIterations();
    stats["barrier_iterations"] = BarrierIterations();
    stats["node_count"] = NodeCount();
    return stats;
}
void CoptBackend::DoWriteProblem(const std::string& name) {
  auto CheckEnds = [&](std::string e) {
    return ends_with(name, e) || ends_with(name, e + ".gz");
  };
  if (CheckEnds(".lp"))
    COPT_CCALL(COPT_WriteLp(lp(), name.c_str()));
  else if (CheckEnds(".mps"))
    COPT_CCALL(COPT_WriteMps(lp(), name.c_str()));
  else if (CheckEnds(".nl"))
    COPT_CCALL(COPT_WriteNL(lp(), name.c_str()));
  else if (CheckEnds(".cbf"))
    COPT_CCALL(COPT_WriteCbf(lp(), name.c_str()));
  else if (CheckEnds(".bin"))
    COPT_CCALL(COPT_WriteBin(lp(), name.c_str()));
  else if (CheckEnds(".par"))
    COPT_CCALL(COPT_WriteParam(lp(), name.c_str()));
  else
    MP_RAISE(
        fmt::format("Unknown export format: '{}'.", name));
}

void CoptBackend::DoWriteSolution(const std::string& name) {
  auto CheckEnds = [&](std::string e) {
    return ends_with(name, e) || ends_with(name, e + ".gz");
  };
  if (CheckEnds(".bas"))
    COPT_CCALL(COPT_WriteBasis(lp(), name.c_str()));
  else if (CheckEnds(".iis"))
    COPT_CCALL(COPT_WriteIIS(lp(), name.c_str()));
  else if (CheckEnds(".sol"))
    COPT_CCALL(COPT_WriteSol(lp(), name.c_str()));
  else
    MP_RAISE(
        fmt::format("Unknown result export format: '{}'.", name));
}


void CoptBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptCopt, lp());
}

void CoptBackend::Solve() {
  if (!storedOptions_.tunebase_.empty())
    DoCOPTTune();
  if (feasrelax())
    DoCOPTFeasRelax();
  else
    COPT_CCALL(COPT_Solve(lp()));
  WindupCOPTSolve();
}

void CoptBackend::DoCOPTTune() {
  COPT_CCALL(COPT_Tune(lp()));
  int n_results = getIntAttr(COPT_INTATTR_TUNERESULTS);
  if (n_results <= 0) {
    AddWarning("TUNING", "No tuning results: COPT found no parameter "
               "settings improving on the defaults.");
    return;
  }
  auto tbc = storedOptions_.tunebase_;
  if (tbc.size()>=4 && 0==tbc.compare(tbc.size()-4, 4, ".par"))
    tbc.resize(tbc.size()-4);
  tbc += "_{}_";
  tbc += ".par";
  std::string tfn;
  for (int k=n_results; k--; ) {   // result 0 is the best one
    tfn = fmt::format(tbc.c_str(), k+1);
    COPT_CCALL(COPT_WriteTuneParam(lp(), k, tfn.c_str()));
  }
  COPT_CCALL(COPT_LoadTuneParam(lp(), 0));      // solve with the best set
  AddToSolverMessage(
        fmt::format("Tuning: wrote {} parameter files, best file: '{}'\n",
                    n_results, tfn));
}

const char* CoptBackend::FixedModelIterationLimitHelp() const {
  return "Iteration limit for the fixed-model solve only. -1 (default) uses simplex "
         "with the solver default limit. Setting this option selects barrier with "
         "crossover and limits barrier iterations. 0 skips the supplementary solve.";
}

CoptBackend::FixedLPResult CoptBackend::SolveFixedLP(ArrayRef<double> x) {
  if (!getIntAttr(COPT_INTATTR_HASMIPSOL)) MP_RAISE("no integer-feasible incumbent");
  if (NumObjs() > 1) MP_RAISE("multiple objectives are unsupported by fixed-model emulation");
  copt_prob* copy = nullptr;
  COPT_CCALL(COPT_CreateCopy(lp(), &copy));
  auto destroy = [](copt_prob* model) { COPT_DeleteProb(&model); };
  std::unique_ptr<copt_prob, decltype(destroy)> fixed(copy, destroy);
  int n = NumVars(), m = NumLinCons();
  {                      // the copy's log follows tech:outlev
    int logging = 0;
    GetSolverOption(COPT_INTPARAM_LOGGING, logging);
    COPT_CCALL(COPT_SetIntParam(copy, COPT_INTPARAM_LOGGING, logging));
  }
  for (const auto& row: lazy_rows_) {   // lazy constraints do not apply to the LP
    COPT_CCALL(COPT_SetRowLower(copy, 1, &row.index, &row.lb));
    COPT_CCALL(COPT_SetRowUpper(copy, 1, &row.index, &row.ub));
  }
  std::vector<char> types(n);
  std::vector<double> lb(n), ub(n);
  COPT_CCALL(COPT_GetColType(lp(), n, nullptr, types.data()));
  COPT_CCALL(COPT_GetColInfo(lp(), COPT_DBLINFO_LB, n, nullptr, lb.data()));
  COPT_CCALL(COPT_GetColInfo(lp(), COPT_DBLINFO_UB, n, nullptr, ub.data()));
  double tol;
  GetSolverOption(COPT_DBLPARAM_INTTOL, tol);
  for (int j = 0; j < n; ++j) {
    if (types[j] == COPT_BINARY || types[j] == COPT_INTEGER) {
      double value = FixedIntegerValue(x[j], lb[j], ub[j], tol);
      char continuous = COPT_CONTINUOUS;
      COPT_CCALL(COPT_SetColType(copy, 1, &j, &continuous));
      COPT_CCALL(COPT_SetColLower(copy, 1, &j, &value));
      COPT_CCALL(COPT_SetColUpper(copy, 1, &j, &value));
    } else if (types[j] != COPT_CONTINUOUS)
      MP_RAISE("unsupported variable type");
  }
  double time_limit, solution_time_limit;
  int barrier_iterations;
  COPT_CCALL(COPT_GetDblParamDef(copy, COPT_DBLPARAM_TIMELIMIT, &time_limit));
  COPT_CCALL(COPT_GetDblParamDef(copy, COPT_DBLPARAM_SOLTIMELIMIT, &solution_time_limit));
  COPT_CCALL(COPT_GetIntParamDef(copy, COPT_INTPARAM_BARITERLIMIT, &barrier_iterations));
  if (fixed_model_time_limit() >= 0) time_limit = fixed_model_time_limit();
  if (fixed_model_iter_limit() >= 0) barrier_iterations = fixed_model_iter_limit();
  COPT_CCALL(COPT_SetDblParam(copy, COPT_DBLPARAM_TIMELIMIT, time_limit));
  COPT_CCALL(COPT_SetDblParam(copy, COPT_DBLPARAM_SOLTIMELIMIT, solution_time_limit));
  COPT_CCALL(COPT_SetIntParam(copy, COPT_INTPARAM_LPMETHOD, fixed_model_iter_limit() < 0 ? 1 : 2));
  COPT_CCALL(COPT_SetIntParam(copy, COPT_INTPARAM_CROSSOVER, 1));
  COPT_CCALL(COPT_SetIntParam(copy, COPT_INTPARAM_BARITERLIMIT, barrier_iterations));
  if (sensitivity())
    COPT_CCALL(COPT_SetIntParam(copy, COPT_INTPARAM_REQSENSITIVITY, 1));
  COPT_CCALL(COPT_Solve(copy));
  int status = 0, basis = 0;
  COPT_CCALL(COPT_GetIntAttr(copy, COPT_INTATTR_LPSTATUS, &status));
  COPT_CCALL(COPT_GetIntAttr(copy, COPT_INTATTR_HASBASIS, &basis));
  if (status != COPT_LPSTATUS_OPTIMAL || !basis)
    MP_RAISE(fmt::format("fixed LP not optimal or no basis (COPT status {})", status));
  FixedLPResult result;
  COPT_CCALL(COPT_GetDblAttr(copy, COPT_DBLATTR_LPOBJVAL, &result.objective));
  result.primal.resize(n);
  result.dual.resize(m);
  result.varstt.resize(n);
  result.constt.resize(m);
  COPT_CCALL(COPT_GetLpSolution(copy, result.primal.data(), nullptr, result.dual.data(), nullptr));
  COPT_CCALL(COPT_GetBasis(copy, result.varstt.data(), result.constt.data()));
  auto convert = [](int status) {
    switch (status) {
    case COPT_BASIS_BASIC: return (int)BasicStatus::bas;
    case COPT_BASIS_LOWER: return (int)BasicStatus::low;
    case COPT_BASIS_UPPER: return (int)BasicStatus::upp;
    case COPT_BASIS_FIXED: return (int)BasicStatus::equ;
    case COPT_BASIS_SUPERBASIC: return (int)BasicStatus::sup;
    default: MP_RAISE(fmt::format("Unknown COPT fixed-LP basis status {}", (int)status));
    }
  };
  for (auto& status : result.varstt) status = convert(status);
  for (auto& status : result.constt) status = convert(status);
  COPT_CCALL(COPT_GetIntAttr(copy, fixed_model_iter_limit() < 0 ?
      COPT_INTATTR_SIMPLEXITER : COPT_INTATTR_BARRIERITER, &result.iterations));
  if (sensitivity())
    fixed_sens_ = ReadSensRaw(copy);
  return result;
}

CoptBackend::SensRaw CoptBackend::ReadSensRaw(copt_prob* prob) {
  SensRaw raw;
  int has = 0;
  COPT_CCALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASSENSITIVITY, &has));
  if (!has)
    return raw;
  int n = 0, m = 0;
  COPT_CCALL(COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &n));
  COPT_CCALL(COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &m));
  auto col = [&](const char* info, std::vector<double>& v) {
    v.resize(n);
    COPT_CCALL(COPT_GetColInfo(prob, info, n, NULL, v.data()));
  };
  auto row = [&](const char* info, std::vector<double>& v) {
    v.resize(m);
    COPT_CCALL(COPT_GetRowInfo(prob, info, m, NULL, v.data()));
  };
  col(COPT_DBLINFO_OBJ, raw.obj);
  col(COPT_DBLINFO_LB, raw.lb);
  col(COPT_DBLINFO_UB, raw.ub);
  col(COPT_DBLINFO_SAOBJLOW, raw.objlo);
  col(COPT_DBLINFO_SAOBJUP, raw.objhi);
  col(COPT_DBLINFO_SALBLOW, raw.lblo);
  col(COPT_DBLINFO_SALBUP, raw.lbhi);
  col(COPT_DBLINFO_SAUBLOW, raw.ublo);
  col(COPT_DBLINFO_SAUBUP, raw.ubhi);
  row(COPT_DBLINFO_LB, raw.rowlb);
  row(COPT_DBLINFO_UB, raw.rowub);
  row(COPT_DBLINFO_SALBLOW, raw.rowlblo);
  row(COPT_DBLINFO_SALBUP, raw.rowlbhi);
  row(COPT_DBLINFO_SAUBLOW, raw.rowublo);
  row(COPT_DBLINFO_SAUBUP, raw.rowubhi);
  raw.valid = true;
  return raw;
}

/// Variable ranges come from the column infos. COPT receives range
/// constraints natively (lb <= body <= ub), so the bound ranges of
/// every row come from the row infos; the side a one-sided row does
/// not have is reported as unbounded, as the Gurobi driver does.
SensRanges CoptBackend::GetSensRanges() {
  SensRanges sensr;
  SensRaw raw;
  if (IsMIP()) {
    if (!HasFixedLP())
      return sensr;
    raw = fixed_sens_;
  } else
    raw = ReadSensRaw(lp());
  if (!raw.valid) {
    AddWarning("SENSITIVITY",
               "COPT did not provide sensitivity information "
               "(no optimal basis, e.g., barrier without crossover)");
    return sensr;
  }
  const double inf = 1e100;
  auto clip = [&](std::vector<double>& v) {
    for (auto& x: v) {
      if (x <= -COPT_INFINITY) x = -inf;
      else if (x >= COPT_INFINITY) x = inf;
    }
  };
  clip(raw.objlo); clip(raw.objhi); clip(raw.lblo); clip(raw.lbhi);
  clip(raw.ublo); clip(raw.ubhi); clip(raw.lb); clip(raw.ub);

  size_t m = raw.rowlb.size();
  std::vector<double> lb(m), ub(m), lblo(m), lbhi(m), ublo(m), ubhi(m),
      rhs(m), rhslo(m), rhshi(m);
  for (size_t i=0; i<m; ++i) {
    bool lbfin = raw.rowlb[i] > -COPT_INFINITY, ubfin = raw.rowub[i] < COPT_INFINITY;
    lb[i] = lbfin ? raw.rowlb[i] : -inf;
    ub[i] = ubfin ? raw.rowub[i] : inf;
    lblo[i] = lbfin ? raw.rowlblo[i] : -inf;
    lbhi[i] = lbfin ? raw.rowlbhi[i] : inf;
    ublo[i] = ubfin ? raw.rowublo[i] : -inf;
    ubhi[i] = ubfin ? raw.rowubhi[i] : inf;
    if (ubfin && !(lbfin && raw.rowlb[i] == raw.rowub[i])) {   // <= and ranges: ub side
      rhs[i] = ub[i]; rhslo[i] = ublo[i]; rhshi[i] = ubhi[i];
    } else {                                                   // >= and ==: lb side
      rhs[i] = lb[i]; rhslo[i] = lblo[i]; rhshi[i] = lbhi[i];
    }
  }
  clip(lblo); clip(lbhi); clip(ublo); clip(ubhi); clip(rhslo); clip(rhshi);

  auto vars = [&](const std::vector<double>& v) {
    return GetValuePresolver().PostsolveGenericDbl( { {v} } ).GetVarValues()();
  };
  auto cons = [&](const std::vector<double>& v) {
    return GetValuePresolver().PostsolveGenericDbl(
          { {}, {{{ CG_Linear, v }}} } ).GetConValues()();
  };
  sensr.varobj = vars(raw.obj);
  sensr.varlb = vars(raw.lb);
  sensr.varub = vars(raw.ub);
  sensr.varobjlo = vars(raw.objlo);
  sensr.varobjhi = vars(raw.objhi);
  sensr.varlblo = vars(raw.lblo);
  sensr.varlbhi = vars(raw.lbhi);
  sensr.varublo = vars(raw.ublo);
  sensr.varubhi = vars(raw.ubhi);
  sensr.conlb = cons(lb);
  sensr.conub = cons(ub);
  sensr.conlblo = cons(lblo);
  sensr.conlbhi = cons(lbhi);
  sensr.conublo = cons(ublo);
  sensr.conubhi = cons(ubhi);
  sensr.conrhs = cons(rhs);
  sensr.conrhslo = cons(rhslo);
  sensr.conrhshi = cons(rhshi);
  return sensr;
}

void CoptBackend::WindupCOPTSolve() { }

void CoptBackend::ReportResults() {
  ReportCOPTResults();
  ConsiderEmulatedFixedModel();
  BaseBackend::ReportResults();
}

void CoptBackend::ReportCOPTResults() {
  SetStatus( GetSolveResult() );
  AddCOPTMessages();
  if (need_multiple_solutions())
    ReportCOPTPool();

  if(storedOptions_.mempeak())
  {
	  double mempeak = getDblAttr(COPT_DBLATTR_MEMPEAK);
	  AddToSolverMessage(fmt::format("Peak memory used: {} MB\n", mempeak));
      ReportDblSuffix({ "peakmem", suf::PROBLEM }, { { mempeak } });
  }
  if (storedOptions_.numericalStats())
  {
      std::vector<std::pair<const char*, std::string>> vec = {
				{ COPT_DBLATTR_MINELEM, "minelem" },                
                { COPT_DBLATTR_MAXELEM, "maxelem" },
		        { COPT_DBLATTR_MINBOUND, "minbound" },
				{ COPT_DBLATTR_MAXBOUND, "maxbound" },
				{ COPT_DBLATTR_MINRHS, "minrhs" },
				{ COPT_DBLATTR_MAXRHS, "maxrhs" },
				{ COPT_DBLATTR_MINCOST, "mincost" },
				{ COPT_DBLATTR_MAXCOST, "maxcost" },
				{ COPT_DBLATTR_MINQELEM, "minqelem" },
				{ COPT_DBLATTR_MAXQELEM, "maxqelem" },
				{ COPT_DBLATTR_MINQLELEM, "minqlelem" },
				{ COPT_DBLATTR_MAXQLELEM, "maxqlelem" },
				{ COPT_DBLATTR_MINQRHS, "minqrhs" },
				{ COPT_DBLATTR_MAXQRHS, "maxqrhs" },
				{ COPT_DBLATTR_MINQCOST, "minqcost" },
				{ COPT_DBLATTR_MAXQCOST, "maxqcost" }
      };

	  for (const auto& [attr, name] : vec) {
		  double val = getDblAttr(attr);
		  ReportDblSuffix({ name, suf::PROBLEM }, { { val } });
	  }
  }

}
std::vector<double> CoptBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
  COPT_CCALL(COPT_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double CoptBackend::getPoolObjective(int i)
{
  double obj;
  COPT_CCALL(COPT_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
bool CoptBackend::IsWithinPoolGap(double obj, double best) const {
  double worsening = (obj - best)
      * (getIntAttr(COPT_INTATTR_OBJSENSE) == COPT_MAXIMIZE ? -1 : 1);
  return worsening <= storedOptions_.poolGapAbs_
      && worsening <= storedOptions_.poolGap_ * std::fabs(best);
}

/// Report the pool solutions found during the search, best first,
/// within sol:poolgap / sol:poolgapabs and up to sol:poollimit of them
void CoptBackend::ReportCOPTPool() {
  if (!IsMIP())
    return;
  int nPool = getIntAttr(COPT_INTATTR_POOLSOLS);
  std::vector<std::pair<double, int> > sols;      // (objective, index)
  for (int i=0; i<nPool; ++i)
    sols.emplace_back(getPoolObjective(i), i);
  if (getIntAttr(COPT_INTATTR_OBJSENSE) == COPT_MAXIMIZE)
    std::sort(sols.rbegin(), sols.rend());
  else
    std::sort(sols.begin(), sols.end());
  int nReported = 0;
  for (const auto& sol: sols) {
    if (nReported >= std::max(1, storedOptions_.nbest_)
        || !IsWithinPoolGap(sol.first, sols.front().first))
      break;
    auto mv = GetValuePresolver().PostsolveSolution(  // only single-objective with pool
          { { getPoolSolution(sol.second) },
            {},                                       // no duals
            std::vector<double>{ sol.first } } );
    ReportIntermediateSolution(
          { mv.GetVarValues()(), mv.GetConValues()(), mv.GetObjValues()() });
    ++nReported;
  }
}


void CoptBackend::AddCOPTMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> CoptBackend::GetSolveResult() {
  namespace sol = mp::sol;

  if (feasrelax_solved_) {       // COPT_FeasRelax() does not set the status attributes
    if (getIntAttr(COPT_INTATTR_HASFEASRELAXSOL))
      return { sol::SOLVED, "feasibility relaxation solved" };
    return { sol::FAILURE, "feasibility relaxation failed" };
  }

  int solstatus = getIntAttr(COPT_INTATTR_HASSOL);
#ifdef COPT_INTATTR_HASLPSOL
  int sollpstatus = getIntAttr(COPT_INTATTR_HASLPSOL);
  int solmipstatus = getIntAttr(COPT_INTATTR_HASMIPSOL);
  if (solstatus != (sollpstatus || solmipstatus))
    AddWarning("COPT HASSOL", "COPT HASSOL attribute not agrees with HAS(LP/MIP)SOL");
  solstatus = (sollpstatus || solmipstatus);
#endif

  switch (getIntAttr(COPT_INTATTR_STATUS))
  {
    case COPT_STATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case COPT_STATUS_LOCAL_OPTIMAL:
      return { sol::OPTIMAL_LOCALLY, "locally optimal solution" };
    case COPT_STATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case COPT_STATUS_LOCAL_INFEASIBLE:
      return { sol::INFEASIBLE_LOCALLY, "locally infeasible problem" };
    case COPT_STATUS_UNBOUNDED:
      if (solstatus)
        return { sol::UNBOUNDED_FEAS, "unbounded problem, feasible solution" };
      return { sol::UNBOUNDED_NO_FEAS, "unbounded problem, no solution" };
    case COPT_STATUS_INF_OR_UNB:
      return { sol::LIMIT_INF_UNB, "infeasible or unbounded problem. " };
    case COPT_STATUS_TIMEOUT:
      if (solstatus)
        return { sol::LIMIT_FEAS_TIME, "time limit, feasible solution" };
      return { sol::LIMIT_NO_FEAS_TIME, "time limit, without a feasible solution" };
    case COPT_STATUS_INTERRUPTED:
      if (solstatus)
        return { sol::LIMIT_FEAS_INTERRUPT, "interrupted, feasible solution" };
      return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted, without a feasible solution" };
    case COPT_STATUS_NODELIMIT:
      if (solstatus)
        return { sol::LIMIT_FEAS_NODES, "node limit, feasible solution" };
      return { sol::LIMIT_NO_FEAS_NODES, "node limit, without a feasible solution" };
    case COPT_STATUS_ITERLIMIT:
      if (solstatus)
        return { sol::LIMIT_FEAS_ITER, "iteration limit, feasible solution" };
      return { sol::LIMIT_NO_FEAS_ITER, "iteration limit, without a feasible solution" };
     case COPT_STATUS_UNFINISHED:
      if (solstatus)
        return { sol::LIMIT_FEAS_ITER, "unfinished, feasible solution" };
      return { sol::LIMIT_NO_FEAS_ITER, "unfinished, without a feasible solution" };
     case COPT_STATUS_UNSTARTED:
       if (solstatus)
         return { sol::UNCERTAIN, "failure, feasible solution" };
       return { sol::FAILURE, "failure (unstarted)" };
     case COPT_STATUS_IMPRECISE:
       return { sol::UNCERTAIN, "solution is imprecise" };
     case COPT_STATUS_NUMERICAL:
       if (solstatus)
         return { sol::UNCERTAIN, "solution returned but error likely" };
       return { sol::NUMERIC, "failure, numeric issues" };
    default:
      return { sol::UNKNOWN, "unknown" };
  }
}


void CoptBackend::FinishOptionParsing() {
  int v = 1;
  GetSolverOption(COPT_INTPARAM_LOGGING, v);       // tech:outlev
  if (!storedOptions_.logFile_.empty())
  {
    /// The log file needs Logging=1; keep the console quiet
    /// unless tech:outlev asked for output
    SetSolverOption(COPT_INTPARAM_LOGTOCONSOLE, v);
    SetSolverOption(COPT_INTPARAM_LOGGING, 1);
    COPT_CCALL(COPT_SetLogFile(lp(), storedOptions_.logFile_.data()));
  }
  set_verbose_mode(v>0);

  /// Native parameters: file, then inline settings, then write
  if (!storedOptions_.paramRead_.empty())
    COPT_CCALL(COPT_ReadParam(lp(), storedOptions_.paramRead_.c_str()));
  for (const auto& prm: storedOptions_.inlineParams_)
    SetNativeParam(prm);
  if (!storedOptions_.paramWrite_.empty())
    COPT_CCALL(COPT_WriteParam(lp(), storedOptions_.paramWrite_.c_str()));
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo lp_values_method[] = {
  { "-1", "Automatic (default)", -1},
  { "1", "Dual simplex", 1},
  { "2", "Barrier. See also alg:gpumode.", 2},
  { "3", "Crossover", 3},
  { "4", "Concurrent (simplex and barrier simultaneously)", 4},
  { "5", "Choose between simplex and barrier automatically", 5},
  { "6", "First-order method (PDLP). See also alg:gpumode.", 6}
};

static const mp::OptionValueInfo concurrentlpmode_values[] = {
    {"-1", "Choose automatically (default)", -1 },
    { "0","CPU only(simplex + barrier).", 0},
    { "1" ,"CPU(simplex + barrier) + GPU/PDLP (on GPU device 0)", 1},
    {"2","CPU(simplex + barrier) + GPU/barrier (on GPU device 0)", 2},
    {"3", "CPU(simplex + barrier) + GPU/PDLP + GPU/barrier (on GPU device 0 and GPU device 1)", 3}
};

static const mp::OptionValueInfo lp_values_gpu[] = {
  {"-1", "Automatic (default). The first-order method (PDLP) will attempt to use the GPU, "
                   "while the barrier method will use the CPU", -1},
  {"0", "Force CPU mode", 0},
  {"1", "Attempt to use the standard GPU mode", 1},
            {"2", "For the barrier method, attempt to use the high-performance GPU mode "
                  "when solving LP problems, which may lead to higher memory usage. "
                  "For the first-order method (PDLP), this is equivalent to gpumode=1 "
                  "(standard GPU mode)", 2}
};

static const mp::OptionValueInfo mip_start_values_[] = {
    { "-1", "Automatic (default)", -1},
    { "0", "Off", 0},
    { "1", "Only load full and feasible MIP starts", 1},
    { "2", " Only load feasible ones (complete partial solutions by solving subMIPs).", 2}
};

static const mp::OptionValueInfo alg_values_level[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Fast", 1},
  { "2", "Normal", 2},
  { "3", "Aggressive", 3}
};

static const mp::OptionValueInfo pre_values_presolve[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Fast", 1},
  { "2", "Normal", 2},
  { "3", "Aggressive", 3},
  { "4", "No limitations: presolve until the model cannot be "
         "modified any more (may be very time-consuming)", 4}
};

static const mp::OptionValueInfo lp_dualprices_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Use Devex pricing algorithm", 0},
  { "1", "Using dual steepest-edge pricing algorithm", 1}
};

static const mp::OptionValueInfo values_loglevel_[] = {
  { "2", "Print basic optimization logs (default)", 2},
  { "3", "Print memory usage information in addition "
          "to basic optimization logs (for MIP problems).", 3}
};

static const mp::OptionValueInfo lp_barorder_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Approximate Minimum Degree (AMD)", 0},
  { "1", "Nested Dissection (ND1)", 1},
  { "2", "Modified Nested Dissection (ND2)", 2}
};

static const mp::OptionValueInfo values_nonconvex[] = {
    { "-1", "Default choice (as of COPT 8, same as 1)", -1},
    { "0", "Report nonconvexity and terminate", 0},
    { "1", "Search for a locally optimal solution", 1},
    { "2", "Search for a globally optimal solution.", 2}
};

static const mp::OptionValueInfo values_iismethod[] = {
  {"-1", "Automatic choice (default)", -1},
  { "0", "Find smaller IIS", 0},
  { "1", "Find IIS quickly", 1},
};

static const mp::OptionValueInfo values_crossover[] = {
  {"-1", "Choose automatically - only run crossover when LP solution is not primal-dual feasible", -1},
  {"0", "No", 0},
  {"1", "Yes", 1}
};


static const mp::OptionValueInfo values_nlpmuupdate[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "Monotonically decreasing", 0},
  {"1", "Adaptive adjustment", 1}
};

static const mp::OptionValueInfo values_nlplinscale[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "No scaling", 0},
  {"1", "Always scale", 1}
};

static const mp::OptionValueInfo values_sdpmethod[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "Primal-dual method", 0},
  {"1", "ADMM", 1},
  {"2", "Dual method", 2}
};

static const mp::OptionValueInfo values_barstart[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "Simple starting point", 0},
  {"1", "Mehrotra starting point", 1},
  {"2", "Modified Mehrotra starting point", 2}
};

static const mp::OptionValueInfo values_level_heavy[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Fast", 1},
  { "2", "Normal", 2},
  { "3", "Aggressive", 3},
  { "4", "Heavy", 4}
};

static const mp::OptionValueInfo values_tunemethod[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "Greedy search", 0},
  {"1", "Broader search", 1}
};

static const mp::OptionValueInfo values_tunemode[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "Solving time", 0},
  {"1", "Optimality gap reached", 1},
  {"2", "Objective value", 2},
  {"3", "Objective bound", 3}
};

static const mp::OptionValueInfo values_tunemeasure[] = {
  {"-1", "Automatic (default)", -1},
  {"0", "Average over the permutations", 0},
  {"1", "Maximum over the permutations", 1}
};

static const mp::OptionValueInfo values_tuneoutlev[] = {
  {"0", "No tuning log", 0},
  {"1", "Summary of improved parameter sets only", 1},
  {"2", "Summary of each attempt (default)", 2},
  {"3", "Detailed log of each attempt", 3}
};

static const mp::OptionValueInfo values_pool_mode[] = {
  { "0", "Just collect solutions during normal solve", 0},
  { "1", "Make some effort at finding additional solutions "
         "(not supported by COPT; same as 0)", 1},
  { "2", "Seek \"poollimit\" best solutions (default; not supported "
         "by COPT, same as 0). The solutions reported are limited by "
         "poolgap and poolgapabs.", 2}
};

static const mp::OptionValueInfo miprepair_values_[] = {
  {"-1", "Only when time left (default)", -1},
  {"0", "Off", 0},
  {"1", "Extend time limit and attempt repair (Fast)", 1},
  {"2", "Extend time limit and attempt repair (Normal)", 2},
  {"3", "Extend time limit and attempt repair (Aggressive)", 3}
};


void CoptBackend::InitCustomOptions() {

  set_option_header(
      "COPT Optimizer Options for AMPL\n"
      "-------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``copt_options``. For example::\n"
      "\n"
      "  ampl: option copt_options 'mipgap=1e-6';\n");

  AddToOptionDescription(
      "tech:writesolution",
      "File name extensions can be "
      "``.sol``, ``.bas``, ``.iis``.");

  AddSolverOption("tech:outlev outlev",
      "0-1: output logging verbosity. "
      "Default = 0 (no logging).",
    COPT_INTPARAM_LOGGING, 0, 1);

  AddSolverOption("tech:loglevel loglevel",
    "Controls the level of detail in the optimization logs:\n"
    "\n.. value-table::\n", COPT_INTPARAM_LOGLEVEL,
	  values_loglevel_, 2);

  AddStoredOption("tech:logfile logfile",
    "Log file name.", storedOptions_.logFile_);

  AddListOption("tech:optionnative optionnative optnative tech:param param",
    "General way to specify values of both documented and "
    "undocumented COPT parameters; value should be a quoted "
    "string (delimited by ' or \") containing a parameter name, a "
    "space (or '='), and the value to be assigned to the parameter, "
    "for example \"RelGap 0.01\".  Can appear more than once. "
    "Cannot be used to query current parameter values.",
    storedOptions_.inlineParams_);

  AddStoredOption("tech:optionnativeread optionnativeread tech:param:read param:read",
    "Name of a COPT parameter file to read (surrounded by 'single' or "
    "\"double\" quotes if the name contains blanks). Each nonempty "
    "line contains a parameter name and its value separated by blanks, "
    "as written by tech:optionnativewrite.",
    storedOptions_.paramRead_);

  AddStoredOption("tech:optionnativewrite optionnativewrite tech:param:write param:write",
    "Name of a COPT parameter file (surrounded by 'single' or \"double\" "
    "quotes if the name contains blanks) to be written with the "
    "parameters that differ from their defaults.",
    storedOptions_.paramWrite_);

  AddStoredOption("tech:mempeak mempeak",
      "Return the peak memory used by the solver in the problem suffix "
      "\"peakmem\". Default = 0 (do not return).", storedOptions_.mempeak_);

  AddStoredOption("tech:numericalranges numericalranges numerical_ranges",
	  "Return the numerical ranges of the problem as problem suffixes: minelem, maxelem, "
	  "minbound, maxbound, minrhs, maxrhs, mincost, maxcost, minqelem, maxqelem, minqlelem, "
      "maxqlelem, minqrhs, maxqrhs, minqcost, maxqcost. See CHANGES.copt.md for details. "
      "Default = 0 (do not return).", 
      storedOptions_.numericalStats_);

  AddStoredOption("tech:tunebase tunebase",
    "Base name for results of running COPT's search for best "
    "parameter settings.  The search is run only when tunebase "
    "is specified.  Results are written to files with names derived "
    "from tunebase by appending \".par\" if \".par\" does not occur in "
    "tunebase and inserting _1_, _2_, ... (for the first, second, "
    "... set of parameter settings) before the right-most \".par\". "
    "The file with _1_ inserted is the best set and the solve "
    "results returned are for this set.  In a subsequent \"solve;\", "
    "you can use \"tech:optionnativeread=...\" to apply the settings "
    "in results file ... .",
    storedOptions_.tunebase_);

  AddSolverOption("tech:tunetimelim tunetimelim lim:tunetime",
    "Time limit (in seconds) on tuning when tunebase is specified; "
    "default 0 ==> chosen by COPT.",
    COPT_DBLPARAM_TUNETIMELIMIT, 0.0, 1e20);

  AddSolverOption("tech:tunetargettime tunetargettime",
    "Target solve time (in seconds) for tuning: the tuner stops "
    "once a parameter set reaches it (default 0.01).",
    COPT_DBLPARAM_TUNETARGETTIME, 0.0, DBL_MAX);

  AddSolverOption("tech:tunetargetrelgap tunetargetrelgap",
    "Target relative MIP gap for tuning (default 1e-4).",
    COPT_DBLPARAM_TUNETARGETRELGAP, 0.0, DBL_MAX);

  AddSolverOption("tech:tunemethod tunemethod",
    "Tuning method:\n"
    "\n.. value-table::\n", COPT_INTPARAM_TUNEMETHOD,
    values_tunemethod, -1);

  AddSolverOption("tech:tunemode tunemode",
    "Tuning criterion:\n"
    "\n.. value-table::\n", COPT_INTPARAM_TUNEMODE,
    values_tunemode, -1);

  AddSolverOption("tech:tunemeasure tunemeasure",
    "How tuning results over several permutations are aggregated:\n"
    "\n.. value-table::\n", COPT_INTPARAM_TUNEMEASURE,
    values_tunemeasure, -1);

  AddSolverOption("tech:tunepermutes tunepermutes",
    "Number of permutations of the model tried for each parameter set; "
    "default 0 ==> chosen by COPT.",
    COPT_INTPARAM_TUNEPERMUTES, 0, INT_MAX);

  AddSolverOption("tech:tuneoutlev tuneoutlev tuneoutput",
    "Amount of tuning output when tunebase is specified:\n"
    "\n.. value-table::\n", COPT_INTPARAM_TUNEOUTPUTLEVEL,
    values_tuneoutlev, 2);

  AddOptionSynonyms_Inline_Front("ams_stub", "sol:stub");

  AddStoredOption("sol:poollimit ams_limit poollimit poolnbest solnlimit",
    "Limit on the number of alternative MIP solutions reported "
    "(see sol:stub, sol:count). Default 10. COPT reports the "
    "solutions found during the search and cannot look for more.",
    storedOptions_.nbest_);

  AddStoredOption("sol:poolmode ams_mode poolmode",
    "Search mode for MIP solutions when sol:stub/sol:count are specified "
    "to request finding several alternative solutions. COPT only "
    "collects the solutions found during the search, so all values "
    "behave alike; the option is accepted for compatibility:\n"
    "\n.. value-table::\n",
    storedOptions_.nPoolMode_, values_pool_mode);

  AddStoredOption("sol:poolgap ams_eps poolgap",
    "Relative tolerance for reporting alternate MIP solutions "
    "(default: no limit): a solution is reported if its objective "
    "is within poolgap * |best objective| of the best one.",
    storedOptions_.poolGap_);

  AddStoredOption("sol:poolgapabs ams_epsabs poolgapabs",
    "Absolute tolerance for reporting alternate MIP solutions "
    "(default: no limit).",
    storedOptions_.poolGapAbs_);

  AddSolverOption("lp:dualprice dualprice",
    "Specifies the dual simplex pricing algorithm:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DUALPRICE,
    lp_dualprices_values_, -1);

  AddSolverOption("lp:dualperturb dualperturb",
    "Whether to allow the objective function perturbation when using "
    "the dual simplex method:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DUALPERTURB,
    values_autonoyes_, -1);

  AddSolverOption("lp:barhomogeneous barhomogeneous",
    "Whether to use homogeneous self-dual form in barrier:\n"
    "\n.. value-table::\n", COPT_INTPARAM_BARHOMOGENEOUS,
    values_autonoyes_, -1);

  AddSolverOption("bar:start barstart",
    "Algorithm for the barrier starting point:\n"
    "\n.. value-table::\n", COPT_INTPARAM_BARSTART,
    values_barstart, -1);

  AddSolverOption("lp:barorder barorder",
    "Barrier ordering algorithm:\n"
    "\n.. value-table::\n", COPT_INTPARAM_BARORDER,
    lp_barorder_values_, -1);

  AddSolverOption("tech:seed seed",
                  "Random number seed, can influence the solution path "
                  "(COPT parameter RandSeed, not in the vendor documentation). "
                  "Default -1.", "RandSeed",
                  -1, INT_MAX);

  
  AddSolverOption("bar:iterlim bariterlim BarIterLimit",
    "Limit on the number of barrier iterations (default 500).",
    COPT_INTPARAM_BARITERLIMIT, 0, INT_MAX);

  AddSolverOption("mip:cutlevel cutlevel",
    "Level of cutting-planes generation:\n"
    "\n.. value-table::\n", COPT_INTPARAM_CUTLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:intfeastol intfeastol inttol",
    "Feasibility tolerance for integer variables (default 1e-06).",
    COPT_DBLPARAM_INTTOL, 1e-9, 0.1);

  AddSolverOption("mip:rootcutlevel rootcutlevel",
    "Level of cutting-planes generation of root node:\n"
    "\n.. value-table::\n", COPT_INTPARAM_ROOTCUTLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:treecutlevel treecutlevel",
    "Level of cutting-planes generation of search tree:\n"
    "\n.. value-table::\n", COPT_INTPARAM_TREECUTLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:repair miprepair repair",
      "Level for repairing the MIP solution in case of "
      "numerical issues:\n"
      "\n.. value-table::\n", COPT_INTPARAM_MIPREPAIR,
      miprepair_values_, -1);
      
  AddSolverOption("mip:rootcutrounds rootcutrounds",
    "Rounds of cutting-planes generation of root node;\n" 
    "default -1 ==> automatic.",
    COPT_INTPARAM_ROOTCUTROUNDS, -1, INT_MAX);

  AddSolverOption("mip:nodecutrounds nodecutrounds",
    "Rounds of cutting-planes generation of search tree node;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_NODECUTROUNDS, -1, INT_MAX);

  AddSolverOption("mip:heurlevel heurlevel",
    "Level of heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_HEURLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:prerootheurlevel prerootheurlevel",
      "Level of pre-root heuristics:\n"
      "\n.. value-table::\n", COPT_INTPARAM_PREROOTHEURLEVEL,
      alg_values_level, -1);


  AddSolverOption("mip:roundingheurlevel roundingheurlevel",
    "Level of rounding heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_ROUNDINGHEURLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:divingheurlevel divingheurlevel",
    "Level of diving heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DIVINGHEURLEVEL,
    alg_values_level, -1);

  AddSolverOption("mip:startmode mipstartmode",
                  "Mode of MIP starts:\n"
                  "\n.. value-table::\n"
                  "COPT 8.0 also documents 3 = accept all starts "
                  "(for hierarchical multi-objective models), but the "
                  "8.0.6 library limits the value to 2.", COPT_INTPARAM_MIPSTARTMODE,
                  mip_start_values_, -1);

  AddSolverOption("mip:fapheurlevel fapheurlevel",
    "Level of fix-and-propagate heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_FAPHEURLEVEL,
    values_level_heavy, -1);

  AddSolverOption("mip:startnodelimit mipstartnodelimit",
    "Node limit for the sub-MIPs completing partial MIP starts "
    "(see mip:startmode); default -1 ==> automatic.",
    COPT_INTPARAM_MIPSTARTNODELIMIT, -1, INT_MAX);

  AddSolverOption("mip:linearizeindicators linearizeindicators",
    "Whether to force the linearization of indicator constraints "
    "into equivalent linear constraints:\n"
    "\n.. value-table::\n", COPT_INTPARAM_LINEARIZEINDICATORS,
    values_01_noyes_0default_, 0);

  AddSolverOption("mip:linearizesos linearizesos",
    "Whether to force the linearization of SOS constraints "
    "into equivalent linear constraints:\n"
    "\n.. value-table::\n", COPT_INTPARAM_LINEARIZESOS,
    values_01_noyes_0default_, 0);

  AddSolverOption("mip:submipheurlevel submipheurlevel",
    "Level of Sub-MIP heuristics:\n"
    "\n.. value-table::\n", COPT_INTPARAM_SUBMIPHEURLEVEL,
    alg_values_level, -1);


  AddSolverOption("mip:strongbranching strongbranching",
    "Level of strong branching:\n"
    "\n.. value-table::\n", COPT_INTPARAM_STRONGBRANCHING,
    alg_values_level, -1); 

    AddSolverOption("mip:conflictanalysis conflictanalysis",
      "Whether to perform conflict analysis:\n"
      "\n.. value-table::\n", COPT_INTPARAM_CONFLICTANALYSIS,
      values_autonoyes_, -1);

  AddSolverOption("mip:gap mipgap",
                  "Relative MIP optimality gap, default 1e-4.\n",
                  COPT_DBLPARAM_RELGAP, 0.0, DBL_MAX);

    AddSolverOption("mip:gapabs mipgapabs",
                    "Absolute MIP optimality gap, default 1e-6.\n",
                    COPT_DBLPARAM_ABSGAP, 0.0, DBL_MAX);

  AddSolverOption("pre:dualize dualize",
    "Whether to dualize the problem before solving it:\n"
    "\n.. value-table::\n", COPT_INTPARAM_DUALIZE,
    values_autonoyes_, -1);

  AddSolverOption("pre:solve presolve",
    "Level of presolving to perform before solving the problem:\n"
    "\n.. value-table::\n", COPT_INTPARAM_PRESOLVE,
    pre_values_presolve, -1);

  AddSolverOption("pre:scale scale",
    "Whether to scale the problem:\n"
    "\n.. value-table::\n"
    "Scaling typically reduces solution times, but it may lead "
    "to larger constraint violations in the original, unscaled "
    "model. Choosing a different scaling option can sometimes "
    "improve performance for particularly numerically difficult "
    "models.",
    COPT_INTPARAM_SCALING, values_autonoyes_, -1);

  AddSolverOption("tech:threads threads",
    "Number of threads to use;\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_THREADS, -1, 128);

  AddSolverOption("tech:barrierthreads barthreads",
      "Number of threads used by the barrier algorithm;\n"
      "default -1 ==> use value in tech:threads.",
    COPT_INTPARAM_BARTHREADS, -1, 128);

  AddSolverOption("tech:crossoverthreads crossoverthreads",
    "Number of threads used by crossover;\n"
    "default -1 ==> use value in tech:threads.",
    COPT_INTPARAM_CROSSOVERTHREADS, -1, 128);

  AddSolverOption("tech:simplexthreads simplexthreads",
    "Number of threads used by dual simplex;\n"
    "default -1 ==> use value in tech:threads.",
    COPT_INTPARAM_SIMPLEXTHREADS, -1, 128);


  AddSolverOption("tech:miptasks miptasks",
    "Number of MIP tasks in parallel (at most 256);\n"
    "default -1 ==> automatic.",
    COPT_INTPARAM_MIPTASKS, -1, 256);

  
  AddSolverOption("lim:mipnlpiterlimit mipnlpiterlimit",
      "Iteration limit for solving NLP problem(s) within the MIP solver; "
      "-1 = no limit. Default 100.",
      COPT_INTPARAM_MIPNLPITERLIMIT, -1, INT_MAX);

  AddSolverOption("lim:nlpiterlimit nlpiterlimit",
    "Iteration limit for the nonlinear solver (default 10000).",
	  COPT_INTPARAM_NLPITERLIMIT, 0, INT_MAX);

  AddSolverOption("lim:time timelim timelimit",
      "limit on solve time (in seconds; default: no limit).",
      COPT_DBLPARAM_TIMELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:soltime soltimelim soltimelimit",
    "Time limit if a primal feasible solution has been found "
    "(in seconds; default: no limit).",
    COPT_DBLPARAM_SOLTIMELIMIT, 0.0, DBL_MAX);

  AddSolverOption("lim:nodes nodelim nodelimit",
    "Maximum MIP nodes to explore; default -1 (no limit).",
    COPT_INTPARAM_NODELIMIT, -1, INT_MAX);

  AddSolverOption("lim:multiobjtime multiobjtimelim multiobjtimelimit obj:multi:timelimit",
    "Time limit (in seconds) for the whole multi-objective solve "
    "(default: no limit).",
    COPT_DBLPARAM_MULTIOBJTIMELIMIT, 0.0, 1e20);

  AddSolverOption("lim:memlimit memlimit",
      "Soft limit (number of MB) on memory allocated; "
      "default 1e20 (no limit).",
      COPT_DBLPARAM_MEMLIMIT, 0.0, Infinity());


  AddSolverOption("alg:method method lp:method lpmethod",
    "Which algorithm to use for non-MIP problems:\n"
    "\n.. value-table::\n", COPT_INTPARAM_LPMETHOD, 
    lp_values_method, -1);

  AddSolverOption("lp:concurrentmode lp:concurrentlpmode concurrentlpmode",
      "The LP concurrent solving mode, only effective when alg:method = 4. "
      "The parameters alg:gpumode and alg:gpudevice are ignored, "
                  "and GPU usage and "
      "device selection are fully controlled by this parameter.",
                  COPT_INTPARAM_CONCURRENTLPMODE,
      concurrentlpmode_values, -1);

  AddSolverOption("alg:gpumode gpumode lp:pdlpgpumode pdlpgpumode bar:gpumode",
    "Whether to use GPU or CPU for barrier and PDLP methods. Note that CUDA "
    "GPU mode is only supported on Windows and Linux:\n"
    "\n.. value-table::\n", COPT_INTPARAM_GPUMODE,
    lp_values_gpu, -1);

  AddSolverOption("alg:gpudevice gpudevice lp:pdlpgpudevice pdlpgpudevice bar:gpudevice",
    "Specify devide ID of GPU to use in case of multiple GPUs "
    "(default -1, choose automatically)",
    COPT_INTPARAM_GPUDEVICE, -1, INT_MAX);

  AddSolverOption("lp:pdlptol pdlptol",
    "Convergence tolerance for PDLP (default 1e-6).",
    COPT_DBLPARAM_PDLPTOL, 1e-12, 1e-4);

  AddSolverOption("alg:nonconvex nonconvex",
                  "Strategy for continuous nonconvex models (default -1):\n"
                  "\n.. value-table::\n",
                  COPT_INTPARAM_NONCONVEX, values_nonconvex, -1);

  AddSolverOption("alg:nlptol nlptol",
    "Convergence tolerance of the nonlinear solver (default 1e-8).",
    COPT_DBLPARAM_NLPTOL, 1e-13, 1e-3);

  AddSolverOption("alg:nlpmuupdate nlpmuupdate",
    "Barrier parameter update strategy of the nonlinear solver:\n"
    "\n.. value-table::\n", COPT_INTPARAM_NLPMUUPDATE,
    values_nlpmuupdate, -1);

  AddSolverOption("alg:nlplinscale nlplinscale",
    "Linear system scaling strategy of the nonlinear solver:\n"
    "\n.. value-table::\n", COPT_INTPARAM_NLPLINSCALE,
    values_nlplinscale, -1);

  AddSolverOption("alg:sdpmethod sdpmethod",
    "Method for semidefinite programs:\n"
    "\n.. value-table::\n", COPT_INTPARAM_SDPMETHOD,
    values_sdpmethod, -1);

  AddSolverOption("alg:feastol feastol",
    "Primal feasibility tolerance (default 1e-6).",
    COPT_DBLPARAM_FEASTOL, 1e-9, 1e-4);

  AddSolverOption("alg:dualfeastol dualfeastol",
    "Tolerance for dual solutions and reduced cost (default 1e-6).",
    COPT_DBLPARAM_DUALTOL, 1e-9, 1e-4);

  AddSolverOption("alg:matrixtol matrixtol",
    "input matrix coefficient tolerance (default 1e-10).",
     COPT_DBLPARAM_MATRIXTOL, 0.0, 1e-7);
    AddSolverOption("alg:qmatrixtol qmatrixtol",
        "input Q matrix coefficient tolerance (default 5e-14).",
    COPT_DBLPARAM_QMATRIXTOL, 0.0, 1e-7);

  AddSolverOption("tech:writedualprob writedualprob",
    "Whether tech:writeprob writes the dual problem instead of "
    "the original one for LP in the LP, MPS and BIN formats:\n"
    "\n.. value-table::\n", COPT_INTPARAM_WRITEDUALPROB,
    values_01_noyes_0default_, 0);

  AddSolverOption("iis:method iismethod alg:iismethod",
    "Which method to use when finding an IIS (irreducible infeasible "
    "set of constraints, including variable bounds):\n"
    "\n.. value-table::\n",
    COPT_INTPARAM_IISMETHOD, values_iismethod, -1);

  AddSolverOption("bar:crossover crossover",
    "Execute crossover to transform a barrier solution to a basic one (default 1).",
    COPT_INTPARAM_CROSSOVER, values_crossover, 1);

  ReplaceOptionDescription("alg:feasrelax",
    "Whether to modify the problem into a feasibility "
    "relaxation problem:\n"
    "\n"
    "| 0 = No (default)\n"
    "| 1 = Yes, minimizing the weighted sum of violations\n"
    "| 2 = Yes, minimizing the weighted sum of squared violations\n"
    "| 3 = Yes, minimizing the weighted count of violations\n"
    "| 4-6 = Same objective as 1-3, but also optimize the "
    "original objective, subject to the violation "
    "objective being minimized.\n"
    "\n"
    "COPT does not support modes 2 and 5 (squared violations) "
    "for models with integer variables. "
    "Weights are given by suffixes .lbpen and .ubpen on variables "
    "and .rhspen on constraints (when nonnegative), else by keywords "
    "alg:lbpen, alg:ubpen, and alg:rhspen, respectively (default values = 1). "
    "Weights < 0 are treated as Infinity, allowing no violation.");
}


double CoptBackend::MIPGap() {
  // TODO Check what happens if not solved
  return getDblAttr(COPT_DBLATTR_BESTGAP);
}
double CoptBackend::BestDualBound() {
  // TODO Check what happens if not solved
  return getDblAttr(COPT_DBLATTR_BESTBND);
}

double CoptBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> CoptBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  COPT_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case COPT_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case COPT_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case COPT_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case COPT_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case COPT_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Copt VBasis value: {}", s));
    }
  }
  return vars;
}

ArrayRef<int> CoptBackend::ConStatii() {
  std::vector<int> cons(NumLinCons());
  COPT_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case COPT_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case COPT_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case COPT_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case COPT_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case COPT_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Copt VBasis value: {}", s));
    }
  }
  return cons;
}

void CoptBackend::SetCoptBasis(
    ArrayRef<int> vst, ArrayRef<int>cst) {
  int index[1];
  std::vector<int> vstt(vst.data(), vst.data() + vst.size());
  for (auto j = vstt.size(); j--; ) {
    auto& s = vstt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = COPT_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = COPT_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = COPT_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = COPT_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = COPT_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!COPT_GetColInfo(lp(), COPT_DBLINFO_LB, 1, index, &lb) && 
        !COPT_GetColInfo(lp(), COPT_DBLINFO_UB, 1, index, &ub))
      { 
        if (lb >= -1e-6)
          s = COPT_BASIS_LOWER;
        else if (ub <= 1e-6)
          s = COPT_BASIS_UPPER;
        else
          s = COPT_BASIS_SUPERBASIC;  // or, leave at 0?
      }
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }
  std::vector<int> cstt(cst.data(), cst.data() + cst.size());
  for (int j = (int)cstt.size(); j--; ) {
    auto& s = cstt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = COPT_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = COPT_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = COPT_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = COPT_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = COPT_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      if(!COPT_GetRowInfo(lp(), COPT_DBLINFO_LB, 1, &j, &lb) &&
          !COPT_GetRowInfo(lp(), COPT_DBLINFO_UB, 1, &j, &ub))
      {
        if (lb >= -1e-6)
          s = COPT_BASIS_LOWER;
        else if (ub <= 1e-6)
          s = COPT_BASIS_UPPER;
        else
          s = COPT_BASIS_SUPERBASIC;  // or, leave at 0?
      }
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  COPT_SetBasis(lp(), vstt.data(), cstt.data());
}

SolutionBasis CoptBackend::GetBasis() {
  if (IsMIP()) return FixedLPBasis(CG_Linear);
  if (!getIntAttr(COPT_INTATTR_HASBASIS)) return {};
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void CoptBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  SetCoptBasis(varstt, constt);
}


void CoptBackend::ComputeIIS() {
  COPT_CCALL(COPT_ComputeIIS(lp()));
  SetStatus(GetSolveResult());   // could be new information
}

IIS CoptBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}


std::vector<int> getIIS(copt_prob* prob, int num, 
    int(*getlb)(copt_prob*, int, const int*, int*),
  int(*getub)(copt_prob*, int, const int*, int*)) {
  std::vector<int> iis_lb(num);
  std::vector<int> iis_ub(num);
  COPT_CCALL(getlb(prob, num, NULL, iis_lb.data()));
  COPT_CCALL(getub(prob, num, NULL, iis_ub.data()));
  for (size_t i = iis_lb.size(); i--; ) {
    if (iis_ub[i]) {
      if (iis_lb[i])
        iis_lb[i] = (int)IISStatus::fix;
      else
        iis_lb[i] = (int)IISStatus::upp;
    }
    else {
      if (iis_lb[i])
        iis_lb[i] = (int)IISStatus::low;
      else
        iis_lb[i] = (int)IISStatus::non;
    }
  }
  return iis_lb;
}

static void ConvertIIS2AMPL(std::vector<int>& ai) {
  for (int i = ai.size(); i--; ) {
    ai[i] = int(ai[i] ? IISStatus::mem : IISStatus::non);
  }
}


ArrayRef<int> CoptBackend::VarsIIS() {
  return getIIS(lp(), NumVars(), COPT_GetColLowerIIS, COPT_GetColUpperIIS);
}
pre::ValueMapInt CoptBackend::ConsIIS() {
  auto iis_lincon =  getIIS(lp(), NumLinCons(), COPT_GetRowLowerIIS, COPT_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  COPT_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  COPT_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };

}


ArrayRef<double> CoptBackend::Ray() {
  if (!getIntAttr(COPT_INTATTR_HASPRIMALRAY))   // e.g., for MIP
    return {};
  auto uray_pres = getVarInfo(COPT_DBLINFO_PRIMALRAY);
  auto mv = GetValuePresolver().PostsolveSolution({ uray_pres });
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> CoptBackend::DRay() {
  if (!getIntAttr(COPT_INTATTR_HASDUALFARKAS))  // e.g., for MIP
    return {};
  auto fd = getConInfo(COPT_DBLINFO_DUALFARKAS);
  auto vm = GetValuePresolver().PostsolveSolution({
                                               {},
                                               {{{CG_Linear, std::move(fd)}}}
    });
  return vm.GetConValues().MoveOut();        // need the vector itself
}


void CoptBackend::AddPrimalDualStart(Solution sol0_unpres) {
  if (IsMIP())
    return;
  auto mv = GetValuePresolver().PresolveSolution(
    { sol0_unpres.primal, sol0_unpres.dual });
  int nvars = NumVars();
  int ncons = NumLinCons();

  auto x0 = mv.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);
  std::vector<double> nv(x0.size());
  std::vector<double> ne(pi0.size());
  COPT_CCALL(COPT_SetLpSolution(lp(),
                                x0.data(),
                                ne.data(),
                                pi0.data(),
                                nv.data()));
}

void CoptBackend::AddMIPStart(
		ArrayRef<double> x0_unpres, ArrayRef<int> sparsity_unpres) {
  auto mv = GetValuePresolver().PresolveSolution({ x0_unpres });
  auto ms = GetValuePresolver().PresolveGenericInt({ sparsity_unpres });
  auto x0 = mv.GetVarValues()();
  auto s0 = ms.GetVarValues()();
  std::vector<int> idx;                 // Create sparse vector
  idx.reserve(x0.size());
  std::vector<double> val;
  val.reserve(x0.size());
  for (int i = 0; i < (int)x0.size(); ++i) {
    if (s0[i]) {
      idx.push_back(i);
      val.push_back(x0[i]);
    }
  }
  COPT_CCALL(COPT_AddMipStart(lp(), val.size(), idx.data(), const_cast<double*>(val.data())));
  /// @todo Only do this if we have expressions
  /// @note AddMIPStart() is called when the model can be non-LP
  COPT_CCALL(COPT_SetNLPrimalStart(lp(), idx.size(), idx.data(), val.data()));
}


} // namespace mp

  // AMPLs
AMPLS_MP_Solver* Open_copt(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::CoptBackend()}, 
    cb);
}

void AMPLSClose_copt(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_copt(AMPLS_MP_Solver* slv) {
  return dynamic_cast<mp::CoptBackend*>(AMPLSGetBackend(slv))->lp();
}

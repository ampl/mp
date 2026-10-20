#include <vector>
#include <climits>
#include <cfloat>
#include <cmath>
#include <algorithm>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "highsmpbackend.h"

extern "C" {
  #include "highsmp-ampls-c-api.h"    // Highs AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptHighs(void* prob) {
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateHighsBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()};
}


namespace mp {

/// Create Highs Model Manager
/// @param gc: the Highs common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return HighsModelMgr
std::unique_ptr<BasicModelManager>
CreateHighsModelMgr(HighsCommon&, Env&, pre::BasicValuePresolver*&);


HighsBackend::HighsBackend() {
  LoadHighsLibrary(false);
  OpenSolver(); 
  pre::BasicValuePresolver* pPre;
  auto data = CreateHighsModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);
}

HighsBackend::~HighsBackend() {
  CloseSolver();
}


const char* HighsBackend::GetBackendName()
  { return "HighsBackend"; }

std::string HighsBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", HIGHS_VERSION_MAJOR, 
    HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH);
}

void HighsBackend::InputStartValues() {
  // Set the accumulated objectives
  // In case of native MO, accObjectives will contain all the problem's objectives - and their meta into
  // and will have the flag hadMultiObjective set to true
  accObjectives().setAllInHighs(lp());
  accObjectives().clear();
  // Changing objective data invalidates a HiGHS start: install it afterwards.
  BaseBackend::InputStartValues();
}
  void HighsBackend::FinishOptionParsing() {
    // If CUDA is specified, scrap the current library, load the cuda one
    // and replay all the options

      if (storedOptions_.onGPU())
      {
    #ifdef __APPLE__
          throw std::runtime_error("GPU support is not available on MacOS");
    #endif
      loader().Highs_destroy(lp());
      LoadHighsLibrary(true);
      OpenSolver();
      ReplaySolverOptions();
    }

    std::string method = storedOptions_.lpmethod_;
    if (method == "pdlp-gpu") method = "pdlp";
    else if (method == "hipdlp-gpu") method = "hipdlp";
    try {
      SetSolverOption("solver", method);
    } catch (const std::exception& e) {
      (void)e;
      MP_RAISE(fmt::format(
          "alg:method={} is not available in this HiGHS library build "
          "(see the HiGHS message above)", storedOptions_.lpmethod_));
    }
    /// Copy env/lp to ModelAPI
    copy_common_info_to_other();
    int v = -1;
    GetSolverOption("output_flag", v);
    set_verbose_mode(v > 0);
  }



bool HighsBackend::IsQCP() const {
  return false; 
}

ArrayRef<double> HighsBackend::PrimalSolution() {
  if (HasFixedLP()) return FixedLP().primal;
  if (IsMIP() && getIntAttr("primal_solution_status") != kHighsSolutionStatusFeasible)
    return {};
  int num_vars = NumVars();
  std::vector<double> x(num_vars);
  // int primal_solution_status;
  // loader().Highs_getIntInfoValue(lp(),
  //                       "primal_solution_status", &primal_solution_status);
  // if ((kHighsSolutionStatusFeasible == primal_solution_status) || (storedOptions_.onGPU()))
  auto e = loader().Highs_getSolution(lp(), x.data(), NULL, NULL, NULL);
  // else
  if (e != kHighsStatusOk && e != kHighsStatusWarning)
    x.clear();
  return x;
}

pre::ValueMapDbl HighsBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> HighsBackend::DualSolution_LP() {
  if (IsMIP()) return HasFixedLP() ? ArrayRef<double>(FixedLP().dual) : ArrayRef<double>{};
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
  int error = loader().Highs_getSolution(lp(), NULL, NULL, NULL, pi.data());
  if (error)
    pi.clear();
  return pi;
}

  ArrayRef<double> HighsBackend::GetObjectiveValues() {
    if (HasFixedLP()) return std::vector<double>{FixedLP().objective};
    if (accObjectives().hadNativeMultiObj()) {
      AddToSolverMessage("Warning: HiGHS does not support returning the objective values for multiple\n"
        "objectives; AMPL will compute the actual values.");
      return std::vector<double>();
    }
    return std::vector<double> { ObjectiveValue() };
}

  void HighsBackend::ObjPriorities(ArrayRef<int> pri) {
    if (pri.size() > 0)
      SetSolverOption("blend_multi_objectives",
        0);
    accObjectives().setPriorities(pri);

}
  void HighsBackend::ObjWeights(ArrayRef<double> w) {
    accObjectives().setWeights(w);
}
  void HighsBackend::ObjAbsTol(ArrayRef<double> a) {
    accObjectives().setAbsTols(a);
  }
  void HighsBackend::ObjRelTol(ArrayRef<double> r) {
    accObjectives().setRelTols(r);
  }

double HighsBackend::ObjectiveValue() const {
  return loader().Highs_getObjectiveValue(lp());
}

int HighsBackend::NodeCount() const {
  return getInt64Attr("mip_node_count");
}

int HighsBackend::SimplexIterations() const {
  return getIntAttr("simplex_iteration_count");
}

int HighsBackend::BarrierIterations() const {
  return getIntAttr("ipm_iteration_count");
}
int HighsBackend::PdlpIterations() const {
  return getIntAttr("pdlp_iteration_count");
}

void HighsBackend::DoWriteProblem(const std::string &file) {
  HIGHS_CCALL(loader().Highs_writeModel(lp(), file.c_str()));
}

void HighsBackend::DoWriteSolution(const std::string &file) {
  HIGHS_CCALL(loader().Highs_writeSolutionPretty(lp(), file.c_str()));
}


void HighsBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptHighs, lp());
}


void HighsBackend::SetupPlateauCallbacks() {
  if (!loader().Highs_setCallback || !loader().Highs_startCallback)
    throw std::runtime_error(
        "This HiGHS build does not provide the callback API "
        "(Highs_setCallback); mip:plateautime is not available.");
  HIGHS_CCALL( loader().Highs_setCallback(
      lp(), &HighsBackend::DoPlateauCallback, this) );
  HIGHS_CCALL( loader().Highs_startCallback(
      lp(), kHighsCallbackMipImprovingSolution) );
    HIGHS_CCALL( loader().Highs_startCallback(
        lp(), kHighsCallbackMipInterrupt) );
    HIGHS_CCALL(loader().Highs_startCallback(
        lp(), kHighsCallbackMipSolution));
    }

void HighsBackend::DoPlateauCallback(
    int callback_type, const char* /*message*/,
    const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
    void* user_callback_data) {
  auto* backend = static_cast<HighsBackend*>(user_callback_data);

  // Check that we're in a correct callback. Probably non necessary.
  if ((callback_type == kHighsCallbackMipImprovingSolution) ||
      (callback_type == kHighsCallbackMipInterrupt) ||
      (callback_type == kHighsCallbackMipSolution))
  {
      double relgap = data_out->mip_gap;
      double absgap = std::fabs(data_out->mip_primal_bound - data_out->mip_dual_bound);
      if (callback_type == kHighsCallbackMipImprovingSolution) {
          // If new solution is reported, send it to the plateau monitor
          if (backend->ReportIncumbentForPlateau(data_out->objective_function_value))
          {
              data_in->user_interrupt = 1;
              return;
          }
      }
	  // Now check if we should stop due to plateau due to mip gap changes
      // Eventual solution improvement is already handled
        if (backend->ReportGapForPlateau(absgap, relgap))
          data_in->user_interrupt = 1;
  } 
}


void HighsBackend::Solve() {
  int status = loader().Highs_run(lp());
  if (status != kHighsStatusOk && status != kHighsStatusWarning)
    throw std::runtime_error(fmt::format("  Error {} while solving with HiGHS").c_str());
    // Mask warnings for pdlp
  if((status == kHighsStatusWarning) && (storedOptions_.lpmethod_!="pdlp-gpu"))
      Print("  Warning code {} while solving with HiGHS\n", status);
  WindupHIGHSSolve();
}

HighsBackend::FixedLPResult HighsBackend::SolveFixedLP(ArrayRef<double> x) {
  if (!loader().Highs_getLp || !loader().Highs_passLp || !loader().Highs_getNumNz)
    MP_RAISE("HiGHS library lacks fixed-model LP API");
  int feasible = 0;
  HIGHS_CCALL(loader().Highs_getIntInfoValue(lp(), "primal_solution_status", &feasible));
  if (feasible != kHighsSolutionStatusFeasible)
    MP_RAISE("no integer-feasible incumbent");
  if (accObjectives().hadNativeMultiObj())
    MP_RAISE("native multiple objectives are unsupported by fixed-model emulation");
  int n = NumVars(), m = NumLinCons(), nz = loader().Highs_getNumNz(lp());
  int sense, nc, nr, nnz;
  double offset, tol;
  GetSolverOption("mip_feasibility_tolerance", tol);
  std::vector<double> cost(n), lb(n), ub(n), lhs(m), rhs(m), val(nz);
  std::vector<int> start(n+1), index(nz), integrality(n);
  HIGHS_CCALL(loader().Highs_getLp(lp(), kHighsMatrixFormatColwise,
      &nc, &nr, &nnz, &sense, &offset, cost.data(), lb.data(), ub.data(),
      lhs.data(), rhs.data(), start.data(), index.data(), val.data(), integrality.data()));
  for (int j = 0; j < n; ++j) {
    if (integrality[j] == kHighsVarTypeInteger)
      lb[j] = ub[j] = FixedIntegerValue(x[j], lb[j], ub[j], tol);
    else if (integrality[j] != kHighsVarTypeContinuous)
      MP_RAISE("unsupported semi-continuous/semi-integer variable");
  }
  auto destroy = [this](void* model) { loader().Highs_destroy(model); };
  std::unique_ptr<void, decltype(destroy)> fixed(loader().Highs_create(), destroy);
  if (!fixed) MP_RAISE("cannot create fixed LP");
  int output = 0, threads = 0;
  GetSolverOption("output_flag", output);
  GetSolverOption("threads", threads);
  HIGHS_CCALL(loader().Highs_setBoolOptionValue(fixed.get(), "output_flag", output));
  HIGHS_CCALL(loader().Highs_setIntOptionValue(fixed.get(), "threads", threads));
  HIGHS_CCALL(loader().Highs_setStringOptionValue(fixed.get(), "solver", "simplex"));
  if (fixed_model_time_limit() >= 0)
    HIGHS_CCALL(loader().Highs_setDoubleOptionValue(fixed.get(), "time_limit", fixed_model_time_limit()));
  if (fixed_model_iter_limit() >= 0)
    HIGHS_CCALL(loader().Highs_setIntOptionValue(fixed.get(), "simplex_iteration_limit", fixed_model_iter_limit()));
  HIGHS_CCALL(loader().Highs_passLp(fixed.get(), n, m, nz,
      kHighsMatrixFormatColwise, sense, offset, cost.data(), lb.data(), ub.data(),
      lhs.data(), rhs.data(), start.data(), index.data(), val.data()));
  HIGHS_CCALL(loader().Highs_run(fixed.get()));
  int status = loader().Highs_getModelStatus(fixed.get());
  if (status != kHighsModelStatusOptimal)
    MP_RAISE(fmt::format("fixed LP not optimal (HiGHS status {})", status));
  FixedLPResult result;
  result.objective = loader().Highs_getObjectiveValue(fixed.get());
  result.primal.resize(n);
  result.dual.resize(m);
  result.varstt.resize(n);
  result.constt.resize(m);
  HIGHS_CCALL(loader().Highs_getSolution(fixed.get(), result.primal.data(), nullptr, nullptr, result.dual.data()));
  HIGHS_CCALL(loader().Highs_getBasis(fixed.get(), result.varstt.data(), result.constt.data()));
  auto convert = [](int status) {
    switch (status) {
    case kHighsBasisStatusBasic: return (int)BasicStatus::bas;
    case kHighsBasisStatusLower: return (int)BasicStatus::low;
    case kHighsBasisStatusUpper: return (int)BasicStatus::upp;
    case kHighsBasisStatusNonbasic:
    case kHighsBasisStatusZero: return (int)BasicStatus::sup;
    default: MP_RAISE(fmt::format("Unknown HiGHS fixed-LP basis status {}", (int)status));
    }
  };
  for (auto& status : result.varstt) status = convert(status);
  for (auto& status : result.constt) status = convert(status);
  HIGHS_CCALL(loader().Highs_getIntInfoValue(fixed.get(), "simplex_iteration_count", &result.iterations));
  return result;
}

void HighsBackend::WindupHIGHSSolve() { }

void HighsBackend::ReportResults() {
   ReportHIGHSResults();
  ConsiderEmulatedFixedModel();
  BaseBackend::ReportResults();
}

void HighsBackend::ReportHIGHSResults() {
  SetStatus( GetSolveResult() );
  AddHIGHSMessages();
}

SolutionBasis HighsBackend::GetBasis() {
  if (IsMIP()) return FixedLPBasis(CG_Linear);
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

void HighsBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarConStatii(varstt, constt);
}


/// We also use this for MIP.
/// Attempting partial MIP start.
void HighsBackend::AddMIPStart(ArrayRef<double> primal, ArrayRef<int> sparsity) {
  Solution start;
  start.primal = primal;
  start.spars_primal = sparsity;
  AddPrimalDualStart(std::move(start));
}

void HighsBackend::AddPrimalDualStart(Solution sol0_unpres) {
  auto mv = GetValuePresolver().PresolveSolution(
    { sol0_unpres.primal, sol0_unpres.dual });
  auto ms = GetValuePresolver().PresolveGenericInt(
        { sol0_unpres.spars_primal } );
  auto x0 = mv.GetVarValues()();
  auto s0 = ms.GetVarValues()();
  auto pi0 = mv.GetConValues()(CG_Linear);
  if (s0.size() < (size_t)NumVars())
    s0.resize(NumVars());
  if (x0.size() < (size_t)NumVars())
    x0.resize(NumVars());
  /// Partial start: if all missing values belong to continuous variables,
  /// pass the start as is; otherwise either let HiGHS complete it
  /// (mip:start=1, mip:maxstartnodes) or complete it here (mip:start=2).
  if (0==*std::min_element(s0.begin(), s0.end())) {
    bool fAllMissingAreRealVars = true;
    for (auto j=s0.size(); j--; ) {
      if (s0[j])
        continue;
      int integr;
      auto res = loader().Highs_getColIntegrality(lp(), j, &integr);
      if (kHighsStatusOk != res)
        break;       // no information, it's an LP
      if (kHighsVarTypeContinuous != integr) {
        fAllMissingAreRealVars = false;
        break;
      }
    }
    if (!fAllMissingAreRealVars) {
      if (2 == storedOptions_.mipstart_ || !loader().Highs_setSparseSolution) {
        if (!CompleteStartByFixing(x0, s0))
          return;
      } else {
        std::vector<int> idx;
        std::vector<double> val;
        for (size_t j = 0; j < s0.size(); ++j)
          if (s0[j]) {
            idx.push_back((int)j);
            val.push_back(x0[j]);
          }
        auto res = loader().Highs_setSparseSolution(
              lp(), (int)idx.size(), idx.data(), val.data());
        if (kHighsStatusError == res)
          AddWarning("MIPSTART",
                     "HiGHS rejected the partial MIP start "
                     "(a value may be outside the variable bounds)");
        return;
      }
    }
  }
  auto res = loader().Highs_setSolution(lp(), x0.data(), NULL, NULL,
                                        pi0.empty() ? NULL : pi0.data());
  if (kHighsStatusError == res)
    AddWarning("WARMSTART",
               "HiGHS rejected the starting point "
               "(a value may be outside the variable bounds)");
}

bool HighsBackend::CompleteStartByFixing(std::vector<double>& x0,
                                         const std::vector<int>& s0) {
  std::vector<double> costs(NumVars());
  std::vector<double> lb(NumVars());
  std::vector<double> ub(NumVars());
  int numnz, ncols;
  HIGHS_CCALL(loader().Highs_getColsByRange(lp(), 0, NumVars()-1, &ncols,
    costs.data(), lb.data(), ub.data(), &numnz, NULL, NULL, NULL));
  HIGHS_CCALL(loader().Highs_changeColsBoundsByMask(
                lp(), s0.data(), x0.data(), x0.data()));
  bool ok = false;
  try {
    auto run = loader().Highs_run(lp());
    int primal_status = kHighsSolutionStatusNone;
    loader().Highs_getIntInfoValue(lp(), "primal_solution_status", &primal_status);
    if (kHighsStatusError != run && kHighsSolutionStatusFeasible == primal_status) {
      std::vector<double> x(NumVars());
      if (kHighsStatusError !=
          loader().Highs_getSolution(lp(), x.data(), NULL, NULL, NULL)) {
        x0 = std::move(x);
        ok = true;
      }
    }
  } catch (...) {
    HIGHS_CCALL(loader().Highs_changeColsBoundsByMask(
                  lp(), s0.data(), lb.data(), ub.data()));
    throw;
  }
  HIGHS_CCALL(loader().Highs_changeColsBoundsByMask(
                lp(), s0.data(), lb.data(), ub.data()));
  if (!ok)
    AddWarning("MIPSTART",
               "Could not complete the partial MIP start by fixing "
               "the provided integer values; start ignored");
  return ok;
}


ArrayRef<int> HighsBackend::VarStatii() {
  std::vector<int> vars(NumVars());
  conStatiii_.resize(NumLinCons());
  HIGHS_CCALL(loader().Highs_getBasis(lp(), vars.data(), conStatiii_.data()));
  for (auto& s : vars) {
    switch (s) {
    case kHighsBasisStatusBasic:
      s = (int)BasicStatus::bas;
      break;
    case kHighsBasisStatusLower:
      s = (int)BasicStatus::low;
      break;
    case kHighsBasisStatusUpper:
      s = (int)BasicStatus::upp;
      break;
    case kHighsBasisStatusNonbasic:
    case kHighsBasisStatusZero:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Highs VBasis value: {}", s));
    }
  }
  return vars;
}

ArrayRef<int> HighsBackend::ConStatii() {
  for (auto& s : conStatiii_) {
    switch (s) {
    case kHighsBasisStatusBasic:
      s = (int)BasicStatus::bas;
      break;
    case kHighsBasisStatusLower:
      s = (int)BasicStatus::low;
      break;
    case kHighsBasisStatusUpper:
      s = (int)BasicStatus::upp;
      break;
    case kHighsBasisStatusNonbasic:
      s = (int)BasicStatus::sup;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Highs VBasis value: {}", s));
    }
  }
  return conStatiii_;
}

void HighsBackend::VarConStatii(ArrayRef<int> vst, ArrayRef<int> cst) {
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  std::vector<int> indicesOfMissing;
  for (size_t j = 0; j<stt.size(); j++) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = kHighsBasisStatusBasic;
      break;
    case BasicStatus::low:
    case BasicStatus::equ:
      s = kHighsBasisStatusLower;
      break;
    case BasicStatus::upp:
      s = kHighsBasisStatusUpper;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = kHighsBasisStatusNonbasic;
      break;
    case BasicStatus::none:
      indicesOfMissing.push_back(j);
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }

  if (indicesOfMissing.size() > 0)
  {
    /// 'none' is assigned to new variables. Compute low/upp/sup:
    /// Depending on where 0.0 is between bounds
    std::vector<double> lb(indicesOfMissing.size());
    std::vector<double> ub(indicesOfMissing.size());
    std::vector<int> di(indicesOfMissing.size(),   // dummy pointers
                        indicesOfMissing.size());
    std::vector<double> dd(indicesOfMissing.size());
    int numnz;
    loader().Highs_getColsBySet(lp(), indicesOfMissing.size(), indicesOfMissing.data(),
      di.data(), dd.data(), lb.data(), ub.data(), &numnz, NULL, NULL, NULL);
    for (size_t i = 0; i < indicesOfMissing.size(); i++) {
        if (lb[i] >= -1e-6)
          stt[indicesOfMissing[i]] = kHighsBasisStatusLower;
        else if (ub[i] <= 1e-6)
          stt[indicesOfMissing[i]] = kHighsBasisStatusUpper;
        else
          stt[indicesOfMissing[i]] = kHighsBasisStatusNonbasic;
      }
  }
  std::vector<int> cstt(cst.data(), cst.data() + cst.size());
  for (auto& s : cstt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = kHighsBasisStatusBasic;
      break;
    case BasicStatus::none:
    case BasicStatus::upp:   
    case BasicStatus::sup:   
    case BasicStatus::low:    
    case BasicStatus::equ:    
    case BasicStatus::btw:    
      s = kHighsBasisStatusNonbasic;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  HIGHS_CCALL(loader().Highs_setBasis(lp(), stt.data(), cstt.data()));
}

ArrayRef<double> HighsBackend::Ray() {
  HighsInt has_ray;
  std::vector<double> uray_pres(NumVars());
  auto res = loader().Highs_getPrimalRay(lp(), &has_ray, uray_pres.data());
  if (res)
    Print("Error while getting primal ray");
  if (res || (!has_ray))
  {
    uray_pres.clear();
    return uray_pres;
  }
  auto mv = GetValuePresolver().PostsolveSolution({ uray_pres });
  auto uray = mv.GetVarValues()();
  return uray;
}

ArrayRef<double> HighsBackend::DRay() {
  HighsInt has_ray;
  std::vector<double> dray_pres(NumLinCons());
  auto res = loader().Highs_getDualRay(lp(), &has_ray, dray_pres.data());
  if (res)
    Print("Error while getting dual ray");
  if (res || (!has_ray))
  {
    dray_pres.clear();
    return dray_pres;
  }
  auto vm = GetValuePresolver().PostsolveSolution({
                                               {},
                                               {{{CG_Linear, std::move(dray_pres)}}}
    });
  return vm.GetConValues().MoveOut();        // need the vector itself
}



void HighsBackend::AddHIGHSMessages() {
  auto pdlp = PdlpIterations();
  if (pdlp > 0)
  {
    AddToSolverMessage(fmt::format("{} PDLP iterations\n", pdlp));
    return;
  }
  auto ni = SimplexIterations();
  if (true)
    AddToSolverMessage(
          fmt::format("{} simplex iterations\n", std::max(0, ni)));
  auto nbi = BarrierIterations();
  if (nbi > -1)
    AddToSolverMessage(
      fmt::format("{} barrier iterations\n", nbi));
  auto nnd = NodeCount();
  if (nnd > -1)
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}


std::map<std::string, std::variant<int, double, std::string>>
HighsBackend::SolutionStats() {
    std::map<std::string, std::variant<int, double, std::string>> stats;
    stats["simplex_iterations"] = SimplexIterations();
    stats["barrier_iterations"] = BarrierIterations();
    stats["node_count"] = NodeCount();
    stats["pdlp_iterations"] = PdlpIterations();
    return stats;
}

std::pair<int, std::string> HighsBackend::GetSolveResult() {
  namespace sol = mp::sol;
  int optstatus = loader().Highs_getModelStatus(lp());
  auto obj = loader().Highs_getObjectiveValue(lp());
  auto inf = loader().Highs_getInfinity(lp());
  bool hasSol = (-inf < obj && obj < inf);
  int primal_solution_status;
  loader().Highs_getIntInfoValue(lp(),
                        "primal_solution_status", &primal_solution_status);
  switch (optstatus) {
  case kHighsModelStatusOptimal:
    return { sol::SOLVED, "optimal solution" };
  case kHighsModelStatusInfeasible:
    return { sol::INFEASIBLE, "infeasible problem" };
  case kHighsModelStatusUnbounded:
    if (hasSol)
      return { sol::UNBOUNDED_FEAS, "unbounded problem, feasible solution" };
    return { sol::UNBOUNDED_NO_FEAS, "unbounded problem, no solution" };
  case kHighsModelStatusUnboundedOrInfeasible:
    return { sol::LIMIT_INF_UNB, "unbounded or infeasible. "
                                "Disable presolve (pre:solve=off) "
                                "or run the IIS finder (iis:find=1) "
                                "for a definitive answer." };
  case kHighsModelStatusModelEmpty:
    return { sol::SOLVED, "empty model" };
  case kHighsModelStatusModelError:
  case kHighsModelStatusLoadError:
    if (kHighsSolutionStatusInfeasible == primal_solution_status)   // HiGHS 7
      return { sol::INFEASIBLE, "infeasible problem" };
    return { sol::FAILURE, "solver error" };
  case kHighsModelStatusPresolveError:
  case kHighsModelStatusSolveError:
  case kHighsModelStatusPostsolveError:
    if (hasSol)
      return { sol::UNCERTAIN, "numeric issue, solution candidate returned" };
    return { sol::NUMERIC, "numeric issue" };
  case kHighsModelStatusTimeLimit:
    if (hasSol)
      return { sol::LIMIT_FEAS_TIME, "time limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_TIME, "time limit, no solution" };
  case kHighsModelStatusIterationLimit:
    if (hasSol)
      return { sol::LIMIT_FEAS_ITER, "iteration limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS_ITER, "iteration limit, no solution" };
  case kHighsModelStatusSolutionLimit:
    if (hasSol) return { sol::LIMIT_FEAS_NUMSOLS, "solution limit, feasible solution" };
    return { sol::LIMIT_NO_FEAS, "solution limit, no solution" };
  case kHighsModelStatusInterrupt:
    if (hasSol)
      return { sol::LIMIT_FEAS_INTERRUPT, "interrupt, feasible solution" };
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupt, no solution" };
  case kHighsModelStatusObjectiveBound:
    if (hasSol)
      return { sol::LIMIT_FEAS_BESTOBJ, "objective bound, feasible solution" };
    return { sol::LIMIT_NO_FEAS_BESTBND, "objective bound, no solution" };
  case kHighsModelStatusObjectiveTarget:
    if (hasSol)
      return { sol::LIMIT_FEAS_BESTOBJ, "objective target, feasible solution" };
    return { sol::LIMIT_NO_FEAS_CUTOFF, "objective target, no solution" };
  default:
    if (hasSol)
      return { sol::UNCERTAIN, "unknown, solution candidate returned" };
    return { sol::UNKNOWN, "unknown" };
  }
  return { sol::UNKNOWN, "not solved" };
}



////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo lp_values_method[] = {
  { "choose", "Automatic (default)", -1},
  { "simplex", "Simplex", 1},
  { "ipm", "Interior Point Method (automatic)", 2},
  { "hipo", "HiPO interior point method (not available in every HiGHS library build; "
            "rejected with a message if unavailable)", 2},
  { "ipx", "IPX interior point method", 2},
  { "pdlp", "cuPDLP-c solver", 3},
  { "pdlp-gpu", "cuPDLP-c solver on NVIDIA GPU. Requires CUDA, not available on MacOS", 4},
  { "hipdlp", "HiPDLP solver", 5},
  { "hipdlp-gpu", "HiPDLP-c solver on NVIDIA GPU. Requires CUDA, not available on MacOS", 6},
  { "qpasm", "QP active set method", 7}
};


static const mp::OptionValueInfo mip_lp_solver_values[] = {
  { "choose", "Automatic (default)", -1},
  { "simplex", "Simplex", 1},
  { "ipm", "Interior Point Method (automatic)", 2},
  { "hipo", "HiPO interior point method", 3},
  { "ipx", "IPX interior point method", 4 }
};

static const mp::OptionValueInfo mip_ipm_solver_values[] = {
  { "choose", "Automatic (default)", -1},
  { "ipx", "IPX interior point method", 2},
  { "hipo", "HiPO interior point method", 3}
};

static const mp::OptionValueInfo off_on_choose_values[] = {
  { "choose", "Automatic (default)", -1},
  { "off", "Off", 1},
  { "on", "On", 2},
};


static const mp::OptionValueInfo pdlperestartmethod_values[] = {
  { "0", "None", 0},
  { "1", "GPU (default)", 1},
  { "2", "CPU", 2},
};

static const mp::OptionValueInfo pdlp_scaling_mode_values[] = {
  { "1", "Ruiz scaling", 0},
  { "2", "L2", 1},
  { "4", "PC (preconditioner-oriented)", 2},
};

static const mp::OptionValueInfo pdlp_step_size_strategy_values[] = {
    {"0", "Fixed", 0},
    {"1", "Adaptive (default)", 1},
    {"2", "Malitsky-Pock", 2},
    {"3", "PID", 3}
};

static const mp::OptionValueInfo pdlp_restart_strategy_values[] = {
  { "0", "Off", 0},
  { "1", "Fixed", 1},
  { "2", "Adaptive (default)", 2},
  {"3", "Halpern", 3}
};

static const mp::OptionValueInfo mipstart_values[] = {
  { "1", "Pass the provided values to HiGHS, which completes a partial "
         "start itself by a limited MIP search (see mip:maxstartnodes; default)", 1},
  { "2", "Complete a partial start in the driver: fix the provided values "
         "of integer variables, solve, unfix, and pass the complete solution "
         "(legacy behavior)", 2}
};

static const mp::OptionValueInfo iis_strategy_values[] = {
    {"0", "Light test only (infeasible bounds)", 0 },
    { "1", "Try dual ray", 1 },
    { "2", "Try elastic LP", 2 },
    { "4", "Find an irreducible (true) IIS by a deletion filter", 4 },
    { "8", "Prioritise columns over rows in the deletion filter", 8 },
    { "16", "Find the IIS of the LP relaxation for a MIP", 16 }
};

static const mp::OptionValueInfo run_crossover_values[] = {
  { "choose", "Run if the results of IPM without crossover is imprecise", -1},
  { "off", "Off", 1},
  { "on", "On (default)", 2},
};

static const mp::OptionValueInfo hipo_system_values[] = {
    { "choose", "Choose automatically (default)", 0},
    { "augmented", "Augmented system", 1},
    { "normaleq", "Normal equations", 2}
};
static const mp::OptionValueInfo hipo_parallel_type_values[] = {
    // tree, node or both
	{ "tree", "Tree parallelism", 1},
	{ "node", "Node parallelism", 2},
	{ "both", "Both tree and node parallelism", 3}
};
static const mp::OptionValueInfo hipo_ordering_values[] = {
	{ "choose", "Choose automatically (default)", 0},
	{ "metis", "Multilevel graph partitioning", 1},
	{ "amd", "Approximate minimum degree", 2},
	{ "rcm", "Reverse Cuthill-McKee", 3}
};
static const mp::OptionValueInfo simplex_strategy_values_[] = {
  { "0", "Choose automatically", 0},
  { "1", "Dual, serial (default)", 1},
  { "2", "Dual, 'SIP' (Single Iteration Parallelism)", 2},
  { "3", "Dual, 'PAMI' (Parallelization Across Multiple Iterations)", 3},
  { "4", "Primal", 4}
};
static const mp::OptionValueInfo simplex_scale_strategy_values_[] = {
  { "0", "Off", 0},
  { "1", "Choose automatically", 1},
  { "2", "Equilibration (default)", 2},
  { "3", "Forced equilibration", 3},
  { "4", "Max value", 4}
};
static const mp::OptionValueInfo simplex_crash_strategy_values_[] = {
  { "0", "Off (default)", 0},
  { "1", "LTSSF", 1},
  { "2", "Bixby", 2}
};
/// HiGHS accepts further experimental crash strategies (3..9), undocumented.
static const mp::OptionValueInfo simplex_edge_weight_strategy_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Dantzig", 0},
  { "1", "Devex", 1},
  { "2", "Steepest", 2}
};

void HighsBackend::InitCustomOptions() {

  set_option_header(
    "HIGHS Optimizer Options for AMPL\n"
    "--------------------------------------------\n"
    "\n"
    "To set these options, assign a string specifying their values to the "
    "AMPL option ``highs_options``. For example::\n"
    "\n"
    "  ampl: option highs_options 'mip:gap=1e-6';\n");

  AddHighsOption("tech:outlev outlev",
    "0*/1: Whether to write HiGHS log lines (chatter) to stdout and to file.",
    "output_flag", 0, 1);

  AddHighsOption("tech:logfile logfile",
    "Log file name.", "log_file");

  std::string c;
  AddStoredOption("alg:method method lp:method lpmethod lp:algorithm solver",
    "Which algorithm to use for non-MIP problems:\n"
    "\n.. value-table::\n", storedOptions_.lpmethod_, lp_values_method);

  AddHighsOption("alg:simplex simplex simplex_strategy",
    "Strategy for simplex solver :\n"
    "\n.. value-table::\n", "simplex_strategy", simplex_strategy_values_, 0);

  AddHighsOption("alg:simplexscale simplexscale simplex_scale_strategy",
    "Simplex scaling strategy :\n"
    "\n.. value-table::\n", "simplex_scale_strategy",
    simplex_scale_strategy_values_, 1);

  AddHighsOption("alg:simplexcrash simplexcrash simplex_crash_strategy",
    "Simplex crash strategy :\n"
    "\n.. value-table::\n", "simplex_crash_strategy",
    simplex_crash_strategy_values_, 0);

  AddHighsOption("alg:simplexdualedge simplexdualedge simplex_dual_edge_weight_strategy",
    "Simplex dual edge weights strategy :\n"
    "\n.. value-table::\n", "simplex_dual_edge_weight_strategy",
    simplex_edge_weight_strategy_values_, 1);

  AddHighsOption("alg:simplexprimaledge simplexprimaledge simplex_primal_edge_weight_strategy",
    "Simplex primal edge weights strategy :\n"
    "\n.. value-table::\n", "simplex_primal_edge_weight_strategy",
    simplex_edge_weight_strategy_values_, 1);

  AddHighsOption("alg:simplexupdatelim simplexupdatelim simplex_update_limit",
    "Limit on the number of simplex UPDATE operations before "
    "refactorization (default 5000).",
    "simplex_update_limit", 0, INT_MAX);

  AddHighsOption("alg:simplexminconcurrency simplexminconcurrency simplex_min_concurrency",
    "Minimum level of concurrency in parallel simplex (1..8, default 1).",
    "simplex_min_concurrency", 1, 8);

  AddHighsOption("alg:simplexmaxconcurrency simplexmaxconcurrency simplex_max_concurrency",
    "Maximum level of concurrency in parallel simplex (1..8, default 8).",
    "simplex_max_concurrency", 1, 8);

  AddHighsOption("pre:solve presolve",
    "Whether to use presolve:\n"
    "\n.. value-table::\n",
    "presolve", off_on_choose_values, c);



  AddHighsOption("pre:userobjectivescale user_objective_scale userobjectivescale",
      "Exponent of power-of-two objective scaling for model (default 0; "
      "negative exponents allowed).",
      "user_objective_scale", std::numeric_limits<HighsInt>::min(),
      std::numeric_limits<HighsInt>::max());


  AddHighsOption("pre:userboundscale user_bound_scale userboundscale",
    "Exponent of power-of-two bound scaling for model (default 0; "
    "negative exponents allowed).",
    "user_bound_scale", std::numeric_limits<HighsInt>::min(),
    std::numeric_limits<HighsInt>::max());

  AddHighsOption("alg:parallel parallel",
    "Parallel option :\n"
    "\n.. value-table::\n", "parallel", off_on_choose_values, c);

  AddHighsOption("lim:time timelim timelimit time_limit",
    "Limit on solve time (in seconds; default: no limit).",
    "time_limit", 0.0, DBL_MAX);

  AddHighsOption("lim:iter iterlim lim:simplexiterationlimit simplexiterationlimit simplex_iteration_limit",
    "Limit on simplex iterations when solving LPs, not counting "
    "MIP subproblems (default: no limit).",
    "simplex_iteration_limit", 0, INT_MAX);

  AddHighsOption("bar:iterlim bariterlim lim:ipmiterationlimit ipmiterationlimit ipm_iteration_limit",
    "Limit on IPM iterations (default: no limit).",
    "ipm_iteration_limit", 0, INT_MAX);

  AddHighsOption("lim:objectivebound objective_bound objectivebound",
    "Objective bound for termination of the dual simplex solver (default: no limit).",
    "objective_bound", -DBL_MAX, DBL_MAX);


  AddHighsOption("lim:objectivetarget objective_target objectivetarget",
    "Objective target for termination of the MIP solver (default: no limit).",
    "objective_target", -DBL_MAX, DBL_MAX);

  /// Not in Highs 1.14:
  // AddHighsOption("lim:pdlpnativetermination pdlp_native_termination pdlpnativetermination",
  //   "Use native termination for PDLP solver:\n"
  //   "\n.. value-table::\n", "pdlp_native_termination", values_01_noyes_0default_, 0);

  AddHighsOption("pdlp:scaling pdlpscaling pre:pdlpscaling_mode pdlp_scaling_mode",
      "Scaling mode for PDLP solver (default 5 = 1 + 4), sum of:\n"
      "\n.. value-table::\n", "pdlp_scaling_mode", pdlp_scaling_mode_values, 5);

  AddHighsOption("lim:pdlpiter pdlpiterlim lim:pdlp_iteration_limit pdlpiterationlimit pdlp_iteration_limit",
    "Iteration limit for PDLP solver (default: no limit).",
    "pdlp_iteration_limit", 0, INT_MAX);

  AddHighsOption("alg:infinitecost infinitecost infinite_cost",
    "Limit on |cost coefficient|: values greater than or equal to this "
    "will be treated as infinite (default: 1e20).",
    "infinite_cost", 1e15, Infinity());

  AddHighsOption("alg:infinitebound infinitebound infinite_bound",
    "Limit on |constraint bound|: values greater than or equal to this "
    "will be treated as infinite (default: 1e20).",
    "infinite_bound", 1e15, Infinity());

  AddHighsOption("alg:infinitecoeff infinitecoeff large_matrix_value",
    "Upper limit on |matrix entries|: values greater than or equal to this "
    "will be treated as infinite (default: 1e15).",
    "large_matrix_value", 1.0, Infinity());

  AddHighsOption("alg:zerocoeff zerocoeff small_matrix_value",
    "Lower limit on |matrix entries|: values less than or equal to this "
    "will be treated as zero (default: 1e-9).",
    "small_matrix_value", 1e-12, Infinity());

  AddHighsOption("alg:feastol feastol primal_feasibility_tolerance",
      "Primal feasibility tolerance (default 1e-7).",
      "primal_feasibility_tolerance", 1e-10, Infinity());

  AddHighsOption("alg:kkttol kkttol alg:kkt_tolerance kkt_tolerance",
    "KKT tolerance (default 1e-7). If changed from its default, "
    "this tolerance is used for all feasibility and optimality (KKT) "
    "measures, overriding the individual tolerances.",
    "kkt_tolerance", 1e-10, Infinity());

  AddHighsOption("alg:dualfeastol dualfeastol dual_feasibility_tolerance",
    "Dual feasibility tolerance (default 1e-7).",
    "dual_feasibility_tolerance", 1e-10, Infinity());

  AddHighsOption("alg:primalrestol primalrestol primal_residual_tolerance",
    "Primal residual tolerance (default 1e-7).",
    "primal_residual_tolerance", 1e-10, Infinity());

  AddHighsOption("alg:opttol opttol optimality_tolerance",
    "Optimality tolerance (default 1e-7).",
    "optimality_tolerance", 1e-10, Infinity());

  AddHighsOption("alg:dualrestol dualrestol dual_residual_tolerance",
    "Dual residual tolerance (default 1e-7).",
    "dual_residual_tolerance", 1e-10, Infinity());

  AddHighsOption("alg:ipmopttol ipmopttol ipm_optimality_tolerance",
    "IPM optimality tolerance (default 1e-8).",
    "ipm_optimality_tolerance", 1e-12, Infinity());

  AddHighsOption("pdlp:restartmethod pdlprestartmethod pdlp:cupdlpc_restart_method pdlpcupdlpcrestartmethod pdlp_cupdlpc_restart_method",
    "Restart mode for PDLP solver:\n"
      "\n.. value-table::\n",
    "pdlp_cupdlpc_restart_method", pdlperestartmethod_values, 1);

  AddHighsOption("pdlp:restart pdlprestart pdlp:restart_strategy pdlprestartstrategy pdlp_restart_strategy",
      "Restart strategy for PDLP solver:\n"
      "\n.. value-table::\n",
      "pdlp_restart_strategy", pdlp_restart_strategy_values, 1);

  AddHighsOption("pdlp:stepsize pdlpstepsize pdlp:step_size_strategy pdlpstepsizestrategy pdlp_step_size_strategy",
      "Stepsize strategy for PDLP solver:\n"
      "\n.. value-table::\n",
      "pdlp_step_size_strategy", pdlp_step_size_strategy_values, 1);

  AddHighsOption("pdlp:opttol pdlpopttol pdlp:optimality_tolerance pdlp_optimality_tolerance pdlp_opt_tol",
    "PDLP optimality tolerance (default 1e-7).",
    "pdlp_optimality_tolerance", 1e-12, Infinity());


  AddHighsOption("pdlp:ruiziter pdlpruiziter pdlp:ruiz_iterations pdlp_ruiz_iterations",
      "Number of Ruiz scaling iterations for PDLP solver (default 10).",
      "pdlp_ruiz_iterations", 0, INT_MAX);

  AddHighsOption("lim:qpiter qpiterlim lim:qp_iteration_limit qp_iteration_limit qpiterationlimit",
      "Iteration limit for QP active set method (default: no limit).",
      "qp_iteration_limit", 0, INT_MAX);

  AddHighsOption("qp:nullspacelim qpnullspacelim qp:nullspace_limit qp_nullspace_limit qpnullspacelimit",
      "Nullspace limit for the active set QP solver (default: 4000).",
      "qp_nullspace_limit", 0, INT_MAX);

  AddHighsOption("qp:regularization qpregularization qp:regularization_value qp_regularization_value qpregularizationvalue",
      "Regularization value added to the Hessian in the active set QP solver (default: 1e-7).",
      "qp_regularization_value", 0.0, DBL_MAX);

  AddHighsOption("bar:crossover crossover run_crossover",
    "Run crossover after IPM to get a basic solution:\n"
      "\n.. value-table::\n", "run_crossover", run_crossover_values, c);

  AddHighsOption("bar:hiposystem hiposystem bar:hipo_system hipo_system",
      "Control which Newton/KKT system is used by HiPO:\n"
      "\n.. value-table::\n",  "hipo_system", hipo_system_values , c);


  AddHighsOption("bar:hipoparallel hipoparallel bar:hipo_parallel_type hipo_parallel_type",
	  "Control the type of parallelism used by HiPO:\n"
      "\n.. value-table::\n", "hipo_parallel_type", hipo_parallel_type_values, c);

  AddHighsOption("bar:hipoordering hipoordering bar:hipo_ordering hipo_ordering",
	  "Control the ordering used by HiPO before factorization:\n"
      "\n.. value-table::\n", "hipo_ordering", hipo_ordering_values, c);

  AddHighsOption("bar:hipoblocksize hipoblocksize hipo_block_size",
    "Block size for dense linear algebra within HiPO (default 128).",
    "hipo_block_size", 0, INT_MAX);

  AddHighsOption("qp:hotstart qphotstart qp_allow_hot_start",
    "0*/1: Allow the active set QP solver to hot start.",
    "qp_allow_hot_start", 0, 1);

  AddHighsOption("tech:threads threads",
    "Maximum number of threads used by HiGHS (parallel simplex, "
    "IPM and MIP); default 0 ==> automatic choice.",
    "threads", 0, INT32_MAX);
  AddHighsOption("mip:symmetry symmetry mip:detsimmetry detsimmetry mip_detect_symmetry",
    "0/1*: Whether MIP symmetry should be detected.",
    "mip_detect_symmetry", 0, 1);

  AddHighsOption("mip:restart miprestart mip_allow_restart",
    "0/1*: Whether MIP restart is permitted.",
    "mip_allow_restart", 0, 1);

  AddStoredOption("mip:start mipstart",
    "How a partial MIP start (initial values for some integer variables, "
    "with alg:start>0) is completed:\n"
    "\n.. value-table::\n", storedOptions_.mipstart_, mipstart_values);

  AddHighsOption("mip:maxstartnodes maxstartnodes mip_max_start_nodes",
    "Maximum number of MIP nodes when completing a partial MIP start "
    "(default 500). See mip:start.",
    "mip_max_start_nodes", 0, INT_MAX);

  AddHighsOption("mip:poolagelimit poolagelimit mip_pool_age_limit",
    "Maximal age of rows in the MIP cutpool before they are deleted "
    "(default 30, max 1000).",
    "mip_pool_age_limit", 0, 1000);

  AddHighsOption("mip:nodecuts nodecuts mip_allow_cut_separation_at_nodes",
    "0/1*: Whether cut separation at nodes other than the root is permitted.",
    "mip_allow_cut_separation_at_nodes", 0, 1);

  AddHighsOption("mip:rootpresolveonly rootpresolveonly mip_root_presolve_only",
    "0*/1: Whether MIP presolve is only applied at the root node.",
    "mip_root_presolve_only", 0, 1);

  AddHighsOption("mip:heur:feasjump heurfeasjump mip_heuristic_run_feasibility_jump",
    "0/1*: Use the feasibility jump heuristic.",
    "mip_heuristic_run_feasibility_jump", 0, 1);

  AddHighsOption("mip:heur:rins heurrins mip_heuristic_run_rins",
    "0/1*: Use the RINS heuristic.",
    "mip_heuristic_run_rins", 0, 1);

  AddHighsOption("mip:heur:rens heurrens mip_heuristic_run_rens",
    "0/1*: Use the RENS heuristic.",
    "mip_heuristic_run_rens", 0, 1);

  AddHighsOption("mip:heur:rootredcost heurrootredcost mip_heuristic_run_root_reduced_cost",
    "0/1*: Use the root reduced cost heuristic.",
    "mip_heuristic_run_root_reduced_cost", 0, 1);

  AddHighsOption("mip:heur:ziround heurziround mip_heuristic_run_zi_round",
    "0*/1: Use the ZI Round heuristic.",
    "mip_heuristic_run_zi_round", 0, 1);

  AddHighsOption("mip:heur:shifting heurshifting mip_heuristic_run_shifting",
    "0*/1: Use the Shifting heuristic.",
    "mip_heuristic_run_shifting", 0, 1);

  AddHighsOption("mip:lifting lifting mip_lifting_for_probing",
    "Level of lifting for probing that is used (default -1).",
    "mip_lifting_for_probing", -1, INT_MAX);

  AddHighsOption("lim:stallnodes stallnodelim stallnodelimit mip_max_stall_nodes",
    "Maximum MIP number of nodes where estimate is above cutoff bound (default: no limit).",
    "mip_max_stall_nodes", 0, INT_MAX);

  AddHighsOption("lim:leavenodes leaveslim mip_max_leaves",
    "Maximum MIP number of leaf nodes (default: no limit).",
    "mip_max_leaves", 0, INT_MAX);

  AddHighsOption("lim:nodes nodelim nodelimit mip_max_nodes",
    "Maximum MIP nodes to explore (default: no limit).",
    "mip_max_nodes", 0, INT_MAX);

  AddHighsOption("lim:sol sollim lim:improvingsols improvingsolslimit mip_max_improving_sols",
    "Stop the MIP solver after this many improving solutions "
    "have been found (default: no limit).",
    "mip_max_improving_sols", 1, INT_MAX);

  AddHighsOption("mip:lpagelimit lpagelimit mip_lp_age_limit",
    "Maximal age of dynamic LP rows before they are removed from the LP relaxation "
    "(default 10, max 32767).",
    "mip_lp_age_limit", 0, 32767);

  /*
  * mip_lp_solver
MIP LP solver: "choose", "simplex", "ipm", "ipx" or "hipo"
Type: string
Default: "choose"
mip_ipm_solver
MIP IPM solver: "choose", "ipx" or "hipo"
Type: string
Default: "choose"
ipm_optimality_tolera
  */

  AddHighsOption("mip:lpsolver miplpsolver mip_lp_solver",
      "Choose the lp solver to be used during mip solves:\n"
      "\n.. value-table::\n", "mip_lp_solver", mip_lp_solver_values, c);

  AddHighsOption("mip:ipmsolver mipipmsolver mip_ipm_solver",
      "Choose the ipm solver to be used during mip solves:\n"
      "\n.. value-table::\n", "mip_ipm_solver", mip_ipm_solver_values, c);



  AddHighsOption("mip:poolsoftlimit poolsoftlimit mip_pool_soft_limit",
    "Soft limit on the number of rows in the cutpool for dynamic age adjustment "
    "(default 10000).",
    "mip_pool_soft_limit", 1, INT_MAX);

  AddHighsOption("mip:pscostreliability pscostreliability mip_pscost_minreliable",
    "Minimal number of observations before pseudo costs are considered reliable "
    "(default 8).",
    "mip_pscost_minreliable", 0, INT_MAX);

  AddHighsOption("mip:mincliquetable mincliquetable mip_min_cliquetable_entries_for_parallelism",
    "Minimal number of entries in the cliquetable before neighborhood queries "
    "of the conflict graph use parallel processing (default 100000).",
    "mip_min_cliquetable_entries_for_parallelism", 0, INT_MAX);

  AddHighsOption("tech:miploglev miploglev mip_report_level",
    "0/1*/2: MIP solver report level",
    "mip_report_level", 0, 2);

  AddHighsOption("tech:miplogfreq miplogfreq mip_min_logging_interval",
    "Minimum interval between MIP log lines, in seconds (default 5).",
    "mip_min_logging_interval", 0.0, DBL_MAX);

  AddHighsOption("tech:logtoconsole logtoconsole log_to_console",
    "0/1*: Whether HiGHS log lines go to the console "
    "(with tech:outlev=1; set 0 to log only to tech:logfile).",
    "log_to_console", 0, 1);

  AddHighsOption("tech:timelesslog timelesslog timeless_log",
    "0*/1: Suppress time-based data in the HiGHS log.",
    "timeless_log", 0, 1);

  AddHighsOption("tech:seed seed random_seed",
    "Random number seed (default 0), affecting perturbations that "
    "may influence the solution path.",
    "random_seed", 0, INT_MAX);

  AddHighsOption("mip:intfeastol intfeastol inttol mip_feasibility_tolerance",
    "Feasibility tolerance for integer variables (default 1e-06).",
    "mip_feasibility_tolerance", 1e-10, Infinity());

  AddHighsOption("mip:heureff heureff mip_heuristic_effort",
    "Fraction of time to spend in MIP heuristics (default 0.05).",
    "mip_heuristic_effort", 0.0, 1.0);

  AddHighsOption("mip:gap mipgap mip:relgaptol relgaptol mip_rel_gap",
    "Tolerance on relative gap, | ub - lb|/|ub | , to determine whether optimality has been reached for a MIP instance "
    "(default 1e-04).",
    "mip_rel_gap", 0.0, Infinity());

  AddHighsOption("mip:gapabs mipgapabs mip:absgaptol absgaptol mip_abs_gap",
    "Tolerance on absolute gap of MIP, |ub-lb|, to determine whether optimality has been reached for a MIP instance "
    "(default 1e-06).",
    "mip_abs_gap", 0.0, Infinity());

  AddHighsOption("bar:centring centring pre:centring run_centring",
    "Perform centring steps or not:\n"
    "\n.. value-table::\n", "run_centring", values_01_noyes_0default_, 0);

  AddHighsOption("bar:maxcentringsteps maxcentringsteps pre:maxcentringsteps max_centring_steps",
    "Maximum number of steps to use when computing the analytic centre "
    "(default 5).",
    "max_centring_steps", 0, INT_MAX);

  AddHighsOption("bar:centringratiotol centringratiotol pre:centringratiotolerance centring_ratio_tolerance centringratiotolerance",
    "Centring stops when the ratio max(x_j*s_j) / min(x_j*s_j) is below "
    "this tolerance (default 100).",
    "centring_ratio_tolerance", 0.0, DBL_MAX);

  AddHighsOption("obj:blend blend_multi_objectives",
    "Whether to blend multiple objectives or apply lexicographical ordering",
    "blend_multi_objectives", values_01_noyes_1default_, 1
  );

  AddHighsOption("iis:strategy iis_strategy iisstrategy",
      "Strategy for IIS calculation (default 0), sum of:\n"
      "\n.. value-table::\n"
      "\nHiGHS presets: 6 = true IIS with row priority, "
      "14 = true IIS with column priority.",
	  "iis_strategy", iis_strategy_values, 0);

  AddHighsOption("lim:iistime iistimelim lim:iis_time_limit iis_time_limit iistimelimit",
	  "Time limit for computing IIS (default: no limit).",
	  "iis_time_limit", 0.0, DBL_MAX);
}

double HighsBackend::MIPGap() {
  if (BarrierIterations() == 0)
    return 0;
  return getDblAttr("mip_gap");
}
double HighsBackend::BestDualBound() {
  if (BarrierIterations() == 0)
    return 0;
  return getDblAttr("mip_dual_bound");
}

double HighsBackend::MIPGapAbs() {
  if (BarrierIterations() == 0)
    return 0;
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}
} // namespace mp


// AMPLs

AMPLS_MP_Solver* Open_highs(CCallbacks cb = {}) {
  AMPLS_MP_Solver* slv =
    AMPLS__internal__Open(
      std::unique_ptr<mp::BasicBackend>{new mp::HighsBackend()},
      cb);
  return slv;
}

void AMPLSClose_highs(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_highs(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::HighsBackend*>(AMPLSGetBackend(slv))->lp();
}


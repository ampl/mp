#include <vector>
#include <atomic>
#include <chrono>
#include <climits>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>

#include "google/protobuf/text_format.h"

#include "ortools/base/version.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/util.h"

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "cpsatbackend.h"

extern "C" {
  #include "cpsat-ampls-c-api.h"    // CP-SAT AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace sat = operations_research::sat;

namespace {

/// Called from the signal handler: must be async-signal-safe,
/// so it only sets an atomic flag. See InterruptWatcher.
bool InterruptCpsat(void* p) {
  static_cast<mp::CpsatData*>(p)->interrupted.store(true);
  return true;
}

/// While alive, stops the CP-SAT search when the interrupt flag is set,
/// or when the optional \a extra_stop returns true (polled likewise).
/// A thread polls the flag because the signal handler cannot call
/// sat::StopSearch() itself (it is not async-signal-safe).
/// The destructor joins the thread, also when the solve throws.
class InterruptWatcher {
public:
  InterruptWatcher(const std::atomic<bool>& flag, sat::Model* model,
                   std::function<bool()> extra_stop = {})
    : thread_([this, &flag, model, extra_stop = std::move(extra_stop)] {
        while (!done_.load()) {
          if (flag.load() || (extra_stop && extra_stop())) {
            sat::StopSearch(model);
            return;
          }
          std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
      }) { }
  ~InterruptWatcher() {
    done_.store(true);
    thread_.join();
  }
private:
  std::atomic<bool> done_{false};
  std::thread thread_;              // after done_: started last
};

/// Feeds CP-SAT's incumbents and bounds to MP's plateau logic, and
/// stops the search when it reports a plateau.
/// CP-SAT calls back from its worker threads, and the InterruptWatcher
/// polls CheckTimeout() from yet another one, while MP's PlateauState is
/// not thread-safe: all calls into the hooks are serialized here.
class PlateauMonitor {
public:
  /// The MP plateau functions (see MIPBackend::Report*ForPlateau()),
  /// each returning true if the plateau is reached
  struct Hooks {
    std::function<bool(double)> incumbent;
    std::function<bool(double, double)> gap;   // absolute, relative
    std::function<bool()> timeout;
  };

  /// Register the callbacks on \a model, which must outlive the solve
  PlateauMonitor(Hooks hooks, sat::Model* model)
    : hooks_(std::move(hooks)), model_(model) {
    model->Add(sat::NewFeasibleSolutionObserver(
        [this](const sat::CpSolverResponse& r) {
      StopIf(OnIncumbent(r.objective_value()));
    }));
    model->Add(sat::NewBestBoundCallback([this](double bound) {
      StopIf(OnBound(bound));
    }));
  }

  /// For polling between callbacks: has the plateau time run out?
  /// Does not stop the search, the caller does.
  bool CheckTimeout() { return Locked(hooks_.timeout); }

  /// Whether a plateau was reached
  bool stopped() const { return stopped_.load(); }

private:
  bool OnIncumbent(double obj) {
    return Locked([&] {
      haveIncumbent_ = true;
      incumbent_ = obj;
      return hooks_.incumbent(obj) || ReportGap();
    });
  }

  bool OnBound(double bound) {
    return Locked([&] {
      haveBound_ = true;
      bound_ = bound;
      return ReportGap();
    });
  }

  /// Call with the lock held. The gap is only known once both
  /// the incumbent and the bound have been seen.
  bool ReportGap() {
    if (!haveIncumbent_ || !haveBound_)
      return false;
    const double absgap = std::fabs(incumbent_ - bound_);
    return hooks_.gap(
          absgap, absgap / std::max(std::fabs(incumbent_), 1e-10));
  }

  template <class F>
  bool Locked(F&& f) {
    bool stop;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop = f();
    }
    if (stop)
      stopped_.store(true);
    return stop;
  }

  void StopIf(bool stop) {
    if (stop)
      sat::StopSearch(model_);
  }

  const Hooks hooks_;
  sat::Model* const model_;
  std::mutex mutex_;                  // guards the hooks and the fields below
  bool haveIncumbent_ = false, haveBound_ = false;
  double incumbent_ = 0.0, bound_ = 0.0;
  std::atomic<bool> stopped_{false};
};

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateCpsatBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CpsatBackend()};
}


namespace mp {

/// Create Cpsat Model Manager
/// @param gc: the Cpsat common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return CpsatModelMgr
std::unique_ptr<BasicModelManager>
CreateCpsatModelMgr(CpsatCommon&, Env&, pre::BasicValuePresolver*&);


CpsatBackend::CpsatBackend() {
  OpenSolver();

  pre::BasicValuePresolver* pPre;
  auto data = CreateCpsatModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy the shared data pointer to ModelAPI
  copy_common_info_to_other();
}

CpsatBackend::~CpsatBackend() {
  CloseSolver();
}

void CpsatBackend::OpenSolver() {
  set_data(new CpsatData());
  // Keep MP's SIGINT handler (see SetInterrupter()): by default CP-SAT
  // replaces it with its own during the solve. That one calls back through
  // a thread_local, which is empty where Windows runs signal handlers,
  // and the process aborts.
  data()->params.set_catch_sigint_signal(false);
}

void CpsatBackend::CloseSolver() {
  delete data();
  set_data(nullptr);
}

const char* CpsatBackend::GetBackendName()
  { return "CpsatBackend"; }

std::string CpsatBackend::GetSolverVersion() {
  return operations_research::OrToolsVersionString();
}


bool CpsatBackend::HasSolution() const {
  const auto st = data()->response.status();
  return st == sat::CpSolverStatus::OPTIMAL ||
         st == sat::CpSolverStatus::FEASIBLE;
}

ArrayRef<double> CpsatBackend::PrimalSolution() {
  std::vector<double> x;
  if (!HasSolution())
    return x;
  const auto* d = data();
  x.resize(d->vars.size());
  for (size_t i = 0; i < x.size(); ++i)
    x[i] = (double)sat::SolutionIntegerValue(d->response, d->vars[i]);
  return x;
}

double CpsatBackend::ObjectiveValue() const {
  return HasSolution() ? data()->response.objective_value()
                       : std::numeric_limits<double>::quiet_NaN();
}

double CpsatBackend::BestDualBound() {
  return data()->response.best_objective_bound();
}

double CpsatBackend::MIPGapAbs() {
  if (!HasSolution())
    return std::numeric_limits<double>::infinity();
  return std::fabs(ObjectiveValue() - BestDualBound());
}

double CpsatBackend::MIPGap() {
  if (!HasSolution())
    return std::numeric_limits<double>::infinity();
  const double obj = ObjectiveValue();
  return std::fabs(obj - BestDualBound()) / std::max(std::fabs(obj), 1e-10);
}


void CpsatBackend::DoWriteProblem(const std::string& file) {
  std::ofstream f(file);
  if (!f)
    throw std::runtime_error("Cannot open file '" + file + "' for writing");
  f << data()->builder.Proto().DebugString();
}

std::map<std::string, std::variant<int, double, std::string>>
CpsatBackend::SolutionStats() {
  const auto& r = data()->response;
  std::map<std::string, std::variant<int, double, std::string>> stats;
  stats["node_count"] = (double)r.num_branches();
  stats["conflicts"] = (double)r.num_conflicts();
  stats["simplex_iterations"] = 0.0;      // not applicable
  stats["barrier_iterations"] = 0.0;
  return stats;
}

void CpsatBackend::SetInterrupter(mp::Interrupter *inter) {
  // Reset before installing the handler, not in Solve(): an interrupt
  // arriving between the two must not be lost.
  data()->interrupted.store(false);
  inter->SetHandler(InterruptCpsat, data());
}

void CpsatBackend::Solve() {
  auto* d = data();
  pool_.clear();
  plateauStopped_ = false;
  // Once interrupted, no further solve (there can be several, e.g.,
  // for multiple objectives): report as not solved.
  if (d->interrupted.load()) {
    d->response = sat::CpSolverResponse();
    return;
  }
  sat::Model model;
  model.Add(sat::NewSatParameters(d->params));

  if (need_multiple_solutions()) {
    model.Add(sat::NewFeasibleSolutionObserver(
        [this, d](const sat::CpSolverResponse& r) {
      std::vector<double> x(d->vars.size());
      for (size_t i = 0; i < x.size(); ++i)
        x[i] = (double)sat::SolutionIntegerValue(r, d->vars[i]);
      pool_.emplace_back(std::move(x), r.objective_value());
    }));
  }

  // If we need the plateau observer, register it here
  // Note that the callbacks are registered in CP-SAT by the constructor:
  // PlateauMonitor(Hooks hooks, sat::Model* model)

  std::optional<PlateauMonitor> plateau;
  std::function<bool()> plateauTimeout;   // polled by the watcher, if set
  if (plateau_active()) {
    plateau.emplace(PlateauMonitor::Hooks{
          [this](double obj) { return ReportIncumbentForPlateau(obj); },
          [this](double abs, double rel) { return ReportGapForPlateau(abs, rel); },
          [this] { return CheckTimeoutForPlateau(); } },
        &model);
    plateauTimeout = [&plateau] { return plateau->CheckTimeout(); };
  }

  // Other watcher that stops the search on interrupt, or when the plateau time runs
  // out between callbacks
  InterruptWatcher watcher(d->interrupted, &model, plateauTimeout);

  d->response = sat::SolveCpModel(d->builder.Build(), &model);

  // If the plateau monitor was used, record whether it stopped the search
  plateauStopped_ = plateau && plateau->stopped();
}

void CpsatBackend::ReportResults() {
  ReportCpsatResults();
  BaseBackend::ReportResults();
}

void CpsatBackend::ReportCpsatResults() {
  SetStatus( GetSolveResult() );
  AddCpsatMessages();
  if (need_multiple_solutions())
    ReportCpsatPool();
}

void CpsatBackend::ReportCpsatPool() {
  for (const auto& s: pool_) {
    auto mv = GetValuePresolver().PostsolveSolution(
          { { s.first }, {}, std::vector<double>{ s.second } } );
    ReportIntermediateSolution(
          { mv.GetVarValues()(), mv.GetConValues()(), mv.GetObjValues()() });
  }
}

void CpsatBackend::AddCpsatMessages() {
  const auto& r = data()->response;
  // With several workers, CP-SAT reports no search counters in the
  // response (they are only in the search log), so don't print zeros
  if (r.num_booleans() > 0 || r.num_conflicts() > 0)
    AddToSolverMessage(fmt::format(
          "{} conflicts, {} branches, {:.3f} s\n",
          r.num_conflicts(), r.num_branches(), r.wall_time()));
  else
    AddToSolverMessage(fmt::format("{:.3f} s\n", r.wall_time()));
}

std::pair<int, std::string> CpsatBackend::GetSolveResult() {
  namespace sol = mp::sol;
  const auto* d = data();
  switch (d->response.status()) {
  case sat::CpSolverStatus::OPTIMAL:
    return { sol::SOLVED, "optimal solution" };
  case sat::CpSolverStatus::FEASIBLE:
    return d->interrupted.load()
        ? std::make_pair((int)sol::LIMIT_FEAS_INTERRUPT,
                         std::string("interrupted with a feasible solution"))
        : plateauStopped_.load()
        ? std::make_pair((int)sol::LIMIT_FEAS,
                         std::string("feasible solution, stopped on a plateau "
                                     "(mip:plateau* options)"))
        : std::make_pair((int)sol::LIMIT_FEAS,
                         std::string("feasible solution, not proven optimal "
                                     "(limit reached)"));
  case sat::CpSolverStatus::INFEASIBLE:
    return { sol::INFEASIBLE, "infeasible problem" };
  case sat::CpSolverStatus::MODEL_INVALID:
    return { sol::FAILURE, "invalid CP-SAT model: "
          + d->response.solution_info() };
  case sat::CpSolverStatus::UNKNOWN:
  default:
    return d->interrupted.load()
        ? std::make_pair((int)sol::LIMIT_NO_FEAS_INTERRUPT,
                         std::string("interrupted, no solution found"))
        : std::make_pair((int)sol::LIMIT_NO_FEAS,
                         std::string("no solution found (limit reached)"));
  }
}


void CpsatBackend::InitOptionParsing() {
    // Todo: experiment with best bounds. Seems to be better with
  ParseOptionString("cvt:pre:boundsbest=1");
}

void CpsatBackend::FinishOptionParsing() {
  auto* d = data();
  d->var_bound = storedOptions_.varbound_;
  d->params.set_log_search_progress(storedOptions_.outlev_ > 0);
  if (!storedOptions_.params_.empty()) {
    if (!google::protobuf::TextFormat::MergeFromString(
          storedOptions_.params_, &d->params))
      throw std::runtime_error(
          "Could not parse 'params' as CP-SAT SatParameters "
          "(protobuf text format): " + storedOptions_.params_);
  }
}


////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo values_linearization_[] = {
  {"0", "No linearization", 0},
  {"1", "Linearize some constraints (default)", 1},
  {"2", "Full linearization", 2} };

void CpsatBackend::InitCustomOptions() {

  set_option_header(
      "CP-SAT Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``cpsat_options``. For example::\n"
      "\n"
      "  ampl: option cpsat_options 'timelim=60 threads=4';\n"
      "\n"
      "CP-SAT works with integer variables and integer coefficients only. "
      "Alldiff, count, numberof, if-then-else, min, max, abs, "
      "and/or/indicator constraints are passed natively.\n");

  AddStoredOption("tech:outlev outlev",
      "0*/1: Whether to write the CP-SAT search log to stdout.",
      storedOptions_.outlev_, 0, 1);

  AddStoredOption("tech:params params parameters",
      "Additional CP-SAT parameters, in the protobuf text format of "
      "SatParameters, applied after all other options. Example: "
      "``params='linearization_level: 2 symmetry_level: 3'``.",
      storedOptions_.params_);

  AddStoredOption("tech:varbound varbound",
      "Magnitude used for infinite variable bounds, which CP-SAT cannot "
      "handle (default 1e12). Larger values may cause integer overflow "
      "in the CP-SAT model.",
      storedOptions_.varbound_, 1.0, 1e18);

  AddSolverOption("tech:threads threads",
      "How many search workers (threads) to use; "
      "default 0 ==> automatic choice.",
      "num_workers", 0, INT_MAX);

  AddSolverOption("tech:seed seed",
      "Random seed (default 1).",
      "random_seed", 0, INT_MAX);

  AddSolverOption("lim:time timelim timelimit",
      "Limit on solve time (in seconds; default: no limit).",
      "max_time_in_seconds", 0.0, DBL_MAX);

  AddSolverOption("lim:dettime dettimelim dettimelimit",
      "Limit on the deterministic time (an abstract, machine-independent "
      "measure, which makes runs reproducible); default: no limit.",
      "max_deterministic_time", 0.0, DBL_MAX);

  AddSolverOption("lim:conflicts conflictlim conflictlimit",
      "Limit on the number of conflicts; default: no limit.",
      "max_number_of_conflicts", 0, INT_MAX);

  AddSolverOption("lim:mem memlimit",
      "Limit on the memory used, in megabytes (default 10000). "
      "The solver stops when it is reached.",
      "max_memory_in_mb", 1, INT_MAX);

  AddSolverOption("mip:gap mipgap",
      "Relative gap between the best objective and the best bound "
      "at which to stop (default 0).",
      "relative_gap_limit", 0.0, DBL_MAX);

  AddSolverOption("mip:gapabs mipgapabs",
      "Absolute gap between the best objective and the best bound "
      "at which to stop (default 1e-4).",
      "absolute_gap_limit", 0.0, DBL_MAX);

  AddSolverOption("pre:solve presolve",
      "Whether to use CP-SAT presolve:\n"
      "\n.. value-table::\n",
      "cp_model_presolve", values_01_noyes_1default_, 1);

  AddSolverOption("alg:linearization linearization_level",
      "Level of linearization of the CP model (default 1):\n"
      "\n.. value-table::\n",
      "linearization_level", values_linearization_, 1);
}

}  // namespace mp


// AMPLs

AMPLS_MP_Solver* Open_cpsat(CCallbacks cb = {}) {
  return AMPLS__internal__Open(
      std::unique_ptr<mp::BasicBackend>{new mp::CpsatBackend()}, cb);
}

void AMPLSClose_cpsat(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_cpsat(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::CpsatBackend*>(AMPLSGetBackend(slv))->data();
}

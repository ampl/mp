#ifndef BACKEND_PLATEAU_STATE_H
#define BACKEND_PLATEAU_STATE_H

#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#include "mp/format.h"

namespace mp {

/// Bookkeeping for the PLATEAU_STOP / PLATEAU_STOP_BOUND features.
/// Tracks time elapsed since the last "sufficient" improvement of the
/// incumbent objective (and, if enabled, the best bound).
struct PlateauState {

  // State
  bool warmingUp_ = true;
  bool stopping_ = false;
  bool active_ = false;
  bool haveIncumbent_ = false;
  bool haveGap_ = false;
  double bestObj_ = 0.0;
  double bestAbsGap_ = 0.0;
  double bestRelGap_ = 0.0;
  std::chrono::steady_clock::time_point lastImprovement_;

  /// Options
  struct Options {
    double warmup_time_ = 0.0;        // mip:plateauwarmup, seconds
    double warmup_absgap_ = 0.0;      // mip::plateauwarmupabsgap; 0=disabled
    double warmup_relgap_ = 0.0;      // mip::plateauwarmuprelgap; 0=disabled
    double plateau_time_ = 0.0;       // mip:plateautime, seconds; 0 = disabled
    double abstol_ = 0.0;             // mip:plateauabstol
    double reltol_ = 0.0;             // mip:plateaureltol
    double absmipgap_tol_ = 0.0;      // mip:plateauabsgaptol: 0=disabled, >0=also track abs gap
    double relmipgap_tol_ = 0.0;      // mip:plateaurelgaptol: 0=disabled, >0=also track rel gap
    int log_ = false;                 // mip:plateaulog: 0=silent, 1=log every callback
  };

  Options opts_dflt_, opts_;

  /// Save current options as default
  void SaveOptions() { opts_dflt_ = opts_; }

  // Per-pass overrides for native multi-objective solving, indexed
  // directly by pass number (as assigned by MOManager -- for Gurobi
  // this is also the native GRBgetmultiobjenv() objective index).
  // Empty (size 0) means no per-pass overrides were captured at all
  // (single-objective solve, or emulated multiobj -- which instead
  // gets a fresh Init() per pass from MIPBackend::SetupPlateau()).
  // std::vector<double> warmupByPass_, warmupAbsGapByPass_, warmupRelGapByPass_,
  //     plateauTimeByPass_, abstolByPass_, reltolByPass_,
  //     absGaptolByPass_, relGaptolByPass_;
  // std::vector<bool> logByPass_;

  void Init(double warmup_time, double warmup_absgap, double warmup_relgap,
            double plateauTimeVal,
            double abstol, double reltol, double absmipgap_tol, double relmipgap_tol,
            bool log) {
    Reset();
    opts_.warmup_time_ = warmup_time;
    opts_.warmup_absgap_ = warmup_absgap;
    opts_.warmup_relgap_ = warmup_relgap;
    opts_.plateau_time_ = plateauTimeVal;
    opts_.abstol_ = abstol;
    opts_.reltol_ = reltol;
    opts_.absmipgap_tol_ = absmipgap_tol;
    opts_.relmipgap_tol_ = relmipgap_tol;
    opts_.log_ = log;
  }

  /// Also resets options
  void Reset() {
    active_ = true;
    warmingUp_ = true;
    stopping_ = false;
    haveIncumbent_ = haveGap_ = false;
    lastImprovement_ = std::chrono::steady_clock::now();

    opts_ = opts_dflt_;
  }

  /// Compares val against the reference point `best` (the value as of
  /// the last reset). Only when the change clears the tolerance does
  /// it count as "progress" and `best`/`lastImprovement_` advance;
  /// otherwise they're left untouched so that a sequence of many
  /// small, sub-threshold changes keeps accumulating against the
  /// same baseline instead of resetting it step by step. Also logs
  /// (if mip:plateau:log) and evaluates the shared stop decision.
  /// @return true if the solve should be terminated now (plateau reached).
  bool Improved(const char* kind, double& best, double val,
                double absTol, double relTol) {
    double old_best = best;
    auto old_last_improvement = lastImprovement_;
    double absDelta = std::fabs(best - val);
    double relDelta = 0.0 != best
                          ? absDelta / std::fabs(best)
                          : (absDelta > 0.0
                                 ? std::numeric_limits<double>::infinity() : 0.0);
    bool improved = absDelta > absTol || relDelta > relTol;
    if (improved) {
      best = val;
      lastImprovement_ = std::chrono::steady_clock::now();
    }
    stopping_ = CheckTimeout();
    if (opts_.log_) {
      double sinceProgress =
          std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() -
                                 old_last_improvement).count();
      LogStatus(kind, val, old_best, improved, stopping_, sinceProgress);
    }
    return stopping_;
  }

  bool IsWarmupDone(double current_absgap, double current_relgap)  {
    if (!warmingUp_ || stopping_) return true;

    if ((opts_.warmup_absgap_ > 0.0 && current_absgap <= opts_.warmup_absgap_) ||
        (opts_.warmup_relgap_ > 0.0 && current_relgap <= opts_.warmup_relgap_))
    {
      warmingUp_ = false;
      lastImprovement_ = std::chrono::steady_clock::now();
      if (opts_.log_)
        fmt::print("    MP Plateau: Warmup gap reached.\n");
      return true;
    }
    // if condition on gap is not reached, check time
    auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration<double>(now - lastImprovement_).count() >=
        opts_.warmup_time_)
    {
      warmingUp_ = false;
      lastImprovement_ = std::chrono::steady_clock::now();
      if (opts_.log_)
        fmt::print("    MP Plateau: Warmup time ({:.1f}s) reached.\n",
                   opts_.warmup_time_);
      return true;

    }
    return false;
  }

  bool CheckTimeout() {
    // Gate on warmup here too (not just in ReportIncumbent/ReportGap):
    // this is also reachable directly via CheckTimeoutForPlateau(),
    // e.g. from a bare periodic callback (Gurobi's GRB_CB_POLLING)
    // that never calls ReportIncumbent/ReportGap in between. Without
    // this, mip:plateauwarmup* would not protect such a callback from
    // triggering a stop before warmup ends. Gap is unknown here, so
    // pass +inf to disable the gap-based early-exit and fall back to
    // the plain elapsed-time check.
    if (stopping_) return true;

    constexpr double kInf = std::numeric_limits<double>::infinity();
    if (!IsWarmupDone(kInf, kInf))
      return false;
    auto now = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration<double>(now - lastImprovement_).count();
    stopping_ = diff >= opts_.plateau_time_;

    if (stopping_) {
      fmt::print("    MP Plateau: stopping after {:.1f}s "
                 "without sufficient improvement (limit={:.1f}s)\n",
                 diff, opts_.plateau_time_);
    }
    return stopping_;
  }

  /// Print a one-line status snapshot; called (if mip:plateaulog=1)
  /// on every ReportIncumbent()/ReportBound() call, i.e. on every
  /// native callback invocation that reaches the plateau logic.
  void LogStatus(const char* kind, double val, double baseline,
                 bool progress, bool stop,
                 double sinceProgress) const {
    // auto now = std::chrono::steady_clock::now();
    fmt::print(
        "    MP Plateau [{}]: current={:.6g} previous={:.6g} progress={} "
        "elapsed={:.1f}s (limit={:.1f}s){}\n",
        kind, val, baseline, progress ? "yes" : "no",
        sinceProgress, opts_.plateau_time_,
        stop ? " -- STOPPING (plateau reached)" : "");
    std::fflush(stdout);
  }

  void SetCurrentObjective(int nobj) {
    Reset();
    if (opts_.log_)
    {
      fmt::print("   MP Plateau: new pass {}, resetting\n",
                 nobj + 1);
      std::fflush(stdout);
    }
  }

  bool ReportIncumbent(double obj) {
    if (!active_) return false;
    if (stopping_) return true;

    // Cannot assume I know the new MIP gap here
    auto kInf = std::numeric_limits<double>::infinity();
    if (!IsWarmupDone(kInf, kInf))
      return false;

    if (!haveIncumbent_) {
      haveIncumbent_ = true;
      bestObj_ = obj;
      auto now = std::chrono::steady_clock::now();
      double sinceProgress =
          std::chrono::duration<double>(now - lastImprovement_).count();
      lastImprovement_ = now  ;
      stopping_ = CheckTimeout();
      if (opts_.log_)
        LogStatus("incumbent", obj, bestObj_,
                  true, stopping_, sinceProgress);
      return stopping_;
    }
    return Improved("incumbent", bestObj_, obj,
                    opts_.abstol_, opts_.reltol_);
  }

  bool ReportGap(double absgap, double relgap) {
    if (!active_) return false;
    if (stopping_) return true;
    if (!IsWarmupDone(absgap, relgap))
      return false;

    if (!haveGap_) {
      haveGap_ = true;
      bestAbsGap_ = absgap;
      bestRelGap_ = relgap;
      auto now = std::chrono::steady_clock::now();
      double sinceProgress =
          std::chrono::duration<double>(now - lastImprovement_).count();
      lastImprovement_ = now;
      stopping_ = CheckTimeout();
      if (opts_.log_) {
        LogStatus("absmipgap", absgap, bestAbsGap_,
                  true, stopping_, sinceProgress);
        LogStatus("relmipgap", relgap, bestRelGap_,
                  true, stopping_, sinceProgress);
      }
      return stopping_;
    }

    // Compare only if abs or gap limit given
    // Any of the two works as a sufficient condition for improvement
    constexpr double kInf = std::numeric_limits<double>::infinity();
    if (opts_.absmipgap_tol_ > 0.0)
      if (Improved("absmipgap", bestAbsGap_, absgap,
                   opts_.absmipgap_tol_, kInf))
        return true;
    if (opts_.relmipgap_tol_ > 0.0)
      if (Improved("relmipgap", bestRelGap_, relgap,
                   opts_.relmipgap_tol_, kInf))
        return true;

    return CheckTimeout();
  }
};

}  // namespace mp

#endif // BACKEND_PLATEAU_STATE_H

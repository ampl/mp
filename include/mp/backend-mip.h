/*
 Abstract MIP solver backend wrapper.

 Copyright (C) 2022 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 */

#ifndef MIPBACKEND_H_
#define MIPBACKEND_H_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#include "mp/common.h"
#include "mp/backend-std.h"

namespace mp {

/// Basis status values of a solution (postsolved)
struct SolutionBasis {
  /// Check if has both vars and cons' statuses
  operator bool() const { return varstt.size() && constt.size(); }
  /// Var and con statuses
  std::vector<int> varstt, constt;
};


/// IIS (postsolved).
/// Elements correspond to IISStatus
struct IIS {
  /// Var and con IIS statuses
  std::vector<int> variis, coniis;
};


/// Sensitivity ranges (postsolved)
struct SensRanges {
  std::vector<double>
      varlblo, varlb, varlbhi,         // varlb/ub, conlb/ub/rhs not needed
      varublo, varub, varubhi,
      varobjlo, varobj, varobjhi,      // varobj is compulsory
      conrhslo, conrhs, conrhshi,      // for rhs-constraints
      conlblo, conlb, conlbhi,
      conublo, conub, conubhi;         // for range constraints
};


/// MIP backend wrapper
///
/// The MIP wrapper provides common functionality relative to MIP solvers;
/// it implements the common suffixes and the logic shared across all MIP
/// solvers
template <class Impl,
          class BaseBackend = StdBackend<Impl> >  ///< parameter for base class
class MIPBackend : public BaseBackend
{
public:
  /// IsMIP().
  /// Basic version: does not consider PL or SOS or Q(C)P.
  bool IsMIP() const override
  { return BackendWithModelManager::HasUnfixedIntVars(); }
  /// IsQP()
  virtual bool IsQP() const { return false; }
  /// IsQCP()
  virtual bool IsQCP() const { return false; }

  /// Always add MIP start if supported:
  /// Gurobi 9.1.2 solves non-convex Q(C)P as MIP.
  /// But model attributes don't work before solve.
  /// Distinguish LP / non-LP instead?
  virtual bool CanBeMIP() const { return true; }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  /////////////// Most are disabled by default ///////////////
  //// To enable, declare ALLOW_STD_FEATURE( name, true ) ////
  /////////// and implement the relevant methods /////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;
  /**
  * Set lazy/user cut attributes.
  * Negative suffix values are "user cuts".
  * Check lazy_/user_cuts() to see which kinds are allowed.
  * Presolve the values if needed.
  **/
  DEFINE_STD_FEATURE( LAZY_USER_CUTS )
  ALLOW_STD_FEATURE( LAZY_USER_CUTS, false )
  virtual void MarkLazyOrUserCuts(ArrayRef<int> ) { }
  /**
  * Get/Set AMPL var/con statii
  **/
  DEFINE_STD_FEATURE( BASIS )
  ALLOW_STD_FEATURE( BASIS, false )
  /// The basis statuses of vars and cons.
  /// MIPBackend handles them in postsolved form (for the NL model)
  /// Impl has to perform value pre- / postsolve if needed
  /// Getter (returns postsolved basis)
  virtual SolutionBasis GetBasis() { return {}; }
  /// Setter (takes unpresolved basis)
  virtual void SetBasis(SolutionBasis )
  { MP_UNSUPPORTED("MIPBackend::SetBasis"); }
  /**
  * General LP warm start, e.g.,
  * set primal+dual initial guesses for continuous case.
  * The specific Backend should
  * presolve the values if needed.
  *
  * @note Only called when both primal and dual initial guesses
  *   are provided. Thus, use MIPSTART for primal-only
  *   starts (MIP, NLP).
  **/
  DEFINE_STD_FEATURE( WARMSTART )
  ALLOW_STD_FEATURE( WARMSTART, false )
  virtual void AddPrimalDualStart(Solution )
  { MP_UNSUPPORTED("MIPBackend::AddPrimalDualStart"); }
  /**
  * Primal-only (MIP / NLP) warm start.
  * Provides solution hints (dense vector),
  * as well as sparsity pattern (dense 0-1 vector),
  * allowing partial MIP / NLP warm start.
  * Presolve the values if needed.
  **/
  DEFINE_STD_FEATURE( MIPSTART )
  ALLOW_STD_FEATURE( MIPSTART, false )
  virtual void AddMIPStart(
      ArrayRef<double> , ArrayRef<int> )
  { MP_UNSUPPORTED("MIPBackend::AddMIPStart"); }
  /**
  * Set branch and bound priority
  **/
  DEFINE_STD_FEATURE( VAR_PRIORITIES )
  ALLOW_STD_FEATURE( VAR_PRIORITIES, false )
  virtual void VarPriorities(ArrayRef<int>)
  { MP_UNSUPPORTED("MIPBackend::VarPriorities"); }
  /**
  * Obtain unbounded/inf rays
  **/
  DEFINE_STD_FEATURE( RAYS )
  ALLOW_STD_FEATURE( RAYS, false )
  virtual ArrayRef<double> Ray() { return {}; }
  virtual ArrayRef<double> DRay() { return {}; }
  /**
  * Compute the IIS and obtain relevant values (postsolved)
  **/
  DEFINE_STD_FEATURE( IIS )
  ALLOW_STD_FEATURE( IIS, false )
  virtual void ComputeIIS() {}
  virtual IIS GetIIS() { return {}; }
  /**
  * Get MIP Gap
  **/
  DEFINE_STD_FEATURE( RETURN_MIP_GAP )
  ALLOW_STD_FEATURE( RETURN_MIP_GAP, false )
  /// Should return AMPLInf() if not available
  virtual double MIPGap() { return MP_DISPATCH( AMPLInf() ); }
  /// Should return AMPLInf() if not available
  virtual double MIPGapAbs() { return MP_DISPATCH( AMPLInf() ); }
  /**
  * Get MIP dual bound
  **/
  DEFINE_STD_FEATURE( RETURN_BEST_DUAL_BOUND )
  ALLOW_STD_FEATURE( RETURN_BEST_DUAL_BOUND, false )
  virtual double BestDualBound()
  { MP_UNSUPPORTED("BestDualBound()"); return 0.0; }
  /**
  * Report sensitivity analysis suffixes
  **/
  DEFINE_STD_FEATURE( SENSITIVITY_ANALYSIS )
  ALLOW_STD_FEATURE( SENSITIVITY_ANALYSIS, false )
  virtual SensRanges GetSensRanges() { return {}; }
  /**
  * FixModel - duals, basis, and sensitivity for MIP
  * No API to overload,
  * Impl should check need_fixed_MIP()
  **/
  DEFINE_STD_FEATURE( FIX_MODEL )
  ALLOW_STD_FEATURE( FIX_MODEL, false )
  /**
  * Stop the MIP search once the incumbent objective (and, if
  * PLATEAU_STOP_BOUND is also enabled, the absolute/relative MIP gap)
  * has not improved enough for mip:plateautime seconds.
  * Implement by registering, in one override of SetupPlateauCallbacks(),
  * the solver's native new-incumbent callback and any other native
  * callback(s) needed to obtain the current gap, and from them call:
  * - ReportIncumbentForPlateau(objVal) with each new incumbent
  *   objective value;
  * - ReportGapForPlateau(absgap, relgap) with the current absolute and
  *   relative MIP gap, whenever available (only meaningful together
  *   with PLATEAU_STOP_BOUND: derive both from whatever
  *   incumbent/bound values the native API provides, e.g.
  *   absgap=|obj-bound|, relgap=absgap/|obj|, unless the solver already
  *   computes them, like HiGHS's mip_gap);
  * - CheckTimeoutForPlateau() from any other periodic/polling callback
  *   that doesn't have a fresh incumbent or gap value to report (e.g.
  *   Gurobi's GRB_CB_POLLING), just to keep checking the shared clock.
  * If any of these return true, terminate the native solve right there.
  **/
  DEFINE_STD_FEATURE( PLATEAU_STOP )
  ALLOW_STD_FEATURE( PLATEAU_STOP, false )
  DEFINE_STD_FEATURE ( PLATEAU_STOP_BOUND )
  ALLOW_STD_FEATURE( PLATEAU_STOP_BOUND, false )
  virtual void SetupPlateauCallbacks()
  { MP_UNSUPPORTED("MIPBackend::SetupPlateauCallbacks"); }


  ////////////////////////////////////////////////////////////////////////////
  /////////////////////// MIP specific derived calculations //////////////////
  ////////////////////////////////////////////////////////////////////////////

  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         INPUT         //////////////////////////
  using BaseBackend::ReadSuffix;
  using BaseBackend::ReadIntSuffix;
  using BaseBackend::ReadDblSuffix;
  using BaseBackend::ReportSuffix;
  using BaseBackend::ReportIntSuffix;
  using BaseBackend::ReportDblSuffix;

  void InputExtras() override {
    BaseBackend::InputExtras();
    InputMIPExtras();
  }

  virtual void InputMIPExtras() {
    if (lazy_user_cuts())
      InputLazyUserCuts();
    InputStartValues();
    if (priorities()) {
      if (auto pri_array = ReadSuffix(suf_varpriority))
        VarPriorities( pri_array );
    }
  }

  virtual void InputLazyUserCuts() {
    auto sufLazyVal = ReadIntSuffix( {"lazy", suf::CON} );
    if (sufLazyVal)
      MarkLazyOrUserCuts(sufLazyVal);
  }

  /// Report a suffix on the problem.
  /// Mainly used for testing.
  virtual void ReportProblemSuffix(const char* suf_name, int val) {
    ReportIntSuffix({suf_name, suf::PROBLEM}, {{val}});
  }

  virtual void InputStartValues() {
    InputPrimalDualStartOrBasis(); /// Always
    if ( CanBeMIP() ) {
      InputMIPStart();
    }
  }

  virtual void InputPrimalDualStartOrBasis() {
    bool useBasis = need_basis_in();
    SolutionBasis basis;
    if (useBasis) {
      basis.varstt = ReadSuffix(suf_varstatus);
      basis.constt = ReadSuffix(suf_constatus);
      useBasis = bool(basis);
    }
    Solution sol0;           // initial guesses
    sol0.primal = this->InitialValues();
    sol0.spars_primal = this->InitialValuesSparsity();
    sol0.dual = this->InitialDualValues();
    bool haveInis = sol0.primal.size() && sol0.dual.size();
    if (haveInis && (
          2<=warmstart() ||
          (1==warmstart() && !useBasis))) {
      AddPrimalDualStart(sol0);
      if ( 2==warmstart() )  // Why should we submit only the warmstart?
        useBasis = false;
      if (debug_mode()) {                 // Report received initials
        ReportSuffix(suf_testvarini, sol0.primal); // Should we check that
        ReportSuffix(suf_testconini, sol0.dual);   // Impl uses them?
      }
    }
    if (useBasis) {
      SetBasis(basis);
      if (debug_mode()) {                    // Report received statuses
        ReportSuffix(suf_testvarstatus, basis.varstt); // Should we check that
        ReportSuffix(suf_testconstatus, basis.constt); // Impl uses them?
      }
    }
  }

  virtual void InputMIPStart() {
    if (warmstart() && this->InitialValues().size() > 0) {
      if (IMPL_HAS_STD_FEATURE( MIPSTART )) {
        AddMIPStart(
              this->InitialValues(),
              this->InitialValuesSparsity() );
        if (debug_mode()) {                    // Report received initials
          ReportSuffix(suf_testMIPini,         // Should we check that
                       this->InitialValues());
        }                                      // Impl uses them?
      }
    }
  }


  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         OUtPUT        //////////////////////////
  void ReportStandardSuffixes() override {
    BaseBackend::ReportStandardSuffixes();
    ReportStandardMIPSuffixes();
  }

  virtual void ReportStandardMIPSuffixes() {
    if (need_basis_out())
      ReportBasis();
    ReportRays();
    CalculateAndReportIIS();
    if (IsMIP())
      CalculateAndReportMIPGap();
    ReportBestDualBound();
    if (sensitivity())
      ReportSensitivity();
  }

  virtual void ReportBasis() {
    /// Rely on solver reporting both vectors only if valid basis exists
    if (auto basis = GetBasis()) {
      ReportSuffix(suf_varstatus, basis.varstt);
      ReportSuffix(suf_constatus, basis.constt);
    }
  }

  virtual void ReportRays() {
    if ( need_ray_primal() &&
         ( this->IsProblemUnbounded() ||
           this->IsProblemIndiffInfOrUnb() )) {
      ReportSuffix(suf_unbdd, Ray() );
    }
    if ( need_ray_dual() &&
         ( this->IsProblemInfeasible() ||
           this->IsProblemIndiffInfOrUnb() )) {
      ReportSuffix(suf_dunbdd, DRay() );
    }
  }

  virtual void CalculateAndReportIIS() {
    if (( this->IsProblemInfeasible() ||
               this->IsProblemIndiffInfOrUnb() ) &&
        GetMIPOptions().exportIIS_) {
      try {
        ComputeIIS();
      } catch (const std::exception& exc) {
        this->AddWarning("IIS_COMPUTE",   // Can add warning before SOL output
                   std::string("Error computing IIS: ")
                       + exc.what());
      }
      this->SetStatus( this->GetSolveResult() );

      if (this->IsProblemInfeasible()) {      // can be unbounded
        auto iis = GetIIS();
        ReportSuffix(sufIISCon, iis.coniis);
        ReportSuffix(sufIISVar, iis.variis);
      }
    }
  }

  virtual void CalculateAndReportMIPGap() {
    std::vector<double> dbl(1);
    if (1 & GetMIPOptions().returnMipGap_) {
      dbl[0] = MP_DISPATCH( MIPGap() );
      ReportSuffix(sufRelMipGapObj, dbl);
      ReportSuffix(sufRelMipGapProb, dbl);
    }
    if (2 & GetMIPOptions().returnMipGap_) {
      dbl[0] = MP_DISPATCH( MIPGapAbs() );
      ReportSuffix(sufAbsMipGapObj, dbl);
      ReportSuffix(sufAbsMipGapProb, dbl);
    }
    if (!(GetMIPOptions().returnMipGap_ & 4)) {
      double absMIPGap = MP_DISPATCH(MIPGapAbs());
      if(absMIPGap > 0. && absMIPGap < MP_DISPATCH(Infinity()))
        BaseBackend::AddToSolverMessage(
              fmt::format("absmipgap={}, relmipgap={}",
                          absMIPGap, MP_DISPATCH(MIPGap())));
    }
  }

  virtual void ReportBestDualBound() {
    if (GetMIPOptions().returnBestDualBound_) {
      std::vector<double> dbl(1, MP_DISPATCH( BestDualBound() ));
      ReportSuffix(sufBestBoundObj, dbl);
      ReportSuffix(sufBestBoundProb, dbl);
    }
  }

  virtual void ReportSensitivity() {
    SensRanges sensr = GetSensRanges();
    ReportSuffix( {"senslbhi", suf::Kind::VAR}, sensr.varlbhi );
    // ReportSuffix( {"senslb", suf::Kind::VAR}, sensr.varlb );
    ReportSuffix( {"senslblo", suf::Kind::VAR}, sensr.varlblo );
    ReportSuffix( {"sensubhi", suf::Kind::VAR}, sensr.varubhi );
    // ReportSuffix( {"sensub", suf::Kind::VAR}, sensr.varub );
    ReportSuffix( {"sensublo", suf::Kind::VAR}, sensr.varublo );
    ReportSuffix( {"sensobjhi", suf::Kind::VAR}, sensr.varobjhi );
    ReportSuffix( {"up", suf::Kind::VAR}, sensr.varobjhi );   // CPLEXASL
    ReportSuffix( {"sensobj", suf::Kind::VAR}, sensr.varobj );
    ReportSuffix( {"current", suf::Kind::VAR}, sensr.varobj );
    ReportSuffix( {"sensobjlo", suf::Kind::VAR}, sensr.varobjlo );
    ReportSuffix( {"down", suf::Kind::VAR}, sensr.varobjlo );
    ReportSuffix( {"sensrhshi", suf::Kind::CON}, sensr.conrhshi );
    ReportSuffix( {"up", suf::Kind::CON}, sensr.conrhshi );
    // ReportSuffix( {"sensrhs", suf::Kind::CON}, sensr.conrhs );
    ReportSuffix( {"sensrhslo", suf::Kind::CON}, sensr.conrhslo );
    ReportSuffix( {"down", suf::Kind::CON}, sensr.conrhslo );
    ReportSuffix( {"senslbhi", suf::Kind::CON}, sensr.conlbhi );
    // ReportSuffix( {"senslb", suf::Kind::CON}, sensr.conlb );
    ReportSuffix( {"senslblo", suf::Kind::CON}, sensr.conlblo );
    ReportSuffix( {"sensubhi", suf::Kind::CON}, sensr.conubhi );
    // ReportSuffix( {"sensub", suf::Kind::CON}, sensr.conub );
    ReportSuffix( {"sensublo", suf::Kind::CON}, sensr.conublo );
  }

  ////////////////////////////////////////////////////////////
  /////////////////// MIP Backend options ////////////////////
  ////////////////////////////////////////////////////////////
private:
  struct Options {
    int lazy_user_cuts_ = 3;
    int basis_=3;
    int warmstart_=3;
    int importPriorities_=1;
    int rays_=3;
    int exportIIS_=0;
    int returnMipGap_=0;
    int returnBestDualBound_=0;
    int solnSens_=0;
    int fixModel_=0;
  };
  Options mipStoredOptions_;


  /// Bookkeeping for the PLATEAU_STOP / PLATEAU_STOP_BOUND features.
  /// Tracks time elapsed since the last "sufficient" improvement of the
  /// incumbent objective (and, if enabled, the best bound).
  struct PlateauState {

      // State
      bool timeoutPassed_ = false;
      bool active_ = false;
      bool haveIncumbent_ = false;
      bool haveGap_ = false;
      double bestObj_ = 0.0;
      double bestAbsGap_ = 0.0;
      double bestRelGap_ = 0.0;
      std::chrono::steady_clock::time_point lastImprovement_;

      // Options
      double warmup_time_ = 0.0;        // mip:plateauwarmup, seconds
      double warmup_absgap_ = 0.0;      // mip::plateauwarmupabsgap; 0=disabled
      double warmup_relgap_ = 0.0;      // mip::plateauwarmuprelgap; 0=disabled
      double plateau_time_ = 0.0;       // mip:plateautime, seconds; 0 = disabled
      double abstol_ = 0.0;             // mip:plateauabstol
      double reltol_ = 0.0;             // mip:plateaureltol
      double absmipgap_tol_ = 0.0;      // mip:plateauabsgaptol: 0=disabled, >0=also track abs gap
      double relmipgap_tol_ = 0.0;      // mip:plateaurelgaptol: 0=disabled, >0=also track rel gap
      int log_ = false;                 // mip:plateaulog: 0=silent, 1=log every callback

      
      
      // Per-pass overrides for native multi-objective solving, indexed
      // directly by pass number (as assigned by MOManager -- for Gurobi
      // this is also the native GRBgetmultiobjenv() objective index).
      // Empty (size 0) means no per-pass overrides were captured at all
      // (single-objective solve, or emulated multiobj -- which instead
      // gets a fresh Init() per pass from MIPBackend::SetupPlateau()).
      std::vector<double> warmupByPass_, warmupAbsGapByPass_, warmupRelGapByPass_,
          plateauTimeByPass_, abstolByPass_, reltolByPass_,
          absGaptolByPass_, relGaptolByPass_;
      std::vector<bool> logByPass_;

      void Init(double warmup_time, double warmup_absgap, double warmup_relgap,
          double plateauTimeVal,
          double abstol, double reltol, double absmipgap_tol, double relmipgap_tol,
          bool log) {
          warmup_time_ = warmup_time;
          warmup_absgap_ = warmup_absgap;
          warmup_relgap_ = warmup_relgap;
          plateau_time_ = plateauTimeVal;
          abstol_ = abstol;
          reltol_ = reltol;
          absmipgap_tol_ = absmipgap_tol;
          relmipgap_tol_ = relmipgap_tol;
          log_ = log;
          Reset();
      }
      void Reset() {
          active_ = true;
          timeoutPassed_ = false;
          haveIncumbent_ = haveGap_ = false;
          lastImprovement_ = std::chrono::steady_clock::now();
          
      }

      /// Compares val against the reference point `best` (the value as of
      /// the last reset). Only when the change clears the tolerance does
      /// it count as "progress" and `best`/`lastImprovement_` advance;
      /// otherwise they're left untouched so that a sequence of many
      /// small, sub-threshold changes keeps accumulating against the
      /// same baseline instead of resetting it step by step. Also logs
      /// (if mip:plateaulog) and evaluates the shared stop decision.
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
          bool stop = CheckTimeout();
          if (log_) {
              double sinceProgress =
                  std::chrono::duration<double>(std::chrono::steady_clock::now() - old_last_improvement).count();
              LogStatus(kind, val, old_best, improved, stop, sinceProgress);
          }
          return stop;
      }

      bool IsWarmupDone(double current_absgap, double current_relgap)  {

		  if (timeoutPassed_) return true;

          if ((warmup_absgap_ > 0.0 && current_absgap <= warmup_absgap_) ||
              (warmup_relgap_ > 0.0 && current_relgap <= warmup_relgap_))
          {
              timeoutPassed_ = true;
              lastImprovement_ = std::chrono::steady_clock::now();
              if (log_) fmt::print("    MP Plateau: Warmum gap reached.\n");
              return true;
          }
          // if condition on gap is not reached, check time
          auto now = std::chrono::steady_clock::now();
          if (std::chrono::duration<double>(now - lastImprovement_).count() >= warmup_time_)
          {
              timeoutPassed_ = true;
              lastImprovement_ = std::chrono::steady_clock::now();
              if (log_) fmt::print("    MP Plateau: Warmup time reached.\n");
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
          constexpr double kInf = std::numeric_limits<double>::infinity();
          if (!IsWarmupDone(kInf, kInf))
              return false;
          auto now = std::chrono::steady_clock::now();
          auto diff = std::chrono::duration<double>(now - lastImprovement_).count();
          auto stop = diff >= plateau_time_;

          if (stop) {
              fmt::print("    MP Plateau: stopping after {:.1f}s without sufficient improvement (limit={:.1f}s)\n",
                  diff, plateau_time_);
          }
          return stop;
      }

      /// Print a one-line status snapshot; called (if mip:plateaulog=1)
      /// on every ReportIncumbent()/ReportBound() call, i.e. on every
      /// native callback invocation that reaches the plateau logic.
      void LogStatus(const char* kind, double val, double baseline,
          bool progress, bool stop, 
          double sinceProgress) const {
          auto now = std::chrono::steady_clock::now();
          fmt::print(
              "    MP Plateau [{}]: current={:.6g} previous={:.6g} progress={} "
              "elapsed={:.1f}s (limit={:.1f}s){}\n",
              kind, val, baseline, progress ? "yes" : "no", 
              sinceProgress, plateau_time_,
              stop ? " -- STOPPING (plateau reached)" : "");
          std::fflush(stdout);
      }


      void SetCurrentObjective(int nobj) {
          Reset();
          if (log_)
          {
              fmt::print("   MP Plateau: new pass {}, resetting\n",
                  nobj + 1);
              std::fflush(stdout);
          }
      }
    bool ReportIncumbent(double obj) {
        if (!active_) return false;

        // Cannot assume I know the new MIP gap here
        if (!IsWarmupDone(std::numeric_limits<double>::infinity(), 1))
            return false;

        if (!haveIncumbent_) {
            haveIncumbent_ = true;
            bestObj_ = obj;
            auto now = std::chrono::steady_clock::now();
            double sinceProgress =
                std::chrono::duration<double>(now - lastImprovement_).count();
            lastImprovement_ = now  ;
            bool stop = CheckTimeout();
            if (log_)
                LogStatus("incumbent", obj, bestObj_, true, stop, sinceProgress);
            return stop;
        }
        return Improved("incumbent", bestObj_, obj, abstol_, reltol_);
    }


    
  bool ReportGap(double absgap, double relgap) {
      if (!active_) return false;
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
          bool stop = CheckTimeout();
          if (log_) {
              LogStatus("absmipgap", absgap, bestAbsGap_, true, stop, sinceProgress);
              LogStatus("relmipgap", relgap, bestRelGap_, true, stop, sinceProgress);
          }
          return stop;
      }

      // Compare only if abs or gap limit given
	  // Any of the two works as a sufficient condition for improvement
      constexpr double kInf = std::numeric_limits<double>::infinity();
      if (absmipgap_tol_ > 0.0)
          Improved("absmipgap", bestAbsGap_, absgap, absmipgap_tol_, kInf);
      if (relmipgap_tol_ > 0.0)
          Improved("relmipgap", bestRelGap_, relgap, relmipgap_tol_, kInf);

      return CheckTimeout();
    }
  };
  PlateauState plateauState_;
protected:
  const Options& GetMIPOptions() const { return mipStoredOptions_; }
  Options& GetMIPOptions() { return mipStoredOptions_; }

  int lazy_user_cuts() const {
    return IMPL_HAS_STD_FEATURE(LAZY_USER_CUTS) ?
                GetMIPOptions().lazy_user_cuts_ : 0;
  }
  /// Whether we need to mark .lazy>0 cuts as 'lazy'
  bool lazy_cuts() const { return 1 & lazy_user_cuts(); }
  /// Whether we need to mark .lazy<0 cuts as 'user'
  bool user_cuts() const { return 2 & lazy_user_cuts(); }

  int basis() const
  { return IMPL_HAS_STD_FEATURE(BASIS) ? GetMIPOptions().basis_ : 0; }
  bool need_basis_in() const { return 1 & basis(); }
  bool need_basis_out() const {
    return IsMIP() ?
          (need_fixed_MIP()) :          // assume the solver did it
          (2 & basis());
  }

  int warmstart() const
  { return IMPL_HAS_STD_FEATURE(WARMSTART) ? GetMIPOptions().warmstart_ : 0; }

  int priorities() const {
    return IMPL_HAS_STD_FEATURE(VAR_PRIORITIES) ?
          GetMIPOptions().importPriorities_ : 0;
  }

  int rays() const
  { return IMPL_HAS_STD_FEATURE(RAYS) ? GetMIPOptions().rays_ : 0; }
  bool need_ray_primal() const { return 1 & rays(); }
  bool need_ray_dual() const { return 2 & rays(); }

  int sensitivity() const {
      return IMPL_HAS_STD_FEATURE(SENSITIVITY_ANALYSIS) ?
                  GetMIPOptions().solnSens_ : 0;
  }

  /// Whether need duals/basis/sens from MIP
  /// Need at least duals when this option is on
  int need_fixed_MIP() const {
    return IMPL_HAS_STD_FEATURE( FIX_MODEL ) ?
          GetMIPOptions().fixModel_ : 0;
  }

  double plateau_time() const
  { return IMPL_HAS_STD_FEATURE(PLATEAU_STOP) ? plateauState_.plateau_time_ : 0.0; }

  /// Whether the plateau stopping logic should be armed for this solve
  bool plateau_active() const
  { return plateau_time() > 0.0; }


  void plateau_set_current_objective(int objn, int nobjs) {
	  // When setting the current objective in a native solve,
	  // also update the per-pass overrides from the options
	  // NOTE: check pass vs objective number for objectives with
	  // equal priority
	  // auto passes = pSetter->GetPassesWithOptions();

	  plateauState_.SetCurrentObjective(objn);
	  // If nobjs==-1, then we are in a MO-emulator solve, so the options
	  // are set while preparing the iteration
	  if (nobjs == -1) return;

	  auto pSetter = this->GetObjOptionSetter();
      double v;
      int iv;
      int pass = objn;

      static const std::unordered_map<std::string, double PlateauState::*> dblOptions = {
          {"plateautime", &PlateauState::plateau_time_},
          {"plateauabstol", &PlateauState::abstol_},
          {"plateaureltol", &PlateauState::reltol_},
          {"plateauwarmup", &PlateauState::warmup_time_},
          {"plateauwarmuprelgap", &PlateauState::warmup_relgap_},
          {"plateauwarmupabsgap", &PlateauState::warmup_absgap_},
          {"plateaurelgaptol", &PlateauState::relmipgap_tol_},
          {"plateauabsgaptol", &PlateauState::absmipgap_tol_}
      };
      static const std::unordered_map<std::string, int PlateauState::*> intOptions = {
         {"plateaulog", &PlateauState::log_}
      };
      for (const auto& [optName, member] : dblOptions) {
          if (pSetter->GetPassOptionValueDbl(pass, optName.c_str(), v))
              plateauState_.*member = v;
      }
      for (const auto& [optName, member] : intOptions)
          if (pSetter->GetPassOptionValueInt(pass, optName.c_str(), iv)) {
              plateauState_.*member = iv;
      }
  }
  /// Whether to print plateau status on every callback report
  bool plateau_log() const
  {  return 0!=plateauState_.log_; }

  /// To be called by the driver's native incumbent callback with the
  /// new incumbent objective value.
  /// Return true if the solve should be terminated now (plateau reached).
  bool ReportIncumbentForPlateau(double objVal) {
      return plateauState_.ReportIncumbent(objVal);
  }

  /// To be called by the driver's native callback with the current
  /// (relative) MIP gap. Uses mip:plateauabsgaptol and mip:plateaurelgaptol as 
  /// tolerances for the gap value(s)
  /// Return true if the solve should be terminated now (plateau reached).
  bool ReportGapForPlateau(double absgap, double relgap) {
      return plateauState_.ReportGap(absgap, relgap);
  }

  /// Can be called in auxiliary callbacks to increase the granularity
  /// of the plateau stopping logic. Returns true if the solve should be
  /// stopped (only checks if the time since the last sufficient improvement 
  /// has exceeded mip:plateautime).
  bool CheckTimeoutForPlateau() {
	  return plateauState_.CheckTimeout();
  }


public:
  void InitStandardOptions() override {
    BaseBackend::InitStandardOptions();
    InitMIPOptions();
  }

  /// Arm the plateau-stop callback(s), if requested and supported,
  /// alongside the standard timer/interrupter setup.
  void SetupTimerAndInterrupter() override {
    BaseBackend::SetupTimerAndInterrupter();
  }
  void SetupPlateau() override {
      if (plateau_active()) {
          SetupPlateauCallbacks();
          if (this->GetMM().IsMOEmulationOn())
            // Emulated MO: this pass's per-objective options are applied
            // separately, while preparing the iteration (SetMultiObjectiveOptions()).
            plateau_set_current_objective(0, -1);
          else
            // Native multiobj (or plain single-objective solve): there is
            // no separate per-iteration setup call, so pass 0's options
            // (if overridden) must be applied right away. Any later
            // passes are picked up by the driver's native per-objective
            // callback invoking plateau_set_current_objective() itself.
            plateau_set_current_objective(0, 0);
      }
  }
  using BaseBackend::AddStoredOption;

  using BaseBackend::debug_mode;


  ////////////////////////////////////////////////////////////////
protected:
  const mp::OptionValueInfo values_01_noyes_0default_[2] = {
    {     "0", "No (default)", 0 },
    {     "1", "Yes.", 1}
  };

  const mp::OptionValueInfo values_01_noyes_1default_[2] = {
    {     "0", "No", 0 },
    {     "1", "Yes (default)", 1}
  };

  const mp::OptionValueInfo values_lpwarmstart_[4] = {
      {     "-1", "Default (equivalent to 2 for PDHG, to 1 otherwise)", -1 },
      {     "0", "Ignore any warm start information (generally).", 0 },
      {     "1", "Use warm start information to solve the original, unpresolved problem.", 1},
      {     "2", "If presolve is enabled, use warm start to solve the presolved problem. "
       "Otherwise, setting 2 prioritizes start vectors (primal/dual), while "
       "setting 1 prioritizes basis statuses.", 2 }
  };

  const mp::OptionValueInfo values_autonoyes_[3] = {
    {     "-1", "Automatic choice (default)", 0 },
    {     "0", "No", 0 },
    {     "1", "Yes.", 1}
  };

  const mp::OptionValueInfo values_autonomodaggr_[4] = {
    {     "-1", "Automatic choice (default)", 0 },
    {     "0", "No", 0 },
    {     "1", "Yes, moderate", 1},
    {     "2", "Yes, aggressive.", 2}
  };

  const mp::OptionValueInfo values_autonoconsaggr_[4] = {
      { "-1", "Automatic choice (default)", -1},
      { "0", "No", 0},
      { "1", "Conservative", 1},
      { "2", "Aggressive.", 2}
  };

  const mp::OptionValueInfo values_basis_[4] = {
    {     "0", "No", 0 },
    {     "1", "Use incoming basis (if provided)", 1},
    {     "2", "Return final basis", 2},
    {     "3", "Both (1 + 2 = default)", 3}
  };

  const mp::OptionValueInfo values_warmstart_[4] = {
      {     "0", "No", 0 },
      {     "1", "Yes (for LP: if there is no incoming alg:basis)", 1},
      {     "2", "Yes (for LP: omitting the incoming alg:basis, if any)", 2},
      {     "3", "Yes (for LP: together with the incoming alg:basis, if any; default).", 3}
  };

  const mp::OptionValueInfo values_rays_[4] = {
    {     "0", "Neither", 0 },
    {     "1", "Just .unbdd", 1},
    {     "2", "Just .dunbdd", 2},
    {     "3", "Both (default)", 3}
  };

  ////////////////////////////////////////////////////////////////
  virtual void InitMIPOptions() {
      if (IMPL_HAS_STD_FEATURE( LAZY_USER_CUTS ))
        AddStoredOption("mip:lazy lazy",
          "Whether to recognize suffix .lazy on constraints: "
          "sum of\n"
          "\n"
          "|  1 - Accept .lazy>0 values (true lazy constraints, if supported)\n"
          "|  2 - Accept .lazy<0 values (user cuts, if supported)\n"
          "\n"
          "Default lazy = 3 ==> accept both.",
          GetMIPOptions().lazy_user_cuts_);

    if (IMPL_HAS_STD_FEATURE( BASIS ))
      AddStoredOption("alg:basis basis",
                      "Whether to use and/or return a basis for LP models "
                      "(variable/constraint suffixes .(s)status):\n"
                      "\n.. value-table::\n"
                      "\n"
                      "See alg:start for interaction with the LP warmstart.\n"
                      "\n"
                      "See also mip:basis and qcp:dual (for some solvers).",
                      GetMIPOptions().basis_, values_basis_);

    if (IMPL_HAS_STD_FEATURE( WARMSTART ))
      AddStoredOption("alg:start warmstart",
                      "Whether to use incoming primal (and dual, for LP) variable values "
                      "in a warmstart:\n "
                      "\n.. value-table::",
                      GetMIPOptions().warmstart_, values_warmstart_);

    if (IMPL_HAS_STD_FEATURE( VAR_PRIORITIES ))
      AddStoredOption("mip:priorities priorities",
        "0/1*: Whether to read the branch and bound priorities from the"
        " .priority suffix.",
        GetMIPOptions().importPriorities_);


    if (IMPL_HAS_STD_FEATURE( RAYS ))
      AddStoredOption("alg:rays rays",
                      "Whether to return suffix .unbdd (unbounded ray) "
                      "if the objective is unbounded "
                      "or suffix .dunbdd (Farkas dual) if the constraints "
                      "are infeasible:\n"
                      "\n.. value-table::\n",
                      GetMIPOptions().rays_, values_rays_);

    if (IMPL_HAS_STD_FEATURE( IIS ))
      AddStoredOption("iis:find iisfind iis alg:iisfind",
                      "Whether to find and export an IIS. "
                      "Default = 0 (don't export).",
                      GetMIPOptions().exportIIS_);

    if (IMPL_HAS_STD_FEATURE( RETURN_MIP_GAP ))
      AddStoredOption("mip:return_gap return_mipgap",
        "Whether to return mipgap suffixes or include mipgap values "
    "(|objectve - .bestbound|) in the solve_message:  sum of\n"
    "\n"
    "| 1 - Return .relmipgap suffix (relative to |obj|)\n"
    "| 2 - Return .absmipgap suffix (absolute mipgap)\n"
    "| 4 - Suppress mipgap values in solve_message.\n"
    "\n"
    "Default = 0.  The suffixes are on the objective and problem. "
    "Returned suffix values are +Infinity if no integer-feasible "
    "solution has been found, in which case no mipgap values are "
    "reported in the solve_message.",
        GetMIPOptions().returnMipGap_);

    if (IMPL_HAS_STD_FEATURE( RETURN_BEST_DUAL_BOUND ))
      AddStoredOption("mip:bestbound bestbound return_bound",
        "Whether to return suffix .bestbound for the "
        "best known MIP dual bound on the objective value:\n"
        "\n.. value-table::\n"
        "The suffix is on the objective and problem and is -Infinity "
        "for minimization problems and +Infinity for maximization "
        "problems if there are no integer variables or if a dual bound "
        "is not available.",
          GetMIPOptions().returnBestDualBound_, values_01_noyes_0default_);

    if (IMPL_HAS_STD_FEATURE( SENSITIVITY_ANALYSIS ))
      AddStoredOption("alg:sens sens solnsens sensitivity",
                      "Whether to return suffixes for solution sensitivities, i.e., "
                      "ranges of values for which the optimal basis remains optimal "
                      "(note that the variable and objective values can change):\n"
                      "\n"
                      "|  0 - No (default)\n"
                      "|  1 - Yes:  suffixes returned on variables are\n"
                      "|    .sensobjlo = smallest objective coefficients\n"
                      "|    .down      = same as .sensobjlo\n"
                      "|    .sensobj   = current objective coefficients\n"
                      "|    .current   = same as .sensobj\n"
                      "|    .sensobjhi = greatest objective coefficients\n"
                      "|    .up        = same as .sensobjhi\n"
                      "|    .senslblo  = smallest variable lower bounds\n"
                      // "|    .senslb   = current variable lower bounds\n"
                      "|    .senslbhi  = greatest variable lower bounds\n"
                      "|    .sensublo  = smallest variable upper bounds\n"
                      // "|    .sensub   = current variable upper bounds\n"
                      "|    .sensubhi  = greatest variable upper bounds;\n\n"
                      " suffixes for all constraints are\n"
                      "|    .senslblo = smallest constraint lower bounds\n"
                      // "|    .senslb   = current constraint lower bounds\n"
                      "|    .senslbhi = greatest constraint lower bounds\n"
                      "|    .sensublo = smallest constraint upper bounds\n"
                      // "|    .sensub   = current constraint upper bounds\n"
                      "|    .sensubhi = greatest constraint upper bounds;\n\n"
                      " suffixes for one-sided constraints only:\n"
                      "|    .sensrhslo = smallest right-hand side values\n"
                      "|    .down      = same as .sensrhslo\n"
                      // "|    .sensrhs   = current right-hand side values\n"
                      "|    .sensrhshi = greatest right-hand side values.\n"
                      "|    .up        = same as .sensrhshi.\n"
                      "\n"
                      "The suffixes correspond to the AMPL solver model, "
                      "command 'solexpand'. For easiest interpretation, "
                      "disable AMPL presolve, 'option presolve 0;'"
                      ,
                    GetMIPOptions().solnSens_);

    if (IMPL_HAS_STD_FEATURE( FIX_MODEL ))
      AddStoredOption("mip:basis fixmodel mip:fix",
                      "Whether to compute duals / basis / sensitivity for MIP models:\n"
                      "\n.. value-table::\n",
                    GetMIPOptions().fixModel_, values_01_noyes_0default_);

    if (IMPL_HAS_STD_FEATURE( PLATEAU_STOP )) {
      AddStoredOption("mip:plateautime plateautime",
        "Stop the MIP search if the incumbent objective (and, if "
        "mip:plateauabsgaptol/mip:plateaurelgaptol are set, the MIP gap) "
        "has not improved by at least mip:plateauabstol or "
        "mip:plateaureltol for this many seconds. Default 0 (disabled).",
        plateauState_.plateau_time_);

      AddStoredOption("mip:plateauabstol plateauabstol",
        "Minimum absolute objective improvement to reset the "
        "mip:plateautime timer. Default 0 (any improvement resets the timer).",
        plateauState_.abstol_);

      AddStoredOption("mip:plateaureltol plateaureltol",
        "Minimum relative objective improvement, as a fraction of the current "
        "value, to reset the mip:plateautime timer. Default 0.",
        plateauState_.reltol_);

      AddStoredOption("mip:plateauwarmup plateauwarmup",
        "Grace period (in seconds) after the solve starts before "
        "mip:plateautime is checked. Default 0.",
        plateauState_.warmup_time_);

      AddStoredOption("mip:plateaulog plateaulog",
          "Whether to print the current plateau status. Default 0 (silent).",
          plateauState_.log_);
    }

    if (IMPL_HAS_STD_FEATURE( PLATEAU_STOP_BOUND )) {
        // Only meaningful for solvers that also report gap: an
        // incumbent-only driver never has a real gap value to offer
        // these checks, so they'd never do anything for it.
        AddStoredOption("mip:plateauwarmuprelgap plateauwarmuprelgap",
            "Relative MIP gap to be reached before mip:plateautime is "
            "checked. Default 0.",
            plateauState_.warmup_relgap_);

        AddStoredOption("mip:plateauwarmupabsgap plateauwarmupabsgap",
            "Absolute MIP gap to be reached before mip:plateautime is "
            "checked. Default 0.",
            plateauState_.warmup_absgap_);

        AddStoredOption("mip:plateaurelgaptol plateaurelgaptol",
            "If set (>0), also track the reported relative MIP gap as progress "
            "for mip:plateautime: the plateau timer resets whenever the relative "
            "gap shrinks by at least this amount. Default 0 (disabled).",
            plateauState_.relmipgap_tol_);

        AddStoredOption("mip:plateauabsgaptol plateauabsgaptol",
            "If set (>0), also track the reported absolute MIP gap as progress "
            "for mip:plateautime: the plateau timer resets whenever the absolute "
            "gap shrinks by at least this amount. Default 0 (disabled).",
            plateauState_.absmipgap_tol_);
    }
  }



  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////// STANDARD MIP SUFFIXES ///////////////////////////
  //////////////////////////////////////////////////////////////////////////////
private:

  const SuffixDef<int> suf_varstatus = { "sstatus", suf::VAR | suf::OUTPUT };
  const SuffixDef<int> suf_constatus = { "sstatus", suf::CON | suf::OUTPUT };
  /// Testing API
  /// Output suffix values to check they were read correctly
  const SuffixDef<int> suf_testvarstatus = { "test_sstatus", suf::VAR | suf::OUTPUT };
  const SuffixDef<int> suf_testconstatus = { "test_sstatus", suf::CON | suf::OUTPUT };

  /// Testing API
  /// Output primal/dual initials to check they were read correctly
  const SuffixDef<double> suf_testvarini = { "test_ini_pri", suf::VAR | suf::OUTPUT };
  const SuffixDef<double> suf_testconini = { "test_ini_dua", suf::CON | suf::OUTPUT };

  /// Testing API
  /// Output MIP initials to check they were read correctly
  const SuffixDef<double> suf_testMIPini = { "test_ini_mip", suf::VAR | suf::OUTPUT };

  const SuffixDef<int> suf_varpriority = { "priority", suf::VAR | suf::INPUT };

  const SuffixDef<double> suf_unbdd  = { "unbdd",   suf::VAR | suf::OUTPUT };
  const SuffixDef<double> suf_dunbdd = { "dunbdd",  suf::CON | suf::OUTPUT };

  const SuffixTable iis_table =
      "\n"
      "0\tnon\tnot in the iis\n"
      "1\tlow\tlower bound in the iis\n"
      "2\tfix\tboth bounds in the iis\n"
      "3\tupp\tupper bound in the iis\n"
      "4\tmem\tmember\n"
      "5\tpmem\tpossible member\n"
      "6\tplow\tpossibly lower bound\n"
      "7\tpupp\tpossibly upper bound\n"
      "8\tbug\n"
      "9\tintvar\tinteger variable\n"
      "10\tsemi\tsemi-continuous or semi-integer variable\n"
      ;
  const SuffixDef<int> sufIISCon = { "iis", suf::CON | suf::OUTPUT, iis_table };
  const SuffixDef<int> sufIISVar = { "iis", suf::VAR | suf::OUTPUT, iis_table };

  const SuffixDef<double> sufRelMipGapObj = { "relmipgap", suf::OBJ | suf::OUTPUT };
  const SuffixDef<double> sufRelMipGapProb = { "relmipgap", suf::PROBLEM | suf::OUTPUT };
  const SuffixDef<double> sufAbsMipGapObj = { "absmipgap", suf::OBJ | suf::OUTPUT };
  const SuffixDef<double> sufAbsMipGapProb = { "absmipgap", suf::PROBLEM | suf::OUTPUT };
  const SuffixDef<double> sufBestBoundObj = { "bestbound", suf::OBJ | suf::OUTPUT };
  const SuffixDef<double> sufBestBoundProb = { "bestbound", suf::PROBLEM | suf::OUTPUT };

};

}  // namespace mp

#endif  // MIPBACKEND_H_

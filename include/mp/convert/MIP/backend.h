/*
 Abstract MIP solver backend wrapper.

 Copyright (C) 2020 AMPL Optimization Inc

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

#include "mp/common.h"
#include "mp/convert/backend.h"

namespace mp {

/// MIP backend wrapper.
/// The MIP wrapper provides common functionality relative to MIP solvers;
/// it implements the common suffixes and the logic shared across all MIP
/// solvers
template <class Impl>
class MIPBackend :
  public BasicBackend<Impl>
{
  using BaseBackend = BasicBackend<Impl>;
  using BaseBackend::debug_mode;

public:
  // Properties
  bool IsMIP() const { UNSUPPORTED( "IsMIP()" ); }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  /////////////// Most are disabled by default ///////////////
  //// To enable, declare ALLOW_STD_FEATURE( name, true ) ////
  /////////// and implement the relevant methods /////////////
  ////////////////////////////////////////////////////////////
  /**
  * Get/Set AMPL var/con statii
  **/
  DEFINE_STD_FEATURE( BASIS )
  ALLOW_STD_FEATURE( BASIS, false )
  std::vector<int> VarStatii() { return {}; }
  void VarStatii(ArrayRef<int> )
  { UNSUPPORTED("MIPBackend::VarStatii"); }
  std::vector<int> ConStatii() { return {}; }
  void ConStatii(ArrayRef<int> )
  { UNSUPPORTED("MIPBackend::ConStatii"); }
  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  DEFINE_STD_FEATURE( WARMSTART )
  ALLOW_STD_FEATURE( WARMSTART, false )
  void InputPrimalDualStart(ArrayRef<double> x0,
                       ArrayRef<double> pi0)
  { UNSUPPORTED("MIPBackend::InputPrimalDualStart"); }
  /**
  * Specifically, MIP warm start
  **/
  DEFINE_STD_FEATURE( MIPSTART )
  ALLOW_STD_FEATURE( MIPSTART, false )
  void AddMIPStart(ArrayRef<double> x0)
  { UNSUPPORTED("MIPBackend::AddMIPStart"); }
  /**
  * Set branch and bound priority
  **/
  DEFINE_STD_FEATURE( VAR_PRIORITIES )
  ALLOW_STD_FEATURE( VAR_PRIORITIES, false )
  void VarPriorities(ArrayRef<int>)
  { UNSUPPORTED("MIPBackend::VarPriorities"); }
  /**
  * Obtain unbounded/inf rays
  **/
  DEFINE_STD_FEATURE( RAYS )
  ALLOW_STD_FEATURE( RAYS, false )
  std::vector<double> Ray() { return {}; }
  std::vector<double> DRay() { return {}; }
  /**
  * Compute the IIS and obtain relevant values
  **/
  DEFINE_STD_FEATURE( IIS )
  ALLOW_STD_FEATURE( IIS, false )
  void ComputeIIS() {}
  /// Elements correspond to IISStatus
  std::vector<int> ConsIIS() { return {}; }
  std::vector<int> VarsIIS() { return {}; }
  /**
  * Get MIP Gap
  **/
  DEFINE_STD_FEATURE( RETURN_MIP_GAP )
  ALLOW_STD_FEATURE( RETURN_MIP_GAP, false )
  double MIPGap() const { return MP_DISPATCH( Infinity() ); }
  double MIPGapAbs() const {
    return std::fabs(
          MP_CONST_DISPATCH(ObjectiveValue()) -
          MP_CONST_DISPATCH(BestDualBound()) );
  }
  /**
  * Get MIP dual bound
  **/
  DEFINE_STD_FEATURE( RETURN_BEST_DUAL_BOUND )
  ALLOW_STD_FEATURE( RETURN_BEST_DUAL_BOUND, false )
  double BestDualBound() const
  { UNSUPPORTED("BestDualBound()"); return 0; }


  ////////////////////////////////////////////////////////////////////////////
  /////////////////////// MIP specific derived calculations //////////////////
  ////////////////////////////////////////////////////////////////////////////

  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         INPUT         //////////////////////////
  using BaseBackend::ReadSuffix;
  using BaseBackend::ReportSuffix;

  void InputStdExtras() {
    BasicBackend<Impl>::InputStdExtras();
    InputMIPExtras();
  }

  void InputMIPExtras() {
    MP_DISPATCH( InputStartValues() );
    if (priorities())
      MP_DISPATCH( VarPriorities( ReadSuffix(suf_varpriority) ) );
  }

  void InputStartValues() {
    if (MP_DISPATCH( IsMIP() )) {
      if (warmstart() &&
          IMPL_HAS_STD_FEATURE( MIPSTART )) {
        MP_DISPATCH( AddMIPStart(
                       MP_DISPATCH( InitialValues() ) ) );
        if (debug_mode()) {                    // Report received initials
          ReportSuffix(suf_testMIPini,         // Should we check that
                       MP_DISPATCH( InitialValues() ));
        }                                      // Impl uses them?
      }
    } else
      MP_DISPATCH( InputPrimalDualStartOrBasis() );
  }

  void InputPrimalDualStartOrBasis() {
    bool useBasis = need_basis_in();
    ArrayRef<int> varstt, constt;
    if (useBasis) {
      varstt = ReadSuffix(suf_varstatus);
      constt = ReadSuffix(suf_constatus);
      useBasis = varstt.size() && constt.size();
    }
    auto X0 = MP_DISPATCH( InitialValues() );
    auto pi0 = MP_DISPATCH( InitialDualValues() );
    bool haveInis = X0.size() && pi0.size();
    if (haveInis && (
          2==warmstart() ||
          (1==warmstart() && !useBasis))) {
      MP_DISPATCH( InputPrimalDualStart(X0, pi0) );
      useBasis = false;
      if (debug_mode()) {                 // Report received initials
        ReportSuffix(suf_testvarini, X0); // Should we check that
        ReportSuffix(suf_testconini, pi0); // Impl uses them?
      }
    }
    if (useBasis) {
      MP_DISPATCH( VarStatii(varstt) );
      MP_DISPATCH( ConStatii(constt) );
      if (debug_mode()) {                    // Report received statuses
        ReportSuffix(suf_testvarstatus, varstt); // Should we check that
        ReportSuffix(suf_testconstatus, constt); // Impl uses them?
      }
    }
  }


  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         OUtPUT        //////////////////////////
  void ReportStandardSuffixes() {
    BaseBackend::ReportStandardSuffixes();
    ReportStandardMIPSuffixes();
  }

  void ReportStandardMIPSuffixes() {
    if (need_basis_out())
      MP_DISPATCH( ReportBasis() );
    MP_DISPATCH( ReportRays() );
    MP_DISPATCH( CalculateAndReportIIS() );
    MP_DISPATCH( CalculateAndReportMIPGap() );
    MP_DISPATCH( ReportBestDualBound() );
  }

  void ReportBasis() {
    ReportSuffix(suf_varstatus,
                              MP_DISPATCH( VarStatii() ));
    ReportSuffix(suf_constatus,
                              MP_DISPATCH( ConStatii() ));
  }

  void ReportRays() {
    if ( need_ray_primal() && MP_DISPATCH( IsProblemUnbounded() )) {
      ReportSuffix(suf_unbdd, MP_DISPATCH( Ray() ));
    }
    if ( need_ray_dual() && MP_DISPATCH( IsProblemInfeasible() )) {
      ReportSuffix(suf_dunbdd, MP_DISPATCH( DRay() ));
    }
  }

  void CalculateAndReportIIS() {
    if (MP_DISPATCH( IsProblemInfOrUnb() ) &&
        mipStoredOptions_.exportIIS_) {
      MP_DISPATCH( ComputeIIS() );

      ReportSuffix(sufIISCon, MP_DISPATCH(ConsIIS()));
      ReportSuffix(sufIISVar, MP_DISPATCH(VarsIIS()));
    }
  }

  void CalculateAndReportMIPGap() {
    if (0 < MP_DISPATCH(NumberOfObjectives()) ) {
      std::vector<double> dbl(1);
      if (1 & mipStoredOptions_.returnMipGap_) {
        dbl[0] = MP_DISPATCH( MIPGap() );
        ReportSuffix(sufRelMipGapObj, dbl);
        ReportSuffix(sufRelMipGapProb, dbl);
      }
      if (2 & mipStoredOptions_.returnMipGap_) {
        dbl[0] = MP_DISPATCH( MIPGapAbs() );
        ReportSuffix(sufAbsMipGapObj, dbl);
        ReportSuffix(sufAbsMipGapProb, dbl);
      }
    }
  }

  void ReportBestDualBound() {
    if (mipStoredOptions_.returnBestDualBound_ &&
        0 < MP_DISPATCH(NumberOfObjectives()) ) {
      std::vector<double> dbl(1, MP_DISPATCH( BestDualBound() ));
      ReportSuffix(sufBestBoundObj, dbl);
      ReportSuffix(sufBestBoundProb, dbl);
    }
  }


  ////////////////////////////////////////////////////////////
  /////////////////// MIP Backend options ////////////////////
  ////////////////////////////////////////////////////////////
private:
  struct Options {
    int basis_=3;
    int warmstart_=1;
    int importPriorities_=1;
    int rays_=3;
    int exportIIS_=0;
    int returnMipGap_=0;
    int returnBestDualBound_=0;
  };
  Options mipStoredOptions_;

protected:
  const Options& GetMIPOptions() const { return mipStoredOptions_; }

  int basis() const
  { return IMPL_HAS_STD_FEATURE(BASIS) ? GetMIPOptions().basis_ : 0; }
  bool need_basis_in() const { return 1 & basis(); }
  bool need_basis_out() const { return 2 & basis(); }

  int warmstart() const
  { return IMPL_HAS_STD_FEATURE(WARMSTART) ? GetMIPOptions().warmstart_ : 0; }

  int priorities() const {
    return IMPL_HAS_STD_FEATURE(VAR_PRIORITIES) ?
          mipStoredOptions_.importPriorities_ : 0;
  }

  int rays() const
  { return IMPL_HAS_STD_FEATURE(RAYS) ? GetMIPOptions().rays_ : 0; }
  bool need_ray_primal() const { return 1 & rays(); }
  bool need_ray_dual() const { return 2 & rays(); }


public:
  void InitStandardOptions() {
    BaseBackend::InitStandardOptions();
    MP_DISPATCH( InitMIPOptions() );
  }

  using BaseBackend::AddStoredOption;


  ////////////////////////////////////////////////////////////////
protected:
  const mp::OptionValueInfo values_01_noyes_0default_[2] = {
      {     "0", "no (default)", 0 },
      {     "1", "yes", 1}
  };

  const mp::OptionValueInfo values_basis_[4] = {
      {     "0", "no", 0 },
      {     "1", "use incoming basis (if provided)", 1},
      {     "2", "return final basis", 2},
      {     "3", "both (1 + 2 = default).", 3}
  };

  const mp::OptionValueInfo values_warmstart_[3] = {
      {     "0", "no", 0 },
      {     "1", "yes (for LP: if there is no incoming mip:basis) (default)", 1},
      {     "2", "yes (for LP: ignoring the incoming mip:basis, if any)", 2}
  };

  const mp::OptionValueInfo values_rays_[4] = {
      {     "0", "neither", 0 },
      {     "1", "just .unbdd", 1},
      {     "2", "just .dunbdd", 2},
      {     "3", "both (default)", 3}
  };

  ////////////////////////////////////////////////////////////////
  void InitMIPOptions() {
    if (IMPL_HAS_STD_FEATURE( BASIS ))
      AddStoredOption("mip:basis basis",
                      "Whether to use or return a basis:\n "
                      "\n.. value-table::\n"
                      "Note that if you provide a valid starting extreme point, "
                      "either through primal/dual status, or through warmstart, "
                      "then Gurobi LP presolve will be disabled. For models where "
                      "presolve greatly reduces the problem size, "
                      "this might hurt performance.",
                      mipStoredOptions_.basis_, values_basis_);

    if (IMPL_HAS_STD_FEATURE( WARMSTART ))
      AddStoredOption("alg:start warmstart",
                      "Whether to use incoming primal (and dual, for LP) variable values "
                      "in a warmstart:\n "
                      "\n.. value-table::\n"
                      "Note that for LP, mip:basis is usually more efficient. "
                      "MIP-specific options can be accessible via mip:start.",
                      mipStoredOptions_.warmstart_, values_warmstart_);

    if (IMPL_HAS_STD_FEATURE(VAR_PRIORITIES))
      AddStoredOption("mip:priorities priorities",
        "0/1*: Whether to read the branch and bound priorities from the"
        " .priority suffix.",
        mipStoredOptions_.importPriorities_);


    if (IMPL_HAS_STD_FEATURE( RAYS ))
      AddStoredOption("alg:rays rays",
                      "Whether to return suffix .unbdd if the objective is unbounded "
                      "or suffix .dunbdd if the constraints are infeasible:\n"
                      "\n.. value-table::\n",
                      mipStoredOptions_.rays_, values_rays_);

    if (IMPL_HAS_STD_FEATURE( IIS ))
      AddStoredOption("alg:iisfind iisfind",
                      "Whether to find and export the IIS. "
                      "Default = 0 (don't export).",
                      mipStoredOptions_.exportIIS_);

    if (IMPL_HAS_STD_FEATURE( RETURN_MIP_GAP ))
      AddStoredOption("mip:return_gap return_mipgap",
        "Whether to return mipgap suffixes or include mipgap values "
    "(|objectve - .bestbound|) in the solve_message:  sum of\n"
    "\n"
    "| 1 - return .relmipgap suffix (relative to |obj|)\n"
    "| 2 - return .absmipgap suffix (absolute mipgap)\n"
    "| 4 - suppress mipgap values in solve_message.\n"
    "\n"
    "Default = 0.  The suffixes are on the objective and problem. "
    "Returned suffix values are +Infinity if no integer-feasible "
    "solution has been found, in which case no mipgap values are "
    "reported in the solve_message.",
        mipStoredOptions_.returnMipGap_);

    if (IMPL_HAS_STD_FEATURE( RETURN_BEST_DUAL_BOUND ))
      AddStoredOption("mip:bestbound bestbound returnbound",
        "Whether to return suffix .bestbound for the "
        "best known MIP dual bound on the objective value:\n"
        "\n.. value-table::\n"
        "The suffix is on the objective and problem and is -Infinity "
        "for minimization problems and +Infinity for maximization "
        "problems if there are no integer variables or if a dual bound "
        "is not available.",
          mipStoredOptions_.returnBestDualBound_, values_01_noyes_0default_);

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
      "0	non	not in the iis\n"
      "1	low	at lower bound\n"
      "2	fix	fixed\n"
      "3	upp	at upper bound\n"
      "4	mem	member\n"
      "5	pmem	possible member\n"
      "6	plow	possibly at lower bound\n"
      "7	pupp	possibly at upper bound\n"
      "8	bug\n";
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

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

#include "mp/common.h"
#include "mp/backend.h"

namespace mp {

/// Basis status values of a solution (unpresolved)
struct SolutionBasis {
  /// Check if has both vars and cons' statuses
  operator bool() const { return varstt.size() && constt.size(); }
  /// Var and con statuses
  std::vector<int> varstt, constt;
};

/// IIS (unpresolved).
/// Elements correspond to IISStatus
struct IIS {
  /// Var and con IIS statuses
  std::vector<int> variis, coniis;
};

/// MIP backend wrapper
///
/// The MIP wrapper provides common functionality relative to MIP solvers;
/// it implements the common suffixes and the logic shared across all MIP
/// solvers
template <class Impl,
          class BaseBackend = Backend<Impl> >  ///< parameter for base class
class MIPBackend : public BaseBackend
{
public:
  /// Properties
  bool IsMIP() const override { return false; }
  virtual bool IsQP() const { return false; }
  virtual bool IsQCP() const { return false; }
  /// Always add MIP start if supported:
  /// Gurobi 9.1.2 solves non-convex Q(C)P as MIP
  /// But model attributes don't work before solve
  /// Count non-linear stuff in the model?
  virtual bool CanBeMIP() const { return true; }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  /////////////// Most are disabled by default ///////////////
  //// To enable, declare ALLOW_STD_FEATURE( name, true ) ////
  /////////// and implement the relevant methods /////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;
  /**
  * Set lazy/user cut attributes
  * Negative suffix values are "user cuts"
  * Check lazy_/user_cuts() to see which kinds are allowed
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
  /// Getter (unpresolved)
  virtual SolutionBasis GetBasis() { return {}; }
  /// Setter (unpresolved)
  virtual void SetBasis(SolutionBasis ) { UNSUPPORTED("MIPBackend::SetBasis"); }
  /**
  * General LP warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  DEFINE_STD_FEATURE( WARMSTART )
  ALLOW_STD_FEATURE( WARMSTART, false )
  virtual void AddPrimalDualStart(Solution )
  { UNSUPPORTED("MIPBackend::AddPrimalDualStart"); }
  /**
  * Specifically, MIP warm start
  **/
  DEFINE_STD_FEATURE( MIPSTART )
  ALLOW_STD_FEATURE( MIPSTART, false )
  virtual void AddMIPStart(ArrayRef<double> )
  { UNSUPPORTED("MIPBackend::AddMIPStart"); }
  /**
  * Set branch and bound priority
  **/
  DEFINE_STD_FEATURE( VAR_PRIORITIES )
  ALLOW_STD_FEATURE( VAR_PRIORITIES, false )
  virtual void VarPriorities(ArrayRef<int>)
  { UNSUPPORTED("MIPBackend::VarPriorities"); }
  /**
  * Obtain unbounded/inf rays
  **/
  DEFINE_STD_FEATURE( RAYS )
  ALLOW_STD_FEATURE( RAYS, false )
  virtual ArrayRef<double> Ray() { return {}; }
  virtual ArrayRef<double> DRay() { return {}; }
  /**
  * Compute the IIS and obtain relevant values
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
  virtual double MIPGap() { return MP_DISPATCH( Infinity() ); }
  virtual double MIPGapAbs() { return MP_DISPATCH( Infinity() ); }
  /**
  * Get MIP dual bound
  **/
  DEFINE_STD_FEATURE( RETURN_BEST_DUAL_BOUND )
  ALLOW_STD_FEATURE( RETURN_BEST_DUAL_BOUND, false )
  virtual double BestDualBound()
  { UNSUPPORTED("BestDualBound()"); return 0.0; }
  /**
  * Report sensitivity analysis suffixes
  **/
  DEFINE_STD_FEATURE( SENSITIVITY_ANALYSIS )
  ALLOW_STD_FEATURE( SENSITIVITY_ANALYSIS, false )
  virtual ArrayRef<double> Senslbhi() const { return {}; }
  virtual ArrayRef<double> Senslblo() const { return {}; }
  virtual ArrayRef<double> Sensobjhi() const { return {}; }
  virtual ArrayRef<double> Sensobjlo() const { return {}; }
  virtual ArrayRef<double> Sensrhshi() const { return {}; }
  virtual ArrayRef<double> Sensrhslo() const { return {}; }
  virtual ArrayRef<double> Sensubhi() const { return {}; }
  virtual ArrayRef<double> Sensublo() const { return {}; }
  /**
  * FixModel - duals, basis, and sensitivity for MIP
  * No API to overload,
  * Impl should check need_fixed_MIP()
  **/
  DEFINE_STD_FEATURE( FIX_MODEL )
  ALLOW_STD_FEATURE( FIX_MODEL, false )


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

  void InputStdExtras() override {
    BaseBackend::InputStdExtras();
    InputMIPExtras();
  }

  virtual void InputMIPExtras() {
    if (lazy_user_cuts())
      InputLazyUserCuts();
    InputStartValues();
    if (priorities())
      VarPriorities( ReadSuffix(suf_varpriority) );
  }

  virtual void InputLazyUserCuts() {
    auto sufLazyVal = ReadIntSuffix( {"lazy", suf::CON} );
    if (sufLazyVal)
      MarkLazyOrUserCuts(sufLazyVal);
  }

  /// Testing API
  virtual void ReportFirstLinearConstraintLazySuffix(int val) {
    if (debug_mode())
      ReportIntSuffix({"test_lin_constr_lazy", suf::PROBLEM}, {{val}});
  }

  /// Testing API
  virtual void ReportFirstQCLazySuffix(int val) {
    if (debug_mode())
      ReportIntSuffix({"test_quad_constr_lazy", suf::PROBLEM}, {{val}});
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
    sol0.dual = this->InitialDualValues();
    bool haveInis = sol0.primal.size() && sol0.dual.size();
    if (haveInis && (
          2==warmstart() ||
          (1==warmstart() && !useBasis))) {
      AddPrimalDualStart(sol0);
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
    if (warmstart() &&
        IMPL_HAS_STD_FEATURE( MIPSTART )) {
      AddMIPStart( this->InitialValues() );
      if (debug_mode()) {                    // Report received initials
        ReportSuffix(suf_testMIPini,         // Should we check that
                     this->InitialValues());
      }                                      // Impl uses them?
    }
  }


  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         OUtPUT        //////////////////////////
  virtual void ReportStandardSuffixes() {
    BaseBackend::ReportStandardSuffixes();
    ReportStandardMIPSuffixes();
  }

  virtual void ReportStandardMIPSuffixes() {
    if (need_basis_out())
      ReportBasis();
    ReportRays();
    CalculateAndReportIIS();
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
    if (( this->IsProblemInfOrUnb() ||
               this->IsProblemIndiffInfOrUnb() ) &&
        GetMIPOptions().exportIIS_) {
      ComputeIIS();

      auto iis = GetIIS();

      ReportSuffix(sufIISCon, iis.coniis);
      ReportSuffix(sufIISVar, iis.variis);
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
  }

  virtual void ReportBestDualBound() {
    if (GetMIPOptions().returnBestDualBound_) {
      std::vector<double> dbl(1, MP_DISPATCH( BestDualBound() ));
      ReportSuffix(sufBestBoundObj, dbl);
      ReportSuffix(sufBestBoundProb, dbl);
    }
  }

  virtual void ReportSensitivity() {
    ReportSuffix( {"senslbhi", suf::Kind::VAR},
                     MP_DISPATCH( Senslbhi() ) );
    ReportSuffix( {"senslblo", suf::Kind::VAR},
                     MP_DISPATCH( Senslblo() ) );
    ReportSuffix( {"sensobjhi", suf::Kind::VAR},
                     MP_DISPATCH( Sensobjhi() ) );
    ReportSuffix( {"sensobjlo", suf::Kind::VAR},
                     MP_DISPATCH( Sensobjlo() ) );
    ReportSuffix( {"sensrhshi", suf::Kind::CON},
                     MP_DISPATCH( Sensrhshi() ) );
    ReportSuffix( {"sensrhslo", suf::Kind::CON},
                     MP_DISPATCH( Sensrhslo() ) );
    ReportSuffix( {"sensubhi", suf::Kind::VAR},
                     MP_DISPATCH( Sensubhi() ) );
    ReportSuffix( {"sensublo", suf::Kind::VAR},
                     MP_DISPATCH( Sensublo() ) );
  }

  ////////////////////////////////////////////////////////////
  /////////////////// MIP Backend options ////////////////////
  ////////////////////////////////////////////////////////////
private:
  struct Options {
    int lazy_user_cuts_ = 3;
    int basis_=3;
    int warmstart_=1;
    int importPriorities_=1;
    int rays_=3;
    int exportIIS_=0;
    int returnMipGap_=0;
    int returnBestDualBound_=0;
    int solnSens_=0;
    int fixModel_=0;
  };
  Options mipStoredOptions_;

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
  bool need_basis_out() const { return 2 & basis(); }

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


public:
  virtual void InitStandardOptions() {
    BaseBackend::InitStandardOptions();
    InitMIPOptions();
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
    {     "1", "Yes (default.)", 1}
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
    {     "3", "Both (1 + 2 = default.)", 3}
  };

  const mp::OptionValueInfo values_warmstart_[3] = {
    {     "0", "No", 0 },
    {     "1", "Yes (for LP: if there is no incoming alg:basis) (default)", 1},
    {     "2", "Yes (for LP: ignoring the incoming alg:basis, if any.)", 2}
  };

  const mp::OptionValueInfo values_rays_[4] = {
    {     "0", "Neither", 0 },
    {     "1", "Just .unbdd", 1},
    {     "2", "Just .dunbdd", 2},
    {     "3", "Both (default).", 3}
  };

  ////////////////////////////////////////////////////////////////
  virtual void InitMIPOptions() {
      if (IMPL_HAS_STD_FEATURE( RETURN_MIP_GAP ))
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
                      "Whether to use or return a basis:\n "
                      "\n.. value-table::\n",
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
                      "Whether to return suffix .unbdd if the objective is unbounded "
                      "or suffix .dunbdd if the constraints are infeasible:\n"
                      "\n.. value-table::\n",
                      GetMIPOptions().rays_, values_rays_);

    if (IMPL_HAS_STD_FEATURE( IIS ))
      AddStoredOption("alg:iisfind iisfind iis",
                      "Whether to find and export the IIS. "
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
                      "ranges of values for which the optimal basis remains optimal:\n"
                        "\n"
                        "|  0 - No (default)\n"
                        "|  1 - Yes:  suffixes return on variables are\n"
                        "|    .sensobjlo = smallest objective coefficient\n"
                        "|    .sensobjhi = greatest objective coefficient\n"
                        "|    .senslblo = smallest variable lower bound\n"
                        "|    .senslbhi = greatest variable lower bound\n"
                        "|    .sensublo = smallest variable upper bound\n"
                        "|    .sensubhi = greatest variable upper bound;"
                              " suffixes for constraints are\n"
                        "|    .sensrhslo = smallest right-hand side value\n"
                        "|    .sensrhshi = greatest right-hand side value.",
                    GetMIPOptions().solnSens_);

    if (IMPL_HAS_STD_FEATURE( FIX_MODEL ))
      AddStoredOption("mip:basis fixmodel mip:fix",
                      "Whether to compute duals / basis / sensitivity for MIP models:\n"
                      "\n.. value-table::\n",
                    GetMIPOptions().fixModel_, values_01_noyes_0default_);
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

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
  DEFINE_STD_FEATURE( BASIS, false )
  std::vector<int> VarStatii() { return {}; }
  void VarStatii(ArrayRef<int> )
  { UNSUPPORTED("MIPBackend::VarStatii"); }
  std::vector<int> ConStatii() { return {}; }
  void ConStatii(ArrayRef<int> )
  { UNSUPPORTED("MIPBackend::ConStatii"); }
  /**
  * Obtain unbounded/inf rays
  **/
  DEFINE_STD_FEATURE( RAYS, false )
  std::vector<double> Ray() { return {}; }
  std::vector<double> DRay() { return {}; }
  /**
  * Compute the IIS and obtain relevant values
  **/
  DEFINE_STD_FEATURE( IIS, false )
  void ComputeIIS() {}
  /// Elements correspond to IISStatus
  std::vector<int> ConsIIS() { return {}; }
  std::vector<int> VarsIIS() { return {}; }
  /**
  * Get MIP Gap
  **/
  DEFINE_STD_FEATURE( ReturnMIPGap, false )
  double MIPGap() const { return MP_DISPATCH( Infinity() ); }
  double MIPGapAbs() const {
    return std::fabs(
          MP_CONST_DISPATCH(ObjectiveValue()) -
          MP_CONST_DISPATCH(BestDualBound()) );
  }
  /**
  * Get MIP dual bound
  **/
  DEFINE_STD_FEATURE( ReturnBestDualBound, false )
  double BestDualBound() const
  { UNSUPPORTED("BestDualBound()"); return 0; }


  ////////////////////////////////////////////////////////////////////////////
  /////////////////////// MIP specific derived calculations //////////////////
  ////////////////////////////////////////////////////////////////////////////

  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         INPUT         //////////////////////////
  using BaseBackend::ReadSuffix;
  using BaseBackend::ReportSuffix;

  void ReadStandardSuffixes() {
    BasicBackend<Impl>::ReadStandardSuffixes();
    ReadStandardMIPSuffixes();
  }

  void ReadStandardMIPSuffixes() {
    if (1 & mipStoredOptions_.basis_)
      MP_DISPATCH( ReadBasis() );
  }

  void ReadBasis() {
    MP_DISPATCH( VarStatii(ReadSuffix(suf_varstatus)) );
    MP_DISPATCH( ConStatii(ReadSuffix(suf_constatus)) );
  }

  //////////////////////// STANDARD MIP SUFFIXES //////////////////////////
  ////////////////////////         OUNPUT        //////////////////////////
  void ReportStandardSuffixes() {
    BaseBackend::ReportStandardSuffixes();
    ReportStandardMIPSuffixes();
  }

  void ReportStandardMIPSuffixes() {
    if (2 & mipStoredOptions_.basis_)
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
    int rays_=3;
    int exportIIS_=0;
    int returnMipGap_=0;
    int returnBestDualBound_=0;
  };
  Options mipStoredOptions_;

public:

  void InitStandardOptions() {
    BaseBackend::InitStandardOptions();
    MP_DISPATCH( InitMIPOptions() );
  }

  using BaseBackend::AddStoredOption;

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

  const mp::OptionValueInfo values_rays_[4] = {
      {     "0", "neither", 0 },
      {     "1", "just .unbdd", 1},
      {     "2", "just .dunbdd", 2},
      {     "3", "both (default)", 3}
  };

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

    if (IMPL_HAS_STD_FEATURE( RAYS ))
      AddStoredOption("mip:rays rays",
                      "Whether to return suffix .unbdd if the objective is unbounded "
                      "or suffix .dunbdd if the constraints are infeasible:\n"
                      "\n.. value-table::\n",
                      mipStoredOptions_.rays_, values_rays_);

    if (IMPL_HAS_STD_FEATURE( IIS ))
      AddStoredOption("mip:iisfind iisfind",
                      "Whether to find and export the IIS. "
                      "Default = 0 (don't export).",
                      mipStoredOptions_.exportIIS_);

    if (IMPL_HAS_STD_FEATURE( ReturnMIPGap ))
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

    if (IMPL_HAS_STD_FEATURE( ReturnBestDualBound ))
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


protected:
  const Options& GetStandardMIPOptions() const { return mipStoredOptions_; }

  bool need_ray_primal() const { return 1 & mipStoredOptions_.rays_; }
  bool need_ray_dual() const { return 2 & mipStoredOptions_.rays_; }



  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////// STANDARD MIP SUFFIXES ///////////////////////////
  //////////////////////////////////////////////////////////////////////////////
private:

  const SuffixDef<int> suf_varstatus = { "sstatus", suf::VAR | suf::OUTPUT };
  const SuffixDef<int> suf_constatus = { "sstatus", suf::CON | suf::OUTPUT };

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

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
  /////////////////// STD FEATURE FLAGS //////////////////////
  //////// Disable optional std features by default //////////
  ////////////////////////////////////////////////////////////
  ALLOW_STD_FEATURE( IIS, false )

  ////////////////////////////////////////////////////////////
  //// Override in the Impl for standard MIP calculations ////
  ////////////////////////////////////////////////////////////
  /**
  * Get/Set AMPL var/con statii
  **/
  std::vector<int> VarStatii() { return {}; }
  void VarStatii(ArrayRef<int> ) { }
  std::vector<int> ConStatii() { return {}; }
  void ConStatii(ArrayRef<int> ) { }
  /**
  * Compute the IIS and relevant values
  **/
  void ComputeIIS();
  /**
  * Get IIS values for constraints / vars
  **/
  std::vector<int> ConsIIS() { return {}; }
  std::vector<int> VarsIIS() { return {}; }
  /**
  * Get MIP Gap
  **/
  double MIPGap() { throw std::runtime_error("Not implemented"); }


  ////////////////////////////////////////////////////////////
  /////////////////// MIP Backend options ////////////////////
  ////////////////////////////////////////////////////////////
  void InitOptions() {
    BaseBackend::InitOptions();
    MP_DISPATCH( InitMIPOptions() );
  }

  using BaseBackend::AddStoredOption;
  void InitMIPOptions() {
    if (IMPL_HAS_STD_FEATURE( IIS ))
      AddStoredOption("iisfind",
                      "Whether to find and export the IIS. "
                      "Default = 0 (don't export).",
                      mipStoredOptions_.exportIIS_);
    AddStoredOption("return_mipgap",
      "Whether to return mipgap suffixes or include mipgap values\n\
		(|objectve - best_bound|) in the solve_message:  sum of\n\
			1 = return relmipgap suffix (relative to |obj|);\n\
			2 = return absmipgap suffix (absolute mipgap);\n\
			4 = suppress mipgap values in solve_message.\n\
		Default = 0.  The suffixes are on the objective and problem.\n\
		Returned suffix values are +Infinity if no integer-feasible\n\
		solution has been found, in which case no mipgap values are\n\
		reported in the solve_message.",
      mipStoredOptions_.returnMipGap_);
  }


  ////////////////////////////////////////////////////////////////////////////
  /////////////////////// MIP specific derived calculations //////////////////
  ////////////////////////////////////////////////////////////////////////////
  void CalculateAndReportDerivedResults() {
    BaseBackend::CalculateAndReportDerivedResults();
    CalculateAndReportDerivedResults_MIP();
  }

  void CalculateAndReportDerivedResults_MIP() {
    CalculateAndReportIIS();
    CalculateAndReportMIPGap();
  }

  using BaseBackend::ReadSuffix;
  using BaseBackend::ReportSuffix;

  void CalculateAndReportIIS() {
    if (MP_DISPATCH( IsProblemInfOrUnb() ) &&
        mipStoredOptions_.exportIIS_) {
      MP_DISPATCH( ComputeIIS() );

      ReportSuffix(sufIISCon, MP_DISPATCH(ConsIIS()));
      ReportSuffix(sufIISVar, MP_DISPATCH(VarsIIS()));
    }
  }

  void CalculateAndReportMIPGap() {
    if (MP_DISPATCH(IsMIP()) &&
        mipStoredOptions_.returnMipGap_) {
      MP_DISPATCH( ComputeMipGap() );

      std::vector<double> dbl(1, MP_DISPATCH( MIPGap() ));
      ReportSuffix(sufRelMipGapObj, dbl);
      ReportSuffix(sufRelMipGapProb, dbl);
    }
  }

  void ComputeMipGap() {}

  //////////////////////// STADDARD MIP SUFFIXES //////////////////////////
  void ReadStandardSuffixes() {
    BasicBackend<Impl>::ReadStandardSuffixes();
    ReadStandardMIPSuffixes();
  }

  void ReadStandardMIPSuffixes() {
    MP_DISPATCH( VarStatii(ReadSuffix(suf_varstatus)) );
    MP_DISPATCH( ConStatii(ReadSuffix(suf_constatus)) );
  }

  void ReportStandardSuffixes() {
    BasicBackend<Impl>::ReportStandardSuffixes();
    ReportStandardMIPSuffixes();
  }

  void ReportStandardMIPSuffixes() {
    ReportSuffix(suf_varstatus,
                              MP_DISPATCH( VarStatii() ));
    ReportSuffix(suf_constatus,
                              MP_DISPATCH( ConStatii() ));
  }

private:
  struct Options {
    int exportIIS_;
    int returnMipGap_;
  };
  Options mipStoredOptions_;


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////// STANDARD SUFFIXES ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  const SuffixDef<int> suf_varstatus = { "sstatus", suf::VAR | suf::OUTPUT };
  const SuffixDef<int> suf_constatus = { "sstatus", suf::CON | suf::OUTPUT };

  const SuffixDef<int> sufIISCon = { "iis", suf::CON | suf::OUTPUT };
  const SuffixDef<int> sufIISVar = { "iis", suf::VAR | suf::OUTPUT };

  const SuffixDef<double> sufRelMipGapObj = { "relmipgap", suf::OBJ | suf::OUTPUT };
  const SuffixDef<double> sufRelMipGapProb = { "relmipgap", suf::PROBLEM | suf::OUTPUT };

};

}  // namespace mp

#endif  // MIPBACKEND_H_

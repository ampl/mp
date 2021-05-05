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
  void CalculateDerivedResults() {
    BasicBackend<Impl>::CalculateDerivedResults();
    if (MP_DISPATCH( IsProblemInfOrUnb() ) &&
        mipStoredOptions_.exportIIS_)
      MP_DISPATCH( ComputeIIS() );
    if (mipStoredOptions_.returnMipGap_)
      MP_DISPATCH( ComputeMipGap() );

  }

  // Override in the Impl for standard MIP calculations
  /**
  * Compute the IIS and relevant values
  **/
  void ComputeIIS();
  /**
  * Get IIS values for constraints 
  **/
  void ConsIIS(std::vector<int>& stt) { stt.clear(); }
  /**
  * Get IIS values for variables
  **/
  void VarsIIS(std::vector<int>& stt) { stt.clear(); }

  void ComputeMipGap() {}
  double MIPGap() { throw std::runtime_error("Not implemented"); }

  const SuffixDef<int> sufIISCon = { "iis", suf::CON | suf::OUTPUT };
  const SuffixDef<int> sufIISVar = { "iis", suf::VAR | suf::OUTPUT };

  const SuffixDef<double> sufRelMipGapObj = { "relmipgap", suf::OBJ | suf::OUTPUT };
  const SuffixDef<double> sufRelMipGapProb = { "relmipgap", suf::PROBLEM | suf::OUTPUT };

  void ReportStandardSuffixes() {
    BasicBackend<Impl>::ReportStandardSuffixes();
    std::vector<int> stt;
    std::vector<double> dbl;    
    if (mipStoredOptions_.exportIIS_)
    {
      MP_DISPATCH(ConsIIS(stt));
      this->DeclareAndReportIntSuffix(sufIISCon, stt);
      MP_DISPATCH(VarsIIS(stt));
      this->DeclareAndReportIntSuffix(sufIISVar, stt);
    }
    if (mipStoredOptions_.returnMipGap_)
    {
      dbl.clear();
      dbl.push_back(MP_DISPATCH( MIPGap() ));
      this->DeclareAndReportDblSuffix(sufRelMipGapObj, dbl);
      this->DeclareAndReportDblSuffix(sufRelMipGapProb, dbl);
    }
  }

private:
  struct Options {
    int exportIIS_;
    int returnMipGap_;
  };
  Options mipStoredOptions_;

};

}  // namespace mp

#endif  // MIPBACKEND_H_

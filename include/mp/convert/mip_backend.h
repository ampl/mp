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
  struct Options {
    int exportIIS_;
  };
public:

  void InitMetaInfoAndOptions() {
    MP_DISPATCH(InitNamesAndVersion());
    InitMIPOptions();
    MP_DISPATCH(InitOptions());
  }

  Options mipstoredOptions_;

  void InitMIPOptions() {
    AddStoredOption("iisfind",
      "Whether to export IIS. "
      "Default = 0 (don't export).",
      mipstoredOptions_.exportIIS_);// , TRUEORFALSE);
  }
  void CalculateDerivedResults() {
    BasicBackend<Impl>::CalculateDerivedResults();
    if (mipstoredOptions_.exportIIS_)
    {
      MP_DISPATCH(ComputeIIS());
    }

  }
  void ComputeIIS() {}
  void ConsIIS(std::vector<int>& stt) { stt.clear(); };
  void VarsIIS(std::vector<int>& stt) { stt.clear(); };

  const SuffixDef<int> suf_coniis = { "iis", suf::CON | suf::OUTPUT };
  const SuffixDef<int> suf_variis = { "iis", suf::VAR | suf::OUTPUT };

  void ReportStandardSuffixes() {
    BasicBackend<Impl>::ReportStandardSuffixes();
    std::vector<int> stt;
    if (mipstoredOptions_.exportIIS_)
    {
      MP_DISPATCH(ConsIIS(stt));
      DeclareAndReportIntSuffix(suf_coniis, stt);
      MP_DISPATCH(VarsIIS(stt));
      DeclareAndReportIntSuffix(suf_variis, stt);
    }
  }

};
}  // namespace mp

#endif  // MIPBACKEND_H_

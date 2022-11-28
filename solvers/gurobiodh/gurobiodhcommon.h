#ifndef GUROBIODHCOMMON_H
#define GUROBIODHCOMMON_H


//#define ODH
/// Common stuff for GurobiBackend and GurobiModelAPI
extern "C" {
  #include "heuristic.h"
}
#include <stdexcept>
#include <cassert>
#include "mp/format.h"

namespace mp {

  /// Information shared by both
  /// `GurobiBackend` and `GurobiModelAPI`
  struct GurobiODHCommonInfo {

    HEURENVptr HEURptr() {
      assert(HEURptr_);
      return HEURptr_;
    }
    void openODH(void* grbenv)
    {
      int status;
      HEURptr_ = HEURopen(&status, grbenv);
      if (status)
        throw std::runtime_error(fmt::format("Failed to open Heuristic environment:\n{}",
          HEURgetlastmsg(HEURptr())));
    }
    const char* getStatusMsg(int status)
    {
      char returnmsg[16][32] = {
      "HEUR_STAT_NOTCALLED",
      "HEUR_STAT_INF",
      "HEUR_STAT_FAIL_INF",
      "HEUR_STAT_ABORT_INF",
      "HEUR_STAT_TLIM_INF",
      "HEUR_STAT_DVSR_INF",
      "HEUR_STAT_FEAS",
      "HEUR_STAT_FAIL_FEAS",
      "HEUR_STAT_ABORT_FEAS",
      "HEUR_STAT_TLIM_FEAS",
      "HEUR_STAT_DVSR_FEAS",
      "HEUR_STAT_OPT",
      "HEUR_STAT_FAIL_OPT",
      "HEUR_STAT_ABORT_OPT",
      "HEUR_STAT_TLIM_OPT",
      "HEUR_STAT_DVSR_OPT"
      };
      return returnmsg[status];
    }
    HEURENVptr HEURptr_ = NULL;
  };


  /// Convenience macro
#define ODH_CCALL( call ) do { if (int e = (call) != 0) \
  throw std::runtime_error( \
    fmt::format("ODH Call failed: '{}':\n{}", #call, HEURgetlastmsg( HEURptr() ) )); } while (0)
}

#endif // GUROBIODHCOMMON_H
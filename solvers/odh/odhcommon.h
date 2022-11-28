#ifndef ODHCOMMON_H
#define ODHCOMMON_H


extern "C" {
  #include "heuristic.h"
}
#include <stdexcept>
#include <cassert>
#include "mp/format.h"

namespace mp {

  struct ODHCommonInfo {

    HEURENVptr HEURptr() {
      assert(HEURptr_);
      return HEURptr_;
    }
    void OpenODH(void* env)
    {
      int status;
      HEURptr_ = HEURopen(&status, env);
      if (status)
        throw std::runtime_error(
          fmt::format("Failed to open Heuristic environment:\n{}\n",
          HEURgetlastmsg(HEURptr())));
    }
    void CloseODH()
    {
      if (HEURptr_)
      {
        int status = HEURclose(&HEURptr_);
        throw std::runtime_error(
          fmt::format("Failed to close Heuristic environment:\n{}\n",
          HEURgetlastmsg(HEURptr())));
      }
    }
    const char* GetODHStatusMsg(int status)
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

#endif // ODHCOMMON_H
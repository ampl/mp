#ifndef ODHCOMMON_H
#define ODHCOMMON_H


extern "C" {
  #include "heuristic.h"
}
#include <stdexcept>
#include <cassert>
#include "mp/format.h"
#include "mp/common.h"
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


    std::pair<int, std::string> ConvertODHStatus() {
      int status = HEURgetstat(HEURptr());
      namespace sol = mp::sol;
      switch (status) {
      case 0:
        return { sol::NOT_SET, "optimizer not called" };
      case 1:
        return { sol::INFEASIBLE, "infeasible problem" };
      case 2:
        return { sol::FAILURE, "failure without a feasible solution" };
      case 3:
        return { sol::INTERRUPTED, "interrupted without a feasible solution" };
      case 4:
        return { sol::LIMIT, "time limit reached without a feasible solution" };
      case 5:
        return { sol::LIMIT, "maximum divisor reached without a feasible solution" };
      case 6:
        return { sol::SOLVED, "feasible solution found" };
      case 7:
        return { sol::FAILURE, "failure with a feasible solution" };
      case 8:
        return { sol::INTERRUPTED, "interrupted with a feasible solution" };
      case 9:
        return { sol::LIMIT, "time limit reached with a feasible solution" };
      case 10:
        return { sol::LIMIT, "maximum divisor reached with a feasible solution" };
      case 11:
        return { sol::SOLVED, "optimal solution" };
      case 12:
        return { sol::FAILURE, "failure with an optimal solution" };
      case 13:
        return { sol::INTERRUPTED, "interrupted with an optimal solution" };
      case 14:
        return { sol::LIMIT, "time limit reached with an optimal solution" };
      case 15:
        return { sol::LIMIT, "maximum divisor reached with an optimal solution" };
      default:
        return { sol::UNKNOWN, "unexpected solution status" };
      }
    }
    
    HEURENVptr HEURptr_ = NULL;
  };


  /// Convenience macro
#define ODH_CCALL( call ) do { if (int e = (call) != 0) \
  throw std::runtime_error( \
    fmt::format("ODH Call failed: '{}':\n{}", #call, HEURgetlastmsg( HEURptr() ) )); } while (0)
}

#endif // ODHCOMMON_H
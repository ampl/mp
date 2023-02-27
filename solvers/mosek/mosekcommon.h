#ifndef MOSEKCOMMON_H
#define MOSEKCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"
#include "mp/format.h"
#include "mp/error.h"

extern "C" {
    #include "mosek.h"
}

namespace mp {

/// Information shared by both `MosekBackend` and `MosekModelAPI`
struct MosekCommonInfo {
  MSKtask_t lp() const { return lp_; }
  void set_lp(MSKtask_t lp) { lp_ = lp; }
  MSKtask_t& lp_ref() { return lp_; }

  MSKenv_t env() const { return env_; }
  void set_env(MSKenv_t env) { env_ = env; }
  MSKenv_t& env_ref() { return env_; }
private:
  MSKtask_t      lp_ = NULL;
  MSKenv_t       env_ = NULL;
};

/// Common API for Mosek classes
class MosekCommon :
    public Backend2ModelAPIConnector<MosekCommonInfo> {
public:
  /// These methods access Mosek options. Used by AddSolverOption()
  void GetSolverOption(MSKiparame key, int& value) const;
  void SetSolverOption(MSKiparame key, int value);
  void GetSolverOption(MSKdparame key, double& value) const;
  void SetSolverOption(MSKdparame key, double value);
  void GetSolverOption(MSKsparame key, std::string& value) const;
  void SetSolverOption(MSKsparame key, const std::string& value);

  // TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  int getIntAttr(MSKiinfitem_enum name) const;
  long long getLongAttr(MSKliinfitem_enum name) const;
  double getDblAttr(MSKdinfitem_enum name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

protected:
};

#define MOSEK_CCALL( call ) do { \
	MSKrescodee e=call; \
	if (e != MSK_RES_OK) handleError(e, #call); \
} while (0)

/// Internal function
void handleError(MSKrescodee e, const char* fname);

} // namespace mp

#endif // MOSEKCOMMON_H

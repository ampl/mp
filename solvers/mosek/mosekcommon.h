#ifndef MOSEKCOMMON_H
#define MOSEKCOMMON_H

#include <map>
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
  /// These methods access Mosek options. Used by AddSolverOption().
  /// A value rejected by MOSEK raises with the AMPL option name.
  void GetSolverOption(MSKiparame key, int& value) const;
  void SetSolverOption(MSKiparame key, int value);
  void GetSolverOption(MSKdparame key, double& value) const;
  void SetSolverOption(MSKdparame key, double value);
  void GetSolverOption(MSKsparame key, std::string& value) const;
  void SetSolverOption(MSKsparame key, const std::string& value);

  /// Remember the AMPL option name of a parameter (for messages)
  void RecordOptionName(MSKparametertypee type, int key,
                        const std::string& ampl_name)
  { optNames_[{type, key}] = ampl_name; }
  /// The MOSEK parameter name, with the AMPL option name if recorded
  std::string ParamName(MSKparametertypee type, int key) const;

  // TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  /// Raise on a failed parameter set
  void CheckParamSet(MSKrescodee code, MSKparametertypee type, int key,
                     const std::string& value) const;

  int getIntAttr(MSKiinfitem_enum name) const;
  long long getLongAttr(MSKliinfitem_enum name) const;
  double getDblAttr(MSKdinfitem_enum name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

private:
  std::map<std::pair<int, int>, std::string> optNames_;
};

#define MOSEK_CCALL( call ) do { \
	MSKrescodee e=call; \
	if (e != MSK_RES_OK) handleError(e, #call); \
} while (0)

/// Internal function
void handleError(MSKrescodee e, const char* fname);

} // namespace mp

#endif // MOSEKCOMMON_H

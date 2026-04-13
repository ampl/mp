#ifndef COPTCOMMON_H
#define COPTCOMMON_H

#include <string>
#include <map>

extern "C" {
  #include "copt.h"
}

#include "mp/error.h"
#include "mp/backend-to-model-api.h"
#include "mp/format.h"

namespace mp {

/// Information inherited by both
/// `CoptBackend` and `CoptModelAPI`
struct CoptCommonInfo {
  copt_env* env() const { return env_; }
  copt_env*& env_ref() { return env_; }
  copt_prob* lp() const { return lp_; }
  copt_prob*& lp_ref() { return lp_; }

  void set_env(copt_env* e) { env_ = e; }
  void set_lp(copt_prob* lp) { lp_ = lp; }



  int num_lin_obj_ {0};          // @note this is not accessible in Backend
                                 // unless via get_other()
  int obj_sense0_ampl_ {-1};     // 0: min, 1: max
  void set_current_objective_options(int options) {
      _current_objective_options = options;
  }
  int current_objective_options() const {
	  return _current_objective_options;
  }
  
private:
  copt_env*      env_ = NULL;
  copt_prob*      lp_ = NULL;
  int _current_objective_options = -1;
};


/// Common API for Copt classes
class CoptCommon :
    public Backend2ModelAPIConnector<CoptCommonInfo> {
public:
  /// These methods access Copt options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  static constexpr double Infinity() { return COPT_INFINITY; }
  static constexpr double MinusInfinity() { return -COPT_INFINITY; }

protected:
  int getIntAttr(const char* name) const;
  double getDblAttr(const char* name) const;
  void setIntAttr(const char* name, int value);
  void setDblAttr(const char* name, double value);

  std::vector<double> getVarInfo(const char* name);
  std::vector<double> getConInfo(const char* name);

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;


private:

  int (*createEnv) (copt_env**) = nullptr;
};


/// Convenience macro
#define COPT_CCALL( call ) do { if (int e = (call) != COPT_RETCODE_OK) { \
  char buf[512] = ""; \
  COPT_GetRetcodeMsg(e, buf, 512); \
  MP_RAISE( \
    fmt::format("  Call failed: '{}' with code {}:\n" \
                "{}", #call, e, buf )); } } while (0)

} // namespace mp

#endif // COPTCOMMON_H

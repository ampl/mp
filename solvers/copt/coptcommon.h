#ifndef COPTCOMMON_H
#define COPTCOMMON_H

#include <string>

extern "C" {
  #include "copt.h"
}

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


private:
  copt_env*      env_ = NULL;
  copt_prob*      lp_ = NULL;
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
  void OpenSolver();
  void CloseSolver();

  int getIntAttr(const char* name) const;
  double getDblAttr(const char* name) const;

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
#define COPT_CCALL( call ) do { if (int e = (call) != COPT_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // COPTCOMMON_H

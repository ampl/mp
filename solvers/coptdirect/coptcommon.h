#ifndef COPTCOMMON_H
#define COPTCOMMON_H

#include <string>

extern "C" {
  #include "copt.h"
}

#include "mp/format.h"

namespace mp {

// Common ancestor for Copt classes
class CoptCommon
{
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

  /// Connection between Backend and ModelAPI
  CoptCommon *other_copt() { return other_; }
  /// Set connection
  void set_other_copt(CoptCommon* o) { other_ = o; }

  void setEnvCreate(int (*createEnvFunction) (copt_env**)) {
    createEnv = createEnvFunction;
  }

protected:
  void OpenSolver();
  void CloseSolver();

  // TODO check if here it makes sense

  int getIntAttr(const char* name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

  copt_env* env() const { return env_; }
  void set_env(copt_env* e) { env_ = e; }
  copt_prob* lp() const { return lp_; }
  void set_lp(copt_prob* lp) { lp_ = lp; }

  void copy_handlers_from_other_copt();
  void copy_handlers_to_other_copt();


private:
  copt_env*      env_ = NULL;
  copt_prob*      lp_ = NULL;
  CoptCommon *other_ = nullptr;

  int (*createEnv) (copt_env**) = nullptr;
  

};


/// Convenience macro
#define COPT_CCALL( call ) do { if (int e = (call) != COPT_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // COPTCOMMON_H

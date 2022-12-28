#ifndef CPLEXCOMMON_H
#define CPLEXCOMMON_H

#include <string>

extern "C" {
  #include <ilcplex/cplex.h>
}

#include "mp/backend-to-model-api.h"
#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `CplexBackend` and `CplexModelAPI`
struct CplexCommonInfo {
  CPXENVptr env() const { return env_; }
  CPXENVptr& env_ref() { return env_; }
  CPXLPptr lp() const { return lp_; }
  CPXLPptr& lp_ref() { return lp_; }

  void set_env(CPXENVptr e) { env_ = e; }
  void set_lp(CPXLPptr lp) { lp_ = lp; }

private:
  CPXENVptr     env_ = NULL;
  CPXLPptr      lp_ = NULL;
};


/// Common API for Cplex classes
class CplexCommon :
    public Backend2ModelAPIConnector<CplexCommonInfo> {
public:
  /// These methods access CPLEX options. Used by AddSolverOption()
  void GetSolverOption(int key, int& value) const;
  void SetSolverOption(int key, int value);
  void GetSolverOption(int key, double& value) const;
  void SetSolverOption(int key, double value);
  void GetSolverOption(int key, std::string& value) const;
  void SetSolverOption(int key, const std::string& value);

  static constexpr double Infinity() { return CPX_INFBOUND; }
  static constexpr double MinusInfinity() { return -CPX_INFBOUND; }

protected:
  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int ModelSense() const;

};


/// Convenience macro
#define CPLEX_CALL( call ) do { if (int e=call) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // CPLEXCOMMON_H

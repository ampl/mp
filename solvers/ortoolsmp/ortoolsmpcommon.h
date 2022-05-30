#ifndef ORTOOLSCOMMON_H
#define ORTOOLSCOMMON_H

#include <string>
#include <limits>

#include "mp/backend-to-model-api.h"

#include "ortools/linear_solver/linear_solver.h"


/// The below would go into actual ...common.h:

#include "mp/format.h"

namespace mp {
/// Information shared by both
/// `OrtoolsBackend` and `OrtoolsModelAPI`
struct OrtoolsCommonInfo {
  operations_research::MPSolver* lp() const { return lp_; }
  void set_lp(operations_research::MPSolver* lp) { lp_ = lp; }
private:
  operations_research::MPSolver*      lp_ = NULL;
};


/// Common API for Ortools classes
class OrtoolsCommon :
    public Backend2ModelAPIConnector<OrtoolsCommonInfo> {
public:
  /// These methods access Ortools options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  // TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return std::numeric_limits<double>::infinity();  }
  static constexpr double MinusInfinity() { return -std::numeric_limits<double>::infinity(); }

protected:
  int getIntAttr(int name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

protected:
  // TODO if desirable, provide function to create the solver's environment
  //int (*createEnv) (solver_env**) = nullptr;
  
};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define ORTOOLS_RETCODE_OK 0
#define ORTOOLS_CCALL( call ) do { if (int e = (call) != ORTOOLS_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // ORTOOLSCOMMON_H

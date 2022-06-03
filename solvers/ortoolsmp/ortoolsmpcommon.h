#ifndef ORTOOLSCOMMON_H
#define ORTOOLSCOMMON_H

#include <string>
#include <limits>

#include "mp/backend-to-model-api.h"

#include "ortools/linear_solver/linear_solver.h"
namespace orr = operations_research;

/// The below would go into actual ...common.h:

#include "mp/format.h"

namespace mp {
/// Information shared by both
/// `OrtoolsBackend` and `OrtoolsModelAPI`
struct OrtoolsCommonInfo {
  orr::MPSolver* lp() const { return lp_; }
  void set_lp(orr::MPSolver* lp) { lp_ = lp; }

 private:
  orr::MPSolver* lp_ = NULL;
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

  double Infinity() { return lp()->infinity(); }
  double MinusInfinity() { return -Infinity(); }

protected:
  std::map<std::string, int> paramNames_ = {
      {"RELATIVE_MIP_GAP",
       orr::MPSolverParameters::RELATIVE_MIP_GAP},
      {"PRIMAL_TOLERANCE",
       orr::MPSolverParameters::PRIMAL_TOLERANCE},
      {"DUAL_TOLERANCE",
       orr::MPSolverParameters::DUAL_TOLERANCE},

      {"PRESOLVE", orr::MPSolverParameters::PRESOLVE},
      {"LP_ALGORITHM", orr::MPSolverParameters::LP_ALGORITHM},
      {"INCREMENTALITY",
       orr::MPSolverParameters::INCREMENTALITY},
      {"SCALING", orr::MPSolverParameters::SCALING}};

  orr::MPSolverParameters params_;

  int getIntAttr(int name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;
};


// This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define ORTOOLS_RETCODE_OK 0
#define ORTOOLS_CCALL( call ) do { if (int e = (call) != ORTOOLS_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // ORTOOLSCOMMON_H

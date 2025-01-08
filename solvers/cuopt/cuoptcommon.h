#ifndef CUOPTCOMMON_H
#define CUOPTCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

//extern "C" {
// TODO Typically import here the solver's C API headers
// Here we use a cplus plus stub instead
   #include "cuopt-solvermodel.h"
//}

#include "mp/format.h"

#include "json.hpp"

// for convenience
using json = nlohmann::json;

#include <iostream>

namespace mp {

/// Information shared by both
/// `CuoptBackend` and `CuoptModelAPI`
struct CuoptCommonInfo {

  json* get_json_prob() { return json_prob_; }
  void set_json_prob(json* json_prob) { json_prob_ = json_prob; }

  bool isMIP() const { return *ismip_; }
  void SetIsMIP(bool ismip) { *ismip_ = ismip; }

private:
  json* json_prob_;
  bool* ismip_ = new bool(false);
};


/// Common API for Cuopt classes
class CuoptCommon :
    public Backend2ModelAPIConnector<CuoptCommonInfo> {
public:
  /// These methods access Cuopt options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  /// TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  int getIntAttr(Solver::ATTRIBS name, Solver::ConsType subtype=Solver::ConsType::CONS_LIN) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;


protected:
  // TODO if desirable, provide function to create the solver's environment
  // with own license
  // int (*createEnv) (solver_env**) = nullptr;
  
};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define CUOPT_RETCODE_OK 0
#define CUOPT_CCALL( call ) do { if (int e = (call) != CUOPT_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // CUOPTCOMMON_H

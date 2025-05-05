#ifndef CUOPTLPCOMMON_H
#define CUOPTLPCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

//extern "C" {
// TODO Typically import here the solver's C API headers

#include <cuopt/linear_programming/cuopt_c.h>

// Here we use a cplus plus stub instead
   #include "cuoptlp-solvermodel.h"
//}

#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `CuoptlpBackend` and `CuoptlpModelAPI`
struct CuoptlpCommonInfo {

  // TODO provide accessors to the solver's in-memory model/environment
  //cuoptlp_env* env() const { return env_; }
  Solver::SolverModel* lp() const { return lp_; }

  // TODO provide accessors to the solver's in-memory model/environment
  //void set_env(cuoptlp_env* e) { env_ = e; }
  void set_lp(Solver::SolverModel* lp) { lp_ = lp; }


private:
  // TODO provide accessors to the solver's in-memory model/environment
  //cuoptlp_env*      env_ = NULL;
  Solver::SolverModel*      lp_ = NULL;

};


/// Common API for Cuoptlp classes
class CuoptlpCommon :
    public Backend2ModelAPIConnector<CuoptlpCommonInfo> {
public:
  /// These methods access Cuoptlp options. Used by AddSolverOption()
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
#define CUOPTLP_RETCODE_OK 0
#define CUOPTLP_CCALL( call ) do { if (int e = (call) != CUOPTLP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // CUOPTLPCOMMON_H

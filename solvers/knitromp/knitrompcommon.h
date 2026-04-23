#ifndef KNITROMPCOMMON_H
#define KNITROMPCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

//extern "C" {
// TODO Typically import here the solver's C API headers
// Here we use a cplus plus stub instead
   #include "knitro.h"
//}

#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `KnitrompBackend` and `KnitrompModelAPI`
struct KnitrompCommonInfo {

  // TODO provide accessors to the solver's in-memory model/environment
    KN_context_ptr lp() const { return lp_; };
    KN_context_ptr& lp_ref() { return lp_; }
  // TODO provide accessors to the solver's in-memory model/environment
  //void set_env(knitromp_env* e) { env_ = e; }
  void set_lp(KN_context_ptr lp) { lp_ = lp; }


private:
    KN_context_ptr      lp_ = NULL;

};


/// Common API for Knitromp classes
class KnitrompCommon :
    public Backend2ModelAPIConnector<KnitrompCommonInfo> {
public:
  /// These methods access Knitromp options. Used by AddSolverOption()
  void GetSolverOption(int key, int& value) const;
  void SetSolverOption(int key, int value);
  void GetSolverOption(int key, double& value) const;
  void SetSolverOption(int key, double value);
  void GetSolverOption(int key, std::string& value) const;
  void SetSolverOption(int key, const std::string& value);
  static constexpr double Infinity() { return KN_INFINITY;  }
  static constexpr double MinusInfinity() { return -KN_INFINITY; }

protected:
  //int getIntAttr(Solver::ATTRIBS name, Solver::ConsType subtype=Solver::ConsType::CONS_LIN) const;
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
#define KNITROMP_RETCODE_OK 0
#define KNITROMP_CCALL( call ) do { if (int e = (call) != KNITROMP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // KNITROMPCOMMON_H

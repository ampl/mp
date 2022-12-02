#ifndef CBCMPCOMMON_H
#define CBCMPCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"
#define CBC_EXTERN_C
  #include "CbcModel.hpp"
  
  #include "CbcSolver.hpp"
#include "Cbc_C_Interface.h"

#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `CbcmpBackend` and `CbcmpModelAPI`
struct CbcmpCommonInfo {
  Cbc_Model* lp() const { return lp_; }
  void set_lp(Cbc_Model* lp) { lp_ = lp; }
private:
  Cbc_Model*      lp_ = NULL;

};


/// Common API for Cbcmp classes
class CbcmpCommon :
    public Backend2ModelAPIConnector<CbcmpCommonInfo> {

  template <typename T> struct Parameter {
    T value_;
    bool isNonDefault_;
    std::function<void(Cbc_Model*, T)> set_;
    std::function<T(Cbc_Model*)> get_;

    public:
    Parameter(T def) : value(def), isNonDefault_(false)
    {}
    
  };
  template <typename T> class OptionsWrapper {
    std::map < std::string, std::function<void(Cbc_Model*, T)>> options;

  public:
    std::function<void(Cbc_Model*, T)> getOption(const char* name) {
      auto o = options.find(name);
      if (o != options.end())
        return o->second;
    }
    void addOption(const char * name, std::function<void(Cbc_Model*, T)> function) {
      options[name] = function;
    }
  };
  
public:
  OptionsWrapper<double> _options;

  void GetCBCParamsList() const;
  /// These methods access Cbcmp options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, std::string value);

  /// TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  int getIntAttr(int name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumSOSCons() const;

protected:
  // TODO if desirable, provide function to create the solver's environment
  // with own license
  // int (*createEnv) (solver_env**) = nullptr;
  
};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define CBCMP_RETCODE_OK 0
#define CBCMP_CCALL( call ) do { if (int e = (call) != CBCMP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // CBCMPCOMMON_H

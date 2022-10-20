#ifndef HIGHSCOMMON_H
#define HIGHSCOMMON_H

#include <string>

extern "C" {
  #include "interfaces/highs_c_api.h"
}

#include "mp/backend-to-model-api.h"
#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `HighsBackend` and `HighsModelAPI`
struct HighsCommonInfo {
  void* lp() const { return lp_; }
  void set_lp(void* lp) { lp_ = lp; }
private:
  void*      lp_ = NULL;
};


/// Common API for Highs classes
class HighsCommon :
    public Backend2ModelAPIConnector<HighsCommonInfo> {
public:
  /// These methods access Highs options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  double myinf = 0;
  double Infinity() {
     if (!myinf) myinf = Highs_getInfinity(lp());
     return myinf;
  }
  double MinusInfinity() { return -Infinity(); }

protected:
  void OpenSolver();
  void CloseSolver();

  int64_t getInt64Attr(const char* name)  const;
  int getIntAttr(const char* name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;
};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define HIGHS_CCALL( call ) do { \
  int e = (call); \
  if (e != kHighsStatusOk && e != kHighsStatusWarning) \
    throw std::runtime_error( \
      fmt::format("  Call failed: '{}' with code {}", #call, e )); \
  if (e == kHighsStatusWarning) \
    std::printf( fmt::format("  WARNING for call: '{}', warning code {}\n", #call, e ).c_str()); \
} while (0)

} // namespace mp

#endif // HIGHSCOMMON_H

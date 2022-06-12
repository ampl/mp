#ifndef MOSEKCOMMON_H
#define MOSEKCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"
#include "mp/format.h"
#include "mp/error.h"

extern "C" {
    #include "mosek.h"
}

namespace mp {

/// Information shared by both `MosekBackend` and `MosekModelAPI`
struct MosekCommonInfo {
  MSKtask_t lp() const { return lp_; }
  void set_lp(MSKtask_t lp) { lp_ = lp; }
  MSKtask_t& lp_ref() { return lp_; }
private:
  MSKtask_t      lp_ = NULL;
};

/// Common API for Mosek classes
class MosekCommon :
    public Backend2ModelAPIConnector<MosekCommonInfo> {
public:
  /// These methods access Mosek options. Used by AddSolverOption()
  void GetSolverOption(MSKiparame key, int& value) const;
  void SetSolverOption(MSKiparame key, int value);
  void GetSolverOption(MSKdparame key, double& value) const;
  void SetSolverOption(MSKdparame key, double value);
  void GetSolverOption(MSKsparame key, std::string& value) const;
  void SetSolverOption(MSKsparame key, const std::string& value);

  // TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  int getIntAttr(MSKiinfitem_enum name) const;
  long long getLongAttr(MSKliinfitem_enum name) const;
  double getDblAttr(MSKdinfitem_enum name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

protected:
};


//MSKrescodee(MSKAPI MSK_getresponseclass) (
  //MSKrescodee r,
 // MSKrescodetypee * rc);



/// Convenience macro
#define MOSEK_RETCODE_OK MSK_RES_ERR_LICENSE
#define MOSEK_CCALL( call ) do { MSKrescodee e= call;\
if (e != MOSEK_RETCODE_OK) {\
MSKrescodee e2; MSKrescodetypee et;\
e2 = MSK_getresponseclass(e, &et);\
if(e2 != MOSEK_RETCODE_OK) MP_RAISE(\
    fmt::format("  Call failed: '{}' with code {}, unknown type", #call, e ));\
if ((int)e2 > (int)MSK_RESPONSE_TRM) MP_RAISE(\
    fmt::format("  Call failed: '{}' with code {}, type {}", #call, e, e2));\
fmt::print("  Warning {} while calling '{}'", #call, e);\
  }} while (0)

} // namespace mp

#endif // MOSEKCOMMON_H

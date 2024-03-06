#ifndef XPRESSMPCOMMON_H
#define XPRESSMPCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

extern "C" {
    #include "xprs.h"
}

#include "mp/format.h"

namespace mp {

struct XpressmpCommonInfo {
  XPRSprob lp() const { return lp_; }
  XPRSprob* lp_ref() { return &lp_; }
  int numIntVars()  const { return numIntVars_; }
  void numIntVars(int num) { numIntVars_= num; }
  
  int numQuadCons() const { return numQuadCons_; }
  void numQuadCons(int num) { numQuadCons_ = num; }
private:
  XPRSprob      lp_ = NULL;
  int numIntVars_ = 0;
  int numQuadCons_ = 0;
};


/// Common API for Xpressmp classes
class XpressmpCommon :
    public Backend2ModelAPIConnector<XpressmpCommonInfo> {
public:
  /// These methods access Xpressmp options. Used by AddSolverOption()
  void GetSolverOption(int key, int& value) const;
  void SetSolverOption(int key, int value);
  void GetSolverOption(int key, double& value) const;
  void SetSolverOption(int key, double value);
  void GetSolverOption(int key, std::string& value) const;
  void SetSolverOption(int key, const std::string& value);

  static constexpr double Infinity() { return XPRS_PLUSINFINITY;  }
  static constexpr double MinusInfinity() { return XPRS_MINUSINFINITY; }

protected:
  int getIntAttr(int attr) const;
  double getDblAttr(int attr) const;
  
  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;
  int NumPWLs() const;
  int NumGenCons() const;

  std::string getErr()  const{
    char errmsg[512];
    if (!XPRSgetlasterror(lp(), errmsg))
      return std::string(errmsg);
    else
      return std::string();
  }
};



#define XPRESSMP_RETCODE_OK 0
#define XPRESSMP_CCALL( call ) do { if (int e = (call) != XPRESSMP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}, message:\n{}\n", #call, e, getErr())); } while (0)

} // namespace mp

#endif // XPRESSMPCOMMON_H

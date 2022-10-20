#include "mp/format.h"
#include "xpressmpcommon.h"

namespace mp {

int XpressmpCommon::getIntAttr(int attr)  const {
  int value = 0;
  XPRESSMP_CCALL(XPRSgetintattrib(lp(), attr, &value));
  return value;
}
double XpressmpCommon::getDblAttr(int attr) const  {
  double value = 0;
  XPRESSMP_CCALL(XPRSgetdblattrib(lp(), attr, &value));
  return value;
}

int XpressmpCommon::NumLinCons() const {
  return getIntAttr(XPRS_ROWS);
}

int XpressmpCommon::NumVars() const {
  return getIntAttr(XPRS_COLS);
}

int XpressmpCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return XPRESSMPgetnumobjs (env_, lp_);
  return 0;
}

int XpressmpCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(XPRESSMP_INTATTR_QCONSTRS);
  return 0;
}

int XpressmpCommon::NumSOSCons() const {
  return getIntAttr(XPRS_SETS);
}

int XpressmpCommon::NumIndicatorCons() const {
  return getIntAttr(XPRS_INDICATORS);
}

void XpressmpCommon::GetSolverOption(int key, int &value) const {
  XPRESSMP_CCALL(XPRSgetintcontrol(lp(), key, &value) );
}

void XpressmpCommon::SetSolverOption(int  key, int value) {
  XPRESSMP_CCALL(XPRSsetintcontrol(lp(), key, value));
}

void XpressmpCommon::GetSolverOption(int  key, double &value) const {
  XPRESSMP_CCALL(XPRSgetdblcontrol(lp(), key, &value));
}

void XpressmpCommon::SetSolverOption(int  key, double value) {
    XPRESSMP_CCALL(XPRSsetdblcontrol(lp(), key, value));
}

void XpressmpCommon::GetSolverOption(int  key, std::string &value) const {
  throw std::runtime_error("Not implemented"); 
}

void XpressmpCommon::SetSolverOption(int  key, const std::string& value) {
  XPRESSMP_CCALL(XPRSsetstrcontrol(lp(), key, value.c_str()));
}




} // namespace mp

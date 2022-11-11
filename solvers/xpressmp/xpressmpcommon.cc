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
  return 1;
}

int XpressmpCommon::NumQPCons() const {
  return getIntAttr(XPRS_QCONSTRAINTS);
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
  int l;
  XPRESSMP_CCALL(XPRSgetstringcontrol(lp(), key, nullptr, 0, &l));
  std::vector<char> s(l);
  XPRESSMP_CCALL(XPRSgetstringcontrol(lp(), key, s.data(), l, &l));
  value.assign(s.data());
}

void XpressmpCommon::SetSolverOption(int  key, const std::string& value) {
  XPRESSMP_CCALL(XPRSsetstrcontrol(lp(), key, value.c_str()));
}




} // namespace mp

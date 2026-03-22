#include "mp/format.h"
#include "coptcommon.h"

namespace mp {

int CoptCommon::getIntAttr(const char* name)  const {
  int value;
  COPT_CCALL(COPT_GetIntAttr(lp(), name, &value));
  return value;
}
double CoptCommon::getDblAttr(const char* name) const  {
  double value;
  COPT_CCALL(COPT_GetDblAttr(lp(), name, &value));
  return value;
}

void CoptCommon::setIntAttr(const char* name, int value)  {
    COPT_CCALL(COPT_SetIntParam(lp(), name, value));
}
void CoptCommon::setDblAttr(const char* name, double value)  {
    COPT_CCALL(COPT_SetDblParam(lp(), name, value));
}

std::vector<double>  CoptCommon::getVarInfo(const char* name) {
  std::vector<double> ret(NumVars());
  COPT_CCALL(COPT_GetColInfo(lp(), name, NumVars(), NULL, ret.data()));
  return ret;
}
std::vector<double>  CoptCommon::getConInfo(const char* name) {
    std::vector<double> ret(NumLinCons());
    COPT_CCALL(COPT_GetRowInfo(lp(), name, NumLinCons(), NULL, ret.data()));
    return ret;
}

int CoptCommon::NumLinCons() const {
  return getIntAttr(COPT_INTATTR_ROWS);
}

int CoptCommon::NumVars() const {
  return getIntAttr(COPT_INTATTR_COLS);
}

int CoptCommon::NumObjs() const {
    return getIntAttr(COPT_INTATTR_MULTIOBJS);
}

int CoptCommon::NumQPCons() const {
  return getIntAttr(COPT_INTATTR_QCONSTRS);
}

int CoptCommon::NumSOSCons() const {
  return getIntAttr(COPT_INTATTR_SOSS);
}

int CoptCommon::NumIndicatorCons() const {
  return getIntAttr(COPT_INTATTR_INDICATORS);
}

void CoptCommon::GetSolverOption(const char* key, int &value) const {
  COPT_CCALL( COPT_GetIntParam(lp(), key, &value) );
}

void CoptCommon::SetSolverOption(const char* key, int value) {
  if (current_objective_options() == -1)
    COPT_CCALL(COPT_SetIntParam(lp(), key, value));
  else
     COPT_CCALL(COPT_MultiObjSetIntParam(lp(), current_objective_options(), key, value));
}

void CoptCommon::GetSolverOption(const char* key, double &value) const {
  COPT_CCALL(COPT_GetDblParam(lp(), key, &value) );
}

void CoptCommon::SetSolverOption(const char* key, double value) {
    if (current_objective_options() == -1)
        COPT_CCALL(COPT_SetDblParam(lp(), key, value));
    else
        COPT_CCALL(COPT_MultiObjSetDblParam(lp(), current_objective_options(), key, value));
}

void CoptCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CoptCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

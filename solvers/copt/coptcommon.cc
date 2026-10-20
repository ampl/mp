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

void CoptCommon::CheckParamRange(const char* key, int value) const {
  int lo, hi;
  if (COPT_GetIntParamMin(lp(), key, &lo) || COPT_GetIntParamMax(lp(), key, &hi))
    return;                      // unknown parameter: let the setter report
  if (value < lo || value > hi)
    MP_RAISE(fmt::format("Rejecting {}={}; must be between {} and {}",
                         key, value, lo, hi));
}

void CoptCommon::CheckParamRange(const char* key, double value) const {
  double lo, hi;
  if (COPT_GetDblParamMin(lp(), key, &lo) || COPT_GetDblParamMax(lp(), key, &hi))
    return;
  if (value < lo || value > hi)
    MP_RAISE(fmt::format("Rejecting {}={}; must be between {} and {}",
                         key, value, lo, hi));
}

void CoptCommon::SetSolverOption(const char* key, int value) {
  CheckParamRange(key, value);
  if (current_objective_options() == -1)
    COPT_CCALL(COPT_SetIntParam(lp(), key, value));
  else
     COPT_CCALL(COPT_MultiObjSetIntParam(lp(), current_objective_options(), key, value));
}

void CoptCommon::GetSolverOption(const char* key, double &value) const {
  COPT_CCALL(COPT_GetDblParam(lp(), key, &value) );
}

void CoptCommon::SetSolverOption(const char* key, double value) {
    CheckParamRange(key, value);
    if (current_objective_options() == -1)
        COPT_CCALL(COPT_SetDblParam(lp(), key, value));
    else
        COPT_CCALL(COPT_MultiObjSetDblParam(lp(), current_objective_options(), key, value));
}

void CoptCommon::SetNativeParam(const std::string& name_value) {
  auto pos = name_value.find_first_of(" =\t");
  if (pos == std::string::npos)
    MP_RAISE(fmt::format("tech:optionnative: expected 'Name value', got '{}'",
                         name_value));
  std::string name = name_value.substr(0, pos);
  auto vpos = name_value.find_first_not_of(" =\t", pos);
  if (vpos == std::string::npos)
    MP_RAISE(fmt::format("tech:optionnative: no value for parameter '{}'", name));
  std::string value = name_value.substr(vpos);
  while (!value.empty() && (value.back()==' ' || value.back()=='\t'))
    value.pop_back();
  int idef; double ddef;
  if (!COPT_GetIntParamDef(lp(), name.c_str(), &idef)) {
    size_t n = 0;
    int v = 0;
    try { v = std::stoi(value, &n); } catch (...) { n = 0; }
    if (!n || n != value.size())
      MP_RAISE(fmt::format("tech:optionnative: integer parameter '{}' "
                           "cannot take value '{}'", name, value));
    SetSolverOption(name.c_str(), v);
  } else if (!COPT_GetDblParamDef(lp(), name.c_str(), &ddef)) {
    size_t n = 0;
    double v = 0.0;
    try { v = std::stod(value, &n); } catch (...) { n = 0; }
    if (!n || n != value.size())
      MP_RAISE(fmt::format("tech:optionnative: double parameter '{}' "
                           "cannot take value '{}'", name, value));
    SetSolverOption(name.c_str(), v);
  } else
    MP_RAISE(fmt::format("tech:optionnative: unknown COPT parameter '{}'", name));
}

void CoptCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CoptCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

#include "mp/format.h"
#include "ortoolsmpcommon.h"
#include <string>

namespace mp {

int OrtoolsCommon::getIntAttr(int name)  const {
  int value = 0;
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  ORTOOLS_CCALL(ORTOOLS_GetIntAttr(lp_, name, &value)); */
  return value;
}
double OrtoolsCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  ORTOOLS_CCALL(ORTOOLS_GetDblAttr(lp_, name, &value)); */
  return value;
}

int OrtoolsCommon::NumLinCons() const {
  return lp()->NumConstraints();
}

int OrtoolsCommon::NumVars() const {
  return lp()->NumVariables();
}

int OrtoolsCommon::NumObjs() const {
  return 1;
}

int OrtoolsCommon::NumQPCons() const {
  return 0;
}

int OrtoolsCommon::NumSOSCons() const {
  return 0;
}

int OrtoolsCommon::NumIndicatorCons() const {
  return 0;
}

void OrtoolsCommon::GetSolverOption(const char* key, int &value) const {
  if (paramNames_.find(key) == paramNames_.end())
    value = optionsManager_.getIntValue(key);
  else
    value = params_.GetIntegerParam(
      static_cast<orr::MPSolverParameters::IntegerParam>(paramNames_.at(key)));
}


void OrtoolsCommon::SetSolverOption(const char* key, int value) {
  if (paramNames_.find(key) == paramNames_.end())
    optionsManager_.set(key, value);
  else
    params_.SetIntegerParam(
      static_cast<orr::MPSolverParameters::IntegerParam>(paramNames_.at(key)),
      value);
}
void OrtoolsCommon::GetSolverOption(const char* key, double &value) const {
  if (paramNames_.find(key) == paramNames_.end())
    value = optionsManager_.getDoubleValue(key);
  else
    value = params_.GetDoubleParam(
      static_cast<orr::MPSolverParameters::DoubleParam>(paramNames_.at(key)));
}

void OrtoolsCommon::SetSolverOption(const char* key, double value) {
  if (paramNames_.find(key) == paramNames_.end())
    optionsManager_.set(key, value);
  else
    params_.SetDoubleParam(
      static_cast<orr::MPSolverParameters::DoubleParam>(paramNames_[key]),
      value);
}

void OrtoolsCommon::GetSolverOption(const char* key, std::string &value) const {
  value = optionsManager_.getStringValue(key);
}

void OrtoolsCommon::SetSolverOption(const char* key, const std::string& value) {
    optionsManager_.set(key, value);
}


} // namespace mp

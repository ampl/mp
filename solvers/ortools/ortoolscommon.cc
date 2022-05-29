#include "mp/format.h"
#include "ortoolscommon.h"

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
  //ORTOOLS_CCALL( ORTOOLS_GetIntParam(lp_, key, &value) );
}

void OrtoolsCommon::SetSolverOption(const char* key, int value) {
  //ORTOOLS_CCALL(ORTOOLS_SetIntParam(lp_, key, value));
}

void OrtoolsCommon::GetSolverOption(const char* key, double &value) const {
  //ORTOOLS_CCALL(ORTOOLS_GetDblParam(lp_, key, &value) );
}

void OrtoolsCommon::SetSolverOption(const char* key, double value) {
 // ORTOOLS_CCALL(ORTOOLS_SetDblParam(lp_, key, value) );
}

void OrtoolsCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void OrtoolsCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

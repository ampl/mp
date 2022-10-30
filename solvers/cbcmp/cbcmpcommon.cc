#include "mp/format.h"
#include "cbcmpcommon.h"


namespace mp {

void CbcmpCommon::GetCBCParamsList() const {
  
  for (auto p : lp()->cbcData->parameters_)
  {
  //   p.printString();
    fmt::print(p.name());

    if (p.type() >= CLP_PARAM_DBL_PRIMALTOLERANCE &&
      p.type() <= CBC_PARAM_DBL_DEXTRA5)
      fmt::print(" (double)");

    if (p.type() >= CLP_PARAM_INT_SOLVERLOGLEVEL &&
      p.type() <= CBC_PARAM_INT_MOREMOREMIPOPTIONS)
      fmt::print("(int)");
    if (p.type() >= CLP_PARAM_STR_DIRECTION &&
      p.type() <= CBC_PARAM_STR_SOSPRIORITIZE)
      fmt::print("(str)");

    fmt::print(": {}\n", p.longHelp());
  }
}

int CbcmpCommon::getIntAttr(int name)  const {
  int value = 0;
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  CBCMP_CCALL(CBCMP_GetIntAttr(lp_, name, &value)); */
  return 0;
}
double CbcmpCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  CBCMP_CCALL(CBCMP_GetDblAttr(lp_, name, &value)); */
  return value;
}

int CbcmpCommon::NumLinCons() const {
  return Cbc_getNumRows(lp());
}

int CbcmpCommon::NumVars() const {
  return Cbc_getNumCols(lp());
}

int CbcmpCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return CBCMPgetnumobjs (env_, lp_);
  return 0;
}

int CbcmpCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(CBCMP_INTATTR_QCONSTRS);
  return 0;
}

int CbcmpCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(CBCMP_INTATTR_SOSS);
  return 0;
}

int CbcmpCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(CBCMP_INTATTR_INDICATORS);
  return 0;
}

void CbcmpCommon::GetSolverOption(const char* key, int &value) const {
  //CBCMP_CCALL( CBCMP_GetIntParam(lp_, key, &value) );
}

void CbcmpCommon::SetSolverOption(const char* key, int value) {
  std::string s = fmt::format("{}", value);
  Cbc_setParameter(lp(), key, s.data());
}

void CbcmpCommon::GetSolverOption(const char* key, double &value) const {
  //CBCMP_CCALL(CBCMP_GetDblParam(lp_, key, &value) );
}

void CbcmpCommon::SetSolverOption(const char* key, double value) {
  std::string s = fmt::format("{}", value);
  Cbc_setParameter(lp(), key, s.data());
}

void CbcmpCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CbcmpCommon::SetSolverOption(const char* key, const std::string& value) {
  Cbc_setParameter(lp(), key, value.data());
}


} // namespace mp

#include "mp/format.h"
#include "cuoptlpcommon.h"

namespace Solver {
  SolverModel* CreateSolverModel() {
    return new SolverModel();
  }
}
namespace mp {

int CuoptlpCommon::getIntAttr(Solver::ATTRIBS name, Solver::ConsType subtype)  const {
  int value = 0;
  return lp()->GetAttribute(name, subtype);
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  CUOPTLP_CCALL(CUOPTLP_GetIntAttr(lp_, name, &value)); */
  return value;
}
double CuoptlpCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  CUOPTLP_CCALL(CUOPTLP_GetDblAttr(lp_, name, &value)); */
  return value;
}

int CuoptlpCommon::NumLinCons() const {
  return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_LIN);
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(CUOPTLP_INTATTR_ROWS);
}

int CuoptlpCommon::NumVars() const {
  return getIntAttr(Solver::NVARS_CONT)+ getIntAttr(Solver::NVARS_INT);
  // TODO Get number of vars using solver API
  //  return getIntAttr(CUOPTLP_INTATTR_COLS);
}

int CuoptlpCommon::NumObjs() const {
  return getIntAttr(Solver::NOBJS);
  // TODO Get number of objectives using solver API
  //return CUOPTLPgetnumobjs (env_, lp_);
}

int CuoptlpCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(CUOPTLP_INTATTR_QCONSTRS);
  return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_QUAD);
}

void CuoptlpCommon::GetSolverOption(const char* key, int &value) const {
  //CUOPTLP_CCALL( CUOPTLP_GetIntParam(lp_, key, &value) );
}

void CuoptlpCommon::SetSolverOption(const char* key, int value) {
  //CUOPTLP_CCALL(CUOPTLP_SetIntParam(lp_, key, value));
}

void CuoptlpCommon::GetSolverOption(const char* key, double &value) const {
  //CUOPTLP_CCALL(CUOPTLP_GetDblParam(lp_, key, &value) );
}

void CuoptlpCommon::SetSolverOption(const char* key, double value) {
 // CUOPTLP_CCALL(CUOPTLP_SetDblParam(lp_, key, value) );
}

void CuoptlpCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CuoptlpCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

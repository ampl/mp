#include "mp/format.h"
#include "cuoptcommon.h"

namespace Solver {
  SolverModel* CreateSolverModel() {
    return new SolverModel();
  }
}
namespace mp {

int CuoptCommon::getIntAttr(Solver::ATTRIBS name, Solver::ConsType subtype)  const {
  int value = 0;
  return 0;//lp()->GetAttribute(name, subtype);
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  CUOPT_CCALL(CUOPT_GetIntAttr(lp_, name, &value)); */
  return value;
}
double CuoptCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  CUOPT_CCALL(CUOPT_GetDblAttr(lp_, name, &value)); */
  return value;
}

int CuoptCommon::NumLinCons() const {
  return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_LIN);
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(CUOPT_INTATTR_ROWS);
}

int CuoptCommon::NumVars() const {
  return getIntAttr(Solver::NVARS_CONT)+ getIntAttr(Solver::NVARS_INT);
  // TODO Get number of vars using solver API
  //  return getIntAttr(CUOPT_INTATTR_COLS);
}

int CuoptCommon::NumObjs() const {
  return getIntAttr(Solver::NOBJS);
  // TODO Get number of objectives using solver API
  //return CUOPTgetnumobjs (env_, lp_);
}

int CuoptCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(CUOPT_INTATTR_QCONSTRS);
  return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_QUAD);
}

void CuoptCommon::GetSolverOption(const char* key, int &value) const {
  //CUOPT_CCALL( CUOPT_GetIntParam(lp_, key, &value) );
}

void CuoptCommon::SetSolverOption(const char* key, int value) {
  //CUOPT_CCALL(CUOPT_SetIntParam(lp_, key, value));
}

void CuoptCommon::GetSolverOption(const char* key, double &value) const {
  //CUOPT_CCALL(CUOPT_GetDblParam(lp_, key, &value) );
}

void CuoptCommon::SetSolverOption(const char* key, double value) {
 // CUOPT_CCALL(CUOPT_SetDblParam(lp_, key, value) );
}

void CuoptCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CuoptCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

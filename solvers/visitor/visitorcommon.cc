#include "mp/format.h"
#include "visitorcommon.h"

namespace Solver {
  SolverModel* CreateSolverModel() {
    return new SolverModel();
  }
}
namespace mp {

int VisitorCommon::getIntAttr(Solver::ATTRIBS name, Solver::TYPE subtype)  const {
  int value = 0;
  return lp()->getAttribute(name, subtype);
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  VISITOR_CCALL(VISITOR_GetIntAttr(lp_, name, &value)); */
  return value;
}
double VisitorCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  VISITOR_CCALL(VISITOR_GetDblAttr(lp_, name, &value)); */
  return value;
}

int VisitorCommon::NumLinCons() const {
  return getIntAttr(Solver::NCONS_TYPE, Solver::CONS_LIN);
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_ROWS);
}

int VisitorCommon::NumVars() const {
  return getIntAttr(Solver::NVARS_CONT)+ getIntAttr(Solver::NVARS_INT);
  // TODO Get number of vars using solver API
  //  return getIntAttr(VISITOR_INTATTR_COLS);
}

int VisitorCommon::NumObjs() const {
  return getIntAttr(Solver::NOBJS);
  // TODO Get number of objectives using solver API
  //return VISITORgetnumobjs (env_, lp_);
}

int VisitorCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_QCONSTRS);
  return getIntAttr(Solver::NCONS_TYPE, Solver::CONS_QUAD);
}

int VisitorCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_SOSS);
  return getIntAttr(Solver::NCONS_TYPE, Solver::CONS_SOS);
}

int VisitorCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_INDICATORS);
  return getIntAttr(Solver::NCONS_TYPE, Solver::CONS_INDIC);
}

void VisitorCommon::GetSolverOption(const char* key, int &value) const {
  //VISITOR_CCALL( VISITOR_GetIntParam(lp_, key, &value) );
}

void VisitorCommon::SetSolverOption(const char* key, int value) {
  //VISITOR_CCALL(VISITOR_SetIntParam(lp_, key, value));
}

void VisitorCommon::GetSolverOption(const char* key, double &value) const {
  //VISITOR_CCALL(VISITOR_GetDblParam(lp_, key, &value) );
}

void VisitorCommon::SetSolverOption(const char* key, double value) {
 // VISITOR_CCALL(VISITOR_SetDblParam(lp_, key, value) );
}

void VisitorCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void VisitorCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

#include "mp/format.h"
#include "visitorcommon.h"

namespace Solver {

  SolverModel* CreateSolverModel() {
    return new SolverModel();
  }
  int addContVar(SolverModel& s) {
    return s.addEntity(Solver::VARS_INT);
  }
  int addIntVar(SolverModel& s) {
    return s.addEntity(Solver::VARS);
  }
  int addLinCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_LIN);
  }
  int addQuadCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_QUAD);
  }
  int addIndicCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_INDIC);
  }
  int addSOSCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_SOS);
  }
  int getSolverIntAttr(SolverModel* s, int attr) {
    return s->getN(attr);
  }
}
namespace mp {

int VisitorCommon::getIntAttr(int name)  const {
  int value = 0;
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  VISITOR_CCALL(VISITOR_GetIntAttr(lp_, name, &value)); */
  return getSolverIntAttr(lp(), name);
}
double VisitorCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  VISITOR_CCALL(VISITOR_GetDblAttr(lp_, name, &value)); */
  return value;
}

int VisitorCommon::NumLinCons() const {
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_ROWS);
  return 0;
}

int VisitorCommon::NumVars() const {
  // TODO Get number of linear constraints using solver API
  //  return getIntAttr(VISITOR_INTATTR_COLS);
  return 0;
}

int VisitorCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return VISITORgetnumobjs (env_, lp_);
  return 0;
}

int VisitorCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_QCONSTRS);
  return 0;
}

int VisitorCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_SOSS);
  return 0;
}

int VisitorCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(VISITOR_INTATTR_INDICATORS);
  return 0;
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

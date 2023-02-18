#include "mp/format.h"
#include "scipcommon.h"

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

int ScipCommon::getIntAttr(int name)  const {
  int value = 0;
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  SCIP_CCALL(SCIP_GetIntAttr(lp_, name, &value)); */
  return getSolverIntAttr(lp(), name);
}
double ScipCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  SCIP_CCALL(SCIP_GetDblAttr(lp_, name, &value)); */
  return value;
}

int ScipCommon::NumLinCons() const {
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(SCIP_INTATTR_ROWS);
  return 0;
}

int ScipCommon::NumVars() const {
  // gets number of active problem variables
  return SCIPgetNVars(lp_);
}

int ScipCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return SCIPgetnumobjs (env_, lp_);
  return 0;
}

int ScipCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(SCIP_INTATTR_QCONSTRS);
  return 0;
}

int ScipCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(SCIP_INTATTR_SOSS);
  return 0;
}

int ScipCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(SCIP_INTATTR_INDICATORS);
  return 0;
}

void ScipCommon::GetSolverOption(const char* key, int &value) const {
  //SCIP_CCALL( SCIP_GetIntParam(lp_, key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, int value) {
  //SCIP_CCALL(SCIP_SetIntParam(lp_, key, value));
}

void ScipCommon::GetSolverOption(const char* key, double &value) const {
  SCIP_CCALL( SCIPgetRealParam(lp_, key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, double value) {
  SCIP_CCALL( SCIPsetRealParam(lp_, key, value) );
}

void ScipCommon::GetSolverOption(const char* key, std::string &value) const {
  char buffer[SCIP_MAXSTRLEN];
  SCIP_CCALL( SCIPsetStringParam(lp_, key, buffer) );
  value = buffer;
}

void ScipCommon::SetSolverOption(const char* key, const std::string& value) {
  SCIP_CCALL( SCIPsetStringParam(lp_, key, value.c_str()) );
}


} // namespace mp

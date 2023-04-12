#include "mp/format.h"
#include "gcgcommon.h"

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

int GcgCommon::getIntAttr(int name)  const {
  int value = 0;
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  GCG_CCALL(GCG_GetIntAttr(lp_, name, &value)); */
  return getSolverIntAttr(lp(), name);
}
double GcgCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  GCG_CCALL(GCG_GetDblAttr(lp_, name, &value)); */
  return value;
}

int GcgCommon::NumLinCons() const {
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(GCG_INTATTR_ROWS);
  return 0;
}

int GcgCommon::NumVars() const {
  // TODO Get number of linear constraints using solver API
  //  return getIntAttr(GCG_INTATTR_COLS);
  return 0;
}

int GcgCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return GCGgetnumobjs (env_, lp_);
  return 0;
}

int GcgCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(GCG_INTATTR_QCONSTRS);
  return 0;
}

int GcgCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(GCG_INTATTR_SOSS);
  return 0;
}

int GcgCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(GCG_INTATTR_INDICATORS);
  return 0;
}

void GcgCommon::GetSolverOption(const char* key, int &value) const {
  //GCG_CCALL( GCG_GetIntParam(lp_, key, &value) );
}

void GcgCommon::SetSolverOption(const char* key, int value) {
  //GCG_CCALL(GCG_SetIntParam(lp_, key, value));
}

void GcgCommon::GetSolverOption(const char* key, double &value) const {
  //GCG_CCALL(GCG_GetDblParam(lp_, key, &value) );
}

void GcgCommon::SetSolverOption(const char* key, double value) {
 // GCG_CCALL(GCG_SetDblParam(lp_, key, value) );
}

void GcgCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void GcgCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

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

void VisitorCommon::OpenSolver() {
  int status = 0;
  // TODO Typically this function creates an instance of the solver environment
  // and an empty model
  void* env_p;
  // Typically try the registered function first (see visitorcommon.h);
  // if not available call the solver's API function directly
  /*
  if (createEnv == nullptr)
    status = VISITOR_CreateEnv(&env_p);
  else
    status = createEnv(&env_p);
    */
  // set_env(env_p); // TODO Set the environment in the VisitorCommon class

  /* Todo catch errors 
  if ( env() == NULL ) {
    // char  errmsg[CPXMESSAGEBUFSIZE]; 
    // CPXgeterrorstring (env(), status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open VISITOR environment.\n{}", status) );
  }
  */

  /* TODO Create problem instance 
  visitor_prob* prob;
  status = VISITOR_CreateProb(env_p, &prob);
 */
  Solver::SolverModel* prob = Solver::CreateSolverModel();
  set_lp(prob); // Assign it
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  /* TODO Typically check call */
  // VISITOR_CCALL(VISITOR_SetIntParam(prob, "Logging", 0));

}

void VisitorCommon::CloseSolver() {
  /* TODO Cleanup: close problem and environment 
  if ( lp() != NULL ) {
    VISITOR_CCALL(VISITOR_DeleteProb(&lp_) );
  }
  if ( env() != NULL ) {
    VISITOR_CCALL(VISITOR_DeleteEnv(&env_) );
  }
  */
}

void VisitorCommon::copy_handlers_from_other_visitor() {
  assert(other_visitor());
  /* TODO Implement the following
  env_ = other_visitor()->env();
  lp_ = other_visitor()->lp(); */
}

void VisitorCommon::copy_handlers_to_other_visitor() {
  assert(other_visitor());
  /* TODO Implement the following
  other_visitor()->set_env(env_);
  other_visitor()->set_lp(lp_); */
}




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

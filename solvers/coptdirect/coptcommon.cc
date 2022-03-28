#include "mp/format.h"
#include "coptcommon.h"

namespace mp {

void CoptCommon::OpenSolver() {
  int status = 0;
  copt_env* env_p;
  if (createEnv == nullptr)
    status = COPT_CreateEnv(&env_p);
  else
    status = createEnv(&env_p);
  set_env(env_p);
  if ( env() == NULL ) {
    // char  errmsg[CPXMESSAGEBUFSIZE]; TODO
    // CPXgeterrorstring (env(), status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open COPT environment.\n{}", status) );
  }
  copt_prob* prob;
  status = COPT_CreateProb(env_p, &prob);
 
  /* Create an empty model */
  set_lp(prob);
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  COPT_CCALL(COPT_SetIntParam(prob, "Logging", 0));

}

void CoptCommon::CloseSolver() {
  if ( lp() != NULL ) {
    COPT_CCALL(COPT_DeleteProb(&lp_) );
  }
  /* Free up the COPT env()ironment, if necessary */
  if ( env() != NULL ) {
    COPT_CCALL(COPT_DeleteEnv(&env_) );
  }
}

void CoptCommon::copy_handlers_from_other_copt() {
  assert(other_copt());
  env_ = other_copt()->env();
  lp_ = other_copt()->lp();
}

void CoptCommon::copy_handlers_to_other_copt() {
  assert(other_copt());
  other_copt()->set_env(env_);
  other_copt()->set_lp(lp_);
}




int CoptCommon::getIntAttr(const char* name)  const {
  int value;
  COPT_CCALL(COPT_GetIntAttr(lp_, name, &value));
  return value;
}
double CoptCommon::getDblAttr(const char* name) const  {
  double value;
  COPT_CCALL(COPT_GetDblAttr(lp_, name, &value));
  return value;
}

int CoptCommon::NumLinCons() const {
  return getIntAttr(COPT_INTATTR_ROWS);
}

int CoptCommon::NumVars() const {
  return getIntAttr(COPT_INTATTR_COLS);
}

int CoptCommon::NumObjs() const {
  return 1;
  //return CPXgetnumobjs (env_, lp_);
}


void CoptCommon::GetSolverOption(const char* key, int &value) const {
  COPT_CCALL( COPT_GetIntParam(lp_, key, &value) );
}

void CoptCommon::SetSolverOption(const char* key, int value) {
  COPT_CCALL(COPT_SetIntParam(lp_, key, value));
}

void CoptCommon::GetSolverOption(const char* key, double &value) const {
  COPT_CCALL(COPT_GetDblParam(lp_, key, &value) );
}

void CoptCommon::SetSolverOption(const char* key, double value) {
  COPT_CCALL(COPT_SetDblParam(lp_, key, value) );
}

void CoptCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CoptCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

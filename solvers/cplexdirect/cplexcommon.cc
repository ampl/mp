#include "mp/format.h"
#include "cplexcommon.h"

namespace mp {

void CplexCommon::OpenSolver() {
  int status;
  set_env( CPXopenCPLEX (&status) );
  if ( env() == NULL ) {
     char  errmsg[CPXMESSAGEBUFSIZE];
     CPXgeterrorstring (env(), status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open CPLEX environment.\n{}", errmsg ) );
  }

  CPLEX_CALL( CPXsetintparam (env(), CPXPARAM_ScreenOutput, CPX_ON) );

  /* Create an empty model */
  set_lp( CPXcreateprob (env(), &status, "amplcplex") );
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
}

void CplexCommon::CloseSolver() {
  if ( lp() != NULL ) {
     CPLEX_CALL( CPXfreeprob (env(), &lp_) );
  }
  /* Free up the CPLEX env()ironment, if necessary */
  if ( env() != NULL ) {
     CPLEX_CALL( CPXcloseCPLEX (&env_) );
  }
}

void CplexCommon::copy_handlers_from_other_cplex() {
  assert(other_cplex());
  env_ = other_cplex()->env();
  lp_ = other_cplex()->lp();
}

void CplexCommon::copy_handlers_to_other_cplex() {
  assert(other_cplex());
  other_cplex()->set_env(env_);
  other_cplex()->set_lp(lp_);
}


int CplexCommon::NumLinCons() const {
  return CPXgetnumrows (env_, lp_);
}

int CplexCommon::NumVars() const {
  return CPXgetnumcols (env_, lp_);
}

int CplexCommon::NumObjs() const {
  return CPXgetnumobjs (env_, lp_);
}


void CplexCommon::GetSolverOption(int key, int &value) const {
  CPLEX_CALL( CPXgetintparam(env_, key, &value) );
}

void CplexCommon::SetSolverOption(int key, int value) {
  CPLEX_CALL( CPXsetintparam(env_, key, value) );
}

void CplexCommon::GetSolverOption(int key, double &value) const {
  CPLEX_CALL( CPXgetdblparam(env_, key, &value) );
}

void CplexCommon::SetSolverOption(int key, double value) {
  CPLEX_CALL( CPXsetdblparam(env_, key, value) );
}

void CplexCommon::GetSolverOption(int key, std::string &value) const {
  char buffer[CPX_STR_PARAM_MAX];
  CPLEX_CALL( CPXgetstrparam(env_, key, buffer) );
  value = buffer;
}

void CplexCommon::SetSolverOption(int key, const std::string& value) {
  CPLEX_CALL( CPXsetstrparam(env_, key, value.c_str()) );
}


} // namespace mp

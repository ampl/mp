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
     CPLEX_CALL( CPXfreeprob (env(), &lp_ref()) );
  }
  /* Free up the CPLEX env()ironment, if necessary */
  if ( env() != NULL ) {
     CPLEX_CALL( CPXcloseCPLEX (&env_ref()) );
  }
}


int CplexCommon::NumLinCons() const {
  return CPXgetnumrows (env(), lp());
}

int CplexCommon::NumVars() const {
  return CPXgetnumcols (env(), lp());
}

int CplexCommon::NumObjs() const {
  return CPXgetnumobjs (env(), lp());
}


void CplexCommon::GetSolverOption(int key, int &value) const {
  CPLEX_CALL( CPXgetintparam(env(), key, &value) );
}

void CplexCommon::SetSolverOption(int key, int value) {
  CPLEX_CALL( CPXsetintparam(env(), key, value) );
}

void CplexCommon::GetSolverOption(int key, double &value) const {
  CPLEX_CALL( CPXgetdblparam(env(), key, &value) );
}

void CplexCommon::SetSolverOption(int key, double value) {
  CPLEX_CALL( CPXsetdblparam(env(), key, value) );
}

void CplexCommon::GetSolverOption(int key, std::string &value) const {
  char buffer[CPX_STR_PARAM_MAX];
  CPLEX_CALL( CPXgetstrparam(env(), key, buffer) );
  value = buffer;
}

void CplexCommon::SetSolverOption(int key, const std::string& value) {
  CPLEX_CALL( CPXsetstrparam(env(), key, value.c_str()) );
}


} // namespace mp

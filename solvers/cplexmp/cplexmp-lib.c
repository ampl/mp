#include "cplexmp/cplexmp-ampls-c-api.h"


AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLloadmodel(int argc, char** argv, CCallbacks cb) {
  const char* nl_filename = argv[1];
  const char* slv_opt = NULL;
  AMPLS_MP_Solver* slv= AMPLSOpenCPLEX(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename);
  return slv;
}

AMPLS_C_EXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseCPLEX(slv);
}

AMPLS_C_EXPORT CPXLPptr AMPLgetCPLEXModel(AMPLS_MP_Solver* slv) {
  return GetCPLEXmodel(slv);
}

AMPLS_C_EXPORT CPXENVptr AMPLgetCPLEXEnv(AMPLS_MP_Solver* slv) {
  return GetCPLEXenv(slv);
}
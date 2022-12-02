#include "cbcmp/cbcmp-ampls-c-api.h"


AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLloadmodel(int argc, char** argv, CCallbacks cb) {
  const char* nl_filename = argv[1];
  const char* slv_opt = NULL;
  AMPLS_MP_Solver* slv = AMPLSOpenCbcmp(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename);
  return slv;
}

AMPLS_C_EXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseCbcmp(slv);
}
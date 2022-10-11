#include "cplexmp/cplexmp-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

extern CCallbacks getCb();

APIEXPORT CPXLPptr AMPLloadmodel(int argc, char** argv, void** slvout) {
  const char* nl_filename = argv[1];
  const char *slv_opt= argv[2];
  CCallbacks cb = getCb();
  AMPLS_MP_Solver* slv= AMPLSOpenCPLEX(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename);
  return slv;
}

APIEXPORT void AMPLwritesolution(AMPLS_MP_Solver* slv, const char* solFileName) {
  AMPLSReportResults(slv, solFileName);
}

APIEXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseCPLEX(slv);
}

APIEXPORT CPXLPptr AMPLgetCPLEXModel(void* slv) {
  return GetCPLEXmodel((AMPLS_MP_Solver*)slv);
}

APIEXPORT CPXENVptr AMPLgetCPLEXEnv(void* slv) {
  return GetCPLEXenv((AMPLS_MP_Solver*)slv);
}
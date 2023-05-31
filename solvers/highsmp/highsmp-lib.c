
#include "interfaces/highs_c_api.h"

#include "highsmp/highsmp-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

APIEXPORT void* AMPLloadmodel(int argc, char** argv, void** slvout) {
  const char* nl_filename = argv[1];
  const char *slv_opt= NULL;
  CCallbacks cb = { NULL,NULL, NULL, NULL };
  AMPLS_MP_Solver* slv = AMPLSOpenHighs(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename, (char**)0);
  return slv;
}

APIEXPORT void* AMPLgetHighsModel(void* slv) {
  return GetHighsmodel(slv);
}

APIEXPORT void AMPLwritesolution(AMPLS_MP_Solver* slv, const char* solFileName) {
  AMPLSReportResults(slv, solFileName);
}

APIEXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseHighs(slv);
}

#include "xpress/xpress-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

APIEXPORT void* AMPLloadmodel(int argc, char** argv, CCallbacks cb) {
  const char* nl_filename = argv[1];
  const char* slv_opt = NULL; // TODO
  AMPLS_MP_Solver* slv = AMPLSOpenXpressmp(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename);
  return slv;
}

APIEXPORT void AMPLwritesolution(AMPLS_MP_Solver* slv, const char* solFileName) {
  AMPLSReportResults(slv, solFileName);
}

APIEXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseXpressmp(slv);
}
APIEXPORT XPRSprob AMPLgetXpressModel(void* slv) {
  return GetXpressmpmodel((AMPLS_MP_Solver*)slv);
}
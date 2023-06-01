#include "scipmp/scipmp-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif


APIEXPORT void* AMPLloadmodel(int argc, char** argv, CCallbacks cb) {
  const char* nl_filename = argv[1];
  const char *slv_opt= argv[2];
  AMPLS_MP_Solver* slv;
  slv = AMPLSOpenScip(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename, (char**)0);
  return slv;
}
APIEXPORT void* AMPLgetScipmodel(void* slv) {
  return GetScipmodel(slv);
}

APIEXPORT void AMPLwritesolution(AMPLS_MP_Solver* slv, const char* solFileName) {
  AMPLSReportResults(slv, solFileName);
}

APIEXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseScip(slv);
}
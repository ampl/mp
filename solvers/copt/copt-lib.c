#include "copt/copt-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

extern CCallbacks getCB(char**);

APIEXPORT  void* AMPLloadmodel(int argc, char** argv) {
  const char* nl_filename = argv[1];
  const char* slv_opt = NULL;
  CCallbacks cb = getCB(argv);
  AMPLS_MP_Solver* slv = AMPLSOpenCopt(slv_opt, cb);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename);
  return slv;
}
APIEXPORT copt_prob* AMPLgetCoptmodel(void* slv) {
  return GetCoptmodel(slv);
}
APIEXPORT void AMPLwritesolution(AMPLS_MP_Solver* slv, const char* solFileName) {
  AMPLSReportResults(slv, solFileName);
}

APIEXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseCopt(slv);
}
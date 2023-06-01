#include "ortoolsmp-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

APIEXPORT void* AMPLloadmodel(int argc, char** argv, void* slvout) {
  const char* nl_filename = argv[1];
  const char *slv_opt= argv[2];
  AMPLS_MP_Solver *slv = AMPLSOpenOrtools(slv_opt);
  int ret = -1;
  ret = AMPLSLoadNLModel(slv, nl_filename, (char**)0);
  void* mdl = GetOrtoolsmodel(slv);
  slvout = &slv;
  return mdl;
}

APIEXPORT void AMPLwritesolution(AMPLS_MP_Solver* slv) {
  AMPLSReportResults(slv);
}

APIEXPORT void AMPLclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseOrtools(slv);
}
#include "gurobi/gurobi-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif


APIEXPORT GRBmodel* AMPLdirectloadmodel(
    int argc, char** argv, void** slvout) {
  const char* nl_filename = argv[1];
  const char *slv_opt= NULL;
  AMPLS_MP_Solver *slv = AMPLSOpenGurobi(slv_opt);
  if (!slv)
    return NULL;
  AMPLSLoadNLModel(slv, nl_filename);
  GRBmodel* mdl = GetGRBmodel(slv);
  *slvout = slv;
  return mdl;
}

APIEXPORT void AMPLdirectwritesolution(void* slv) {
  AMPLSReportResults(slv);
}

APIEXPORT void AMPLdirectclosesolver(void* slv) {
  AMPLSCloseGurobi(slv);
}

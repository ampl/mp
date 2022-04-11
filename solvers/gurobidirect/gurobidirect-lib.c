#include "gurobidirect/gurobi-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

APIEXPORT GRBmodel* AMPLdirectloadmodel(int argc, char** argv, void* slvout) {
  const char* nl_filename = argv[1];
  const char *slv_opt= argv[2];
  AMPLS_MP_Solver slv;
  int ret = -1;
  ret = AMPLSOpenGurobi(&slv, slv_opt);
  ret = AMPLSLoadNLModel(&slv, nl_filename);
  GRBmodel* mdl = GetGRBmodel(&slv);
  slvout = &slv;
  return mdl;
}

APIEXPORT void AMPLdirectwritesolution(AMPLS_MP_Solver* slv) {
  AMPLSReportResults(&slv);
}

APIEXPORT void AMPLdirectclosesolver(AMPLS_MP_Solver* slv) {
  AMPLSCloseGurobi(&slv);
}
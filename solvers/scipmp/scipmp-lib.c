#include "scipmp/scipmp-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_scip(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_scip(cb);
}
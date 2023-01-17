#include "cplexmp/cplexmp-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_cplexmp(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_cplexmp(NULL, cb);
}
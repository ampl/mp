#include "gcgmp/gcgmp-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_gcg(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_gcg(cb);
}
#include "copt/copt-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_copt(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_copt(cb);
}
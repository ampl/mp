#include "cpsat/cpsat-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_cpsat(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_cpsat(cb);
}

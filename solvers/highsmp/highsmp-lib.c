#include "highsmp/highsmp-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_highs(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_highs(cb);
}



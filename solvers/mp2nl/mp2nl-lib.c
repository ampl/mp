#include "mp2nl/mp2nl-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_mp2nl(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_MP2NL(cb);
}

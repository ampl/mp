#include "cbcmp/cbcmp-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_cbcmp(int argc, char** argv)
{
  return Open_cbcmp(NULL);
}
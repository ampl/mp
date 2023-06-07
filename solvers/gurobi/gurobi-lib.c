#include "gurobi/gurobi-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_gurobi(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_gurobi(cb);
}
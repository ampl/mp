#include "xpress/xpress-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_xpress(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_xpress(NULL, cb);
}

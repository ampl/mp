#include "mosek/mosek-ampls-c-api.h"

#ifdef _WIN32
#define APIEXPORT __declspec(dllexport)
#else
#define APIEXPORT  __attribute__((visibility("default")))
#endif

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_mosek(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_mosek(cb);
}

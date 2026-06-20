#include "mp2nl/mp2nl-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_mp2nl(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_MP2NL(cb);
}


#ifdef MP_LINK_WITH_SHARED_LIB
AMPLS_C_EXPORT int mp2nl_main(int argc, char** argv)
{
	return main(argc, argv);
}
#endif

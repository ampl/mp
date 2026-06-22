//#include "interfaces/highs_c_api.h"
#include "highsmp/highsmp-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_highs(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_highs(cb);
}


#ifdef MP_LINK_WITH_SHARED_LIB

int main(int argc, char** argv);

AMPLS_C_EXPORT int highs_main(int argc, char** argv)
{
	return main(argc, argv);
}
#endif

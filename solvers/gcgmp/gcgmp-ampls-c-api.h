#ifndef GCGAMPLSCAPI_H
#define GCGAMPLSCAPI_H
/*
  * C API for MP/Gcg
  */

 //#include "gcg.h"

 #include "mp/ampls-c-api.h"


 /// Initialize AMPLS gcg.

 /// @param slv_opt: a string of solver options
 /// (normally provided in the <solver>_options string).
 /// Can be NULL.
 /// @return pointer to struct AMPLS_MP_Solver to be populated.
 void*  AMPLSOpenGcg(const char* slv_opt, CCallbacks cb);

 /// Shut down solver instance
 void AMPLSCloseGcg(AMPLS_MP_Solver* slv);

 /// Extract the Gcg model handle
 void* GetGcgmodel(AMPLS_MP_Solver* slv);


 #endif // GCGAMPLSCAPI_H
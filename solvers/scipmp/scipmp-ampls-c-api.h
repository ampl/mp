#ifndef SCIPAMPLSCAPI_H
#define SCIPAMPLSCAPI_H
/*
  * C API for MP/Scip
  */

 //#include "scip.h"

 #include "mp/ampls-c-api.h"


 /// Initialize AMPLS scip.

 /// @param slv_opt: a string of solver options
 /// (normally provided in the <solver>_options string).
 /// Can be NULL.
 /// @return pointer to struct AMPLS_MP_Solver to be populated.
 void*  AMPLSOpenScip(const char* slv_opt, CCallbacks cb);

 /// Shut down solver instance
 void AMPLSCloseScip(AMPLS_MP_Solver* slv);

 /// Extract the Scip model handle
 void* GetScipmodel(AMPLS_MP_Solver* slv);


 #endif // SCIPAMPLSCAPI_H

#ifndef KNITROMPAMPLSCAPI_H
#define KNITROMPAMPLSCAPI_H
/*
 * C API for MP/Knitromp
 */

//#include "knitromp.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Knitromp-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS knitromp.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
void*  AMPLSOpenKnitromp(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseKnitromp(AMPLS_MP_Solver* slv);

/// Extract the Knitromp model handle
void* GetKnitrompmodel(AMPLS_MP_Solver* slv);


#endif // KNITROMPAMPLSCAPI_H



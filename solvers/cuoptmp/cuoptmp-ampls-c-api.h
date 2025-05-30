#ifndef CUOPTMPAMPLSCAPI_H
#define CUOPTMPAMPLSCAPI_H
/*
 * C API for MP/cuoptmp
 */

#include "mp/ampls-c-api.h"

/*
 * Below are cuoptmp-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS cuoptmp.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
void*  AMPLSOpenCuoptmp(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseCuoptmp(AMPLS_MP_Solver* slv);

/// Extract the cuoptmp model handle
void* GetCuoptmpmodel(AMPLS_MP_Solver* slv);


#endif // CUOPTMPAMPLSCAPI_H



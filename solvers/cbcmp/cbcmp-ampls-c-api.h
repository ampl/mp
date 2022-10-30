#ifndef CBCMPAMPLSCAPI_H
#define CBCMPAMPLSCAPI_H
/*
 * C API for MP/Cbcmp
 */

//#include "cbcmp.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Cbcmp-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS cbcmp.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
AMPLS_MP_Solver*  AMPLSOpenCbcmp(const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseCbcmp(AMPLS_MP_Solver* slv);

/// Extract the Cbcmp model handle
void* GetCbcmpmodel(AMPLS_MP_Solver* slv);


#endif // CBCMPAMPLSCAPI_H

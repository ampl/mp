#ifndef HIGHSAMPLSCAPI_H
#define HIGHSAMPLSCAPI_H
/*
 * C API for MP/Highs
 */

//#include "highs.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Highs-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS highs.
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
int AMPLSOpenHighs(AMPLS_MP_Solver* slv, const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseHighs(AMPLS_MP_Solver* slv);

/// Extract the Highs model handle
void* GetHighsmodel(AMPLS_MP_Solver* slv);


#endif // HIGHSAMPLSCAPI_H

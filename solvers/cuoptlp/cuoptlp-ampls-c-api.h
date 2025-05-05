#ifndef CUOPTLPAMPLSCAPI_H
#define CUOPTLPAMPLSCAPI_H
/*
 * C API for MP/Cuoptlp
 */

//#include "cuoptlp.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Cuoptlp-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS cuoptlp.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
void*  AMPLSOpenCuoptlp(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseCuoptlp(AMPLS_MP_Solver* slv);

/// Extract the Cuoptlp model handle
void* GetCuoptlpmodel(AMPLS_MP_Solver* slv);


#endif // CUOPTLPAMPLSCAPI_H



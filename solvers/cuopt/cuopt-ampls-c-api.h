#ifndef CUOPTAMPLSCAPI_H
#define CUOPTAMPLSCAPI_H
/*
 * C API for MP/Cuopt
 */

//#include "cuopt.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Cuopt-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS cuopt.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
void*  AMPLSOpenCuopt(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseCuopt(AMPLS_MP_Solver* slv);

/// Extract the Cuopt model handle
void* GetCuoptmodel(AMPLS_MP_Solver* slv);


#endif // CUOPTAMPLSCAPI_H



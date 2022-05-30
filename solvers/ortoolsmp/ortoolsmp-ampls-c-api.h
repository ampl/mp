#ifndef ORTOOLSAMPLSCAPI_H
#define ORTOOLSAMPLSCAPI_H
/*
 * C API for MP/Ortools
 */

#include "mp/ampls-c-api.h"

/*
 * Below are Ortools-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS ortools.
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
int AMPLSOpenOrtools(AMPLS_MP_Solver* slv, const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseOrtools(AMPLS_MP_Solver* slv);

/// Extract the Ortools model handle
void* GetOrtoolsmodel(AMPLS_MP_Solver* slv);


#endif // ORTOOLSAMPLSCAPI_H

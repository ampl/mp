#ifndef COPTAMPLSCAPI_H
#define COPTAMPLSCAPI_H
/*
 * C API for MP/Copt
 */

#include "copt.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Copt-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS copt.
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
int AMPLSOpenCopt(AMPLS_MP_Solver* slv, const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseCopt(AMPLS_MP_Solver* slv);

/// Extract the Copt model handle
copt_prob* GetCoptmodel(AMPLS_MP_Solver* slv);


#endif // COPTAMPLSCAPI_H

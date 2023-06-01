#ifndef GCGAMPLSCAPI_H
#define GCGAMPLSCAPI_H
/*
 * C API for MP/Gcg
 */

#include "gcg/gcg.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Gcg-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS gcg.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see ret_val->warnings_and_or_errors_
AMPLS_MP_Solver* Open_gcg(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
AMPLS_C_EXPORT void AMPLSClose_gcg(AMPLS_MP_Solver* slv);

/// Extract the Gurobi model handle
AMPLS_C_EXPORT SCIP* AMPLSGetModel_gcg(AMPLS_MP_Solver* slv);


#endif // GCGAMPLSCAPI_H



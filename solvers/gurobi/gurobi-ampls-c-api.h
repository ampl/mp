#ifndef GUROBIAMPLSCAPI_H
#define GUROBIAMPLSCAPI_H

/*
 * C API for MP/Gurobi
 */

#include "gurobi_c.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Gurobi-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS Gurobi.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see ret_val->warnings_and_or_errors_
AMPLS_MP_Solver* Open_gurobi(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
AMPLS_C_EXPORT void AMPLSClose_gurobi(AMPLS_MP_Solver* slv);

/// Extract the Gurobi model handle
AMPLS_C_EXPORT GRBmodel* AMPLSGetModel_gurobi(AMPLS_MP_Solver* slv);


#endif // GUROBIAMPLSCAPI_H

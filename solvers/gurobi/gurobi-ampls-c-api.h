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
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
AMPLS_MP_Solver* AMPLSOpenGurobi(const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseGurobi(AMPLS_MP_Solver* slv);

/// Extract the Gurobi model handle
GRBmodel* GetGRBmodel(AMPLS_MP_Solver* slv);


#endif // GUROBIAMPLSCAPI_H

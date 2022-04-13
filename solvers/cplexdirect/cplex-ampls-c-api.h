#ifndef CPLEXAMPLSCAPI_H
#define CPLEXAMPLSCAPI_H

/*
 * C API for MP/Gurobi
 */

#include <ilcplex/cplex.h>

#include "mp/ampls-c-api.h"

/*
 * Below are CPLEX-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS CPLEX.
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
int AMPLSOpenCPLEX(AMPLS_MP_Solver* slv, const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseCPLEX(AMPLS_MP_Solver* slv);

/// Extract the Gurobi model handle
CPXLPptr GetCPLEXmodel(AMPLS_MP_Solver* slv);


#endif // CPLEXAMPLSCAPI_H

#ifndef CPLEXAMPLSCAPI_H
#define CPLEXAMPLSCAPI_H

/*
 * C API for MP/CPLEX
 */

#include <ilcplex/cplex.h>

#include "mp/ampls-c-api.h"

/*
 * Below are CPLEX-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */
DECLARE_SOLVER_API_FUNCTIONS(cplexmp)

AMPLS_C_EXPORT void* AMPLSGetEnv_cplexmp(AMPLS_MP_Solver* slv);

#endif // CPLEXAMPLSCAPI_H

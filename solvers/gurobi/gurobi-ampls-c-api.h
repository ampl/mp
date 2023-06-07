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
DECLARE_SOLVER_API_FUNCTIONS(gurobi)

#endif // GUROBIAMPLSCAPI_H

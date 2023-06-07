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
DECLARE_SOLVER_API_FUNCTIONS(copt)


#endif // COPTAMPLSCAPI_H

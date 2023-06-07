#ifndef MOSEKAMPLSCAPI_H
#define MOSEKAMPLSCAPI_H
/*
 * C API for MP/Mosek
 */

//#include "mosek.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Mosek-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */
DECLARE_SOLVER_API_FUNCTIONS(mosek)

#endif // MOSEKAMPLSCAPI_H

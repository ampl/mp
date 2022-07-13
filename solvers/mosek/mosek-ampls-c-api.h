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

/// Initialize AMPLS mosek.
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
int AMPLSOpenMosek(AMPLS_MP_Solver* slv, const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseMosek(AMPLS_MP_Solver* slv);

/// Extract the Mosek model handle
MSKtask_t GetMosekmodel(AMPLS_MP_Solver* slv);


#endif // MOSEKAMPLSCAPI_H

#ifndef XPRESSMPAMPLSCAPI_H
#define XPRESSMPAMPLSCAPI_H
/*
 * C API for MP/Xpressmp
 */

#include "xprs.h"
#include "mp/ampls-c-api.h"

/*
 * Below are Xpressmp-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS xpressmp.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
AMPLS_MP_Solver* Open_xpress(const char* slv_opt,
  CCallbacks cb);

/// Shut down solver instance
AMPLS_C_EXPORT void AMPLSClose_xpress(AMPLS_MP_Solver* slv);

/// Extract the Xpressmp model handle
AMPLS_C_EXPORT XPRSprob  AMPLSGetModel_xpress(AMPLS_MP_Solver* slv);


#endif // XPRESSMPAMPLSCAPI_H

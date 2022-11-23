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
AMPLS_MP_Solver* AMPLSOpenXpressmp(const char* slv_opt,
  CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseXpressmp(AMPLS_MP_Solver* slv);

/// Extract the Xpressmp model handle
XPRSprob  GetXpressmpmodel(AMPLS_MP_Solver* slv);


#endif // XPRESSMPAMPLSCAPI_H

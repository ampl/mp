#ifndef VISITORAMPLSCAPI_H
#define VISITORAMPLSCAPI_H
/*
 * C API for MP/Visitor
 */

//#include "visitor.h"

#include "mp/ampls-c-api.h"

/*
 * Below are Visitor-specific AMPLS API functions.
 * They complement the 'public' AMPLS API defined in ampls-c-api.h.
 */

/// Initialize AMPLS visitor.

/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return pointer to struct AMPLS_MP_Solver to be populated.
void*  AMPLSOpenVisitor(const char* slv_opt, CCallbacks cb);

/// Shut down solver instance
void AMPLSCloseVisitor(AMPLS_MP_Solver* slv);

/// Extract the Visitor model handle
void* GetVisitormodel(AMPLS_MP_Solver* slv);


#endif // VISITORAMPLSCAPI_H



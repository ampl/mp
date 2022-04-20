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
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options string).
/// Can be NULL.
/// @return 0 on success, otherwise see slv->warnings_and_or_errors_
int AMPLSOpenVisitor(AMPLS_MP_Solver* slv, const char* slv_opt);

/// Shut down solver instance
void AMPLSCloseVisitor(AMPLS_MP_Solver* slv);

/// Extract the Visitor model handle
Solver::SolverModel* GetVisitormodel(AMPLS_MP_Solver* slv);


#endif // VISITORAMPLSCAPI_H

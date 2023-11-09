/**
 * C API: extern "C" wrapper for the NLSOL class
 *
 */

#ifndef NLSOLC_H
#define NLSOLC_H

#include "api/c/nl-feeder2-c.h"
#include "api/c/sol-handler2-c.h"
#include "api/c/nl-writer2-misc-c.h"


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// extern "C" wrapper of mp::NLSOL.
///
/// Manager for solving optimization models via NL files.
/// It performs zero-overhead model/solution transmission.
/// In particular, it does not store any intermediate
/// model representation.
///
/// Usage: see the C API example and the below API.
typedef struct NLSOL_C {
  void* p_;

} NLSOL_C;

//////////// Use the following NLSOL_C API to operate //////////////

/// Construct
NLSOL_C NLW2_MakeNLSOL_C(NLFeeder2_C* , SOLHandler2_C* , NLUtils_C* );

/// Destroy
void NLW2_DestroyNLSOL_C(NLSOL_C* );

/// Set solver, such as "gurobi", "highs", "ipopt"
void NLSOL_C_SetSolver(NLSOL_C* , const char* solver);

/// Set solver options, such as "outlev=1 lim:time=500"
void NLSOL_C_SetSolverOptions(NLSOL_C* , const char* sopts);

/// Solve.
/// @param filestub: filename stub to be used
/// for input files (.nl, .col., .row, etc.),
/// and output files (.sol).
/// @return true if all ok.
int NLSOL_C_Solve(NLSOL_C* , const char* filestub);

/// Get error message.
const char* NLSOL_C_GetErrorMessage(NLSOL_C* );

/// Substep: write NL and any accompanying files.
int NLSOL_C_WriteNLFile(NLSOL_C* , const char* filestub);

/// Substep: invoke chosen solver for \a filestub.
int NLSOL_C_InvokeSolver(NLSOL_C* , const char* filestub);

/// Substep: read solution.
/// @param filename: complete file name,
/// normally (stub).sol.
int NLSOL_C_ReadSolution(NLSOL_C* , const char* filename);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLSOLC_H

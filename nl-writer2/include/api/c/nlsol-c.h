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
/// To create / destroy,
/// use NLW2_MakeNLSOL_C() / NLW2_DestroyNLSOL_C().
///
/// To manipulate, use NLW2_NLSOL_C_SetSolver(),
/// NLW2_NLSOL_C_SetSolverOptions(),
/// etc, see the below API.
///
/// @see the C API example.
typedef struct NLW2_NLSOL_C {
  /// Internal data
  void *p_nlf_, *p_solh_, *p_utl_, *p_nlsol_;

} NLW2_NLSOL_C;

//////////// Use the following NLW2_NLSOL_C API to operate //////////////

/// Construct
NLW2_NLSOL_C NLW2_MakeNLSOL_C(
    NLW2_NLFeeder2_C* , NLW2_SOLHandler2_C* , NLW2_NLUtils_C* );

/// Destroy
void NLW2_DestroyNLSOL_C(NLW2_NLSOL_C* );

/// Set solver, such as "gurobi", "highs", "ipopt"
void NLW2_NLSOL_C_SetSolver(NLW2_NLSOL_C* , const char* solver);

/// Set solver options, such as "outlev=1 lim:time=500"
void NLW2_NLSOL_C_SetSolverOptions(NLW2_NLSOL_C* , const char* sopts);

/// Solve.
/// @param filestub: filename stub to be used
/// for input files (.nl, .col., .row, etc.),
/// and output files (.sol).
/// @return true if all ok.
int NLW2_NLSOL_C_Solve(NLW2_NLSOL_C* , const char* filestub);

/// Get error message.
const char* NLW2_NLSOL_C_GetErrorMessage(NLW2_NLSOL_C* );

/// Substep: write NL and any accompanying files.
/// @return true if all ok.
int NLW2_NLSOL_C_WriteNLFile(NLW2_NLSOL_C* , const char* filestub);

/// Substep: invoke chosen solver for \a filestub.
/// @return true if all ok.
int NLW2_NLSOL_C_InvokeSolver(NLW2_NLSOL_C* , const char* filestub);

/// Substep: read solution.
/// @param filename: complete file name,
///   normally (stub).sol.
/// @return true if all ok.
int NLW2_NLSOL_C_ReadSolution(NLW2_NLSOL_C* , const char* filename);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLSOLC_H

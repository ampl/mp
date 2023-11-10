/**
 * NL Writer C API implementation
 *
 */

#include <cstdlib>
#include <cassert>

#include "api/c/nl-feeder2-c.h"
#include "api/c/sol-handler2-c.h"
#include "api/c/nl-writer2-misc-c.h"
#include "api/c/nlsol-c.h"

#include "mp/nl-writer2-misc.h"

///////////////////////// NLUtils_C ///////////////////////////

NLUtils_C NLW2_MakeNLUtils_C_Default() {
  NLUtils_C result;

  result.p_api_data_ = new mp::NLUtils;

  result.p_user_data_ = NULL;

  return result;
}

void NLW2_DestroyNLUtils_C_Default(NLUtils_C* pu) {
  delete (mp::NLUtils*)pu->p_api_data_;
}


//////////// NLSOL_C API //////////////

/// Construct
NLSOL_C NLW2_MakeNLSOL_C(NLFeeder2_C* , SOLHandler2_C* , NLUtils_C* ) {
  NLSOL_C result;
  return result;
}

/// Destroy
void NLW2_DestroyNLSOL_C(NLSOL_C* ) { assert(0); }

/// Set solver, such as "gurobi", "highs", "ipopt"
void NLW2_SetSolver_C(NLSOL_C* , const char* solver)
{ assert(0); }

/// Set solver options, such as "outlev=1 lim:time=500"
void NLW2_SetSolverOptions_C(NLSOL_C* , const char* sopts)
{ assert(0); }

/// Solve.
/// @param filestub: filename stub to be used
/// for input files (.nl, .col., .row, etc.),
/// and output files (.sol).
/// @return true if all ok.
int NLW2_Solve_C(NLSOL_C* , const char* filestub)
{ assert(0); }

/// Get error message.
const char* NLW2_GetErrorMessage_C(NLSOL_C* )
{ assert(0); }

/// Substep: write NL and any accompanying files.
int NLW2_WriteNLFile_C(NLSOL_C* , const char* filestub)
{ assert(0); }

/// Substep: invoke chosen solver for \a filestub.
int NLW2_InvokeSolver_C(NLSOL_C* , const char* filestub)
{ assert(0); }

/// Substep: read solution.
/// @param filename: complete file name,
/// normally (stub).sol.
int NLW2_ReadSolution_C(NLSOL_C* , const char* filename)
{ assert(0); }

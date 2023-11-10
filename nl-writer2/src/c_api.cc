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

#include "api/c/nl-feeder2-c-impl.h"
#include "api/c/sol-handler2-c-impl.h"
#include "api/c/nl-writer2-misc-c-impl.h"
#include "mp/nlsol.h"

///////////////////////// NLUtils_C ///////////////////////////

NLUtils_C NLW2_MakeNLUtils_C_Default() {
  NLUtils_C result;

  result.p_user_data_ = NULL;

  // TODO set default log/openf/myexit...

  return result;
}

void NLW2_DestroyNLUtils_C_Default(NLUtils_C* ) { }


//////////// NLSOL_C API //////////////

/// Construct.
///
/// Note that the argument objects are stored by value.
NLSOL_C NLW2_MakeNLSOL_C(
    NLFeeder2_C* pnlf, SOLHandler2_C* psolh, NLUtils_C* putl) {
  NLSOL_C result;

  result.p_nlf_ = new mp::NLFeeder2_C_Impl(pnlf);
  result.p_solh_ = new mp::SOLHandler2_C_Impl(psolh);
  result.p_utl_ = new mp::NLUtils_C_Impl(putl);
  result.p_nlsol_
      = new mp::NLSOL<mp::NLFeeder2_C_Impl,
      mp::SOLHandler2_C_Impl>(
        *(mp::NLFeeder2_C_Impl*)result.p_nlf_,
        *(mp::SOLHandler2_C_Impl*)result.p_solh_,
        *(mp::NLUtils_C_Impl*)result.p_utl_);

  return result;
}

/// Destroy
void NLW2_DestroyNLSOL_C(NLSOL_C* pnls) {
  delete (mp::NLUtils_C_Impl*)pnls->p_utl_;
  delete (mp::SOLHandler2_C_Impl*)pnls->p_solh_;
  delete (mp::NLFeeder2_C_Impl*)pnls->p_nlf_;
  delete (mp::NLSOL<mp::NLFeeder2_C_Impl,
      mp::SOLHandler2_C_Impl>*)pnls->p_nlsol_;

  assert(0);
}

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

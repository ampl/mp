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
/// model/solution representation.
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
  void *p_utl_, *p_nlsol_;
} NLW2_NLSOL_C;

//////////// Use the following NLW2_NLSOL_C API to operate //////////////

/// Construct.
/// @param p_utils: utils object.
NLW2_NLSOL_C NLW2_MakeNLSOL_C(NLW2_NLUtils_C* p_utils);

/// Destroy.
void NLW2_DestroyNLSOL_C(NLW2_NLSOL_C* p_nlsol);

/// Set file stub [OPTIONAL].
///
/// Used for filename base of .nl, .col, row, etc. input files,
/// as well as .sol output files.
///
/// If not provided, a temporary filename is used;
/// then, .nl is deleted upon object desruction.
void NLW2_NLSOL_C_SetFileStub(NLW2_NLSOL_C* , const char* stub);

/// Retrieve file stub.
const char* NLW2_NLSOL_C_GetFileStub(NLW2_NLSOL_C* );

/// Get error message.
/// Nonempty iff error occurred.
const char* NLW2_NLSOL_C_GetErrorMessage(NLW2_NLSOL_C* );

/// Write NL and any accompanying files.
///
/// @param nlf_c: NL feeder.
///
/// @return true if all ok, otherwise see
///   GetErrorMessage().
/// @return Non-0 if all ok.
int NLW2_NLSOL_C_LoadModel(NLW2_NLSOL_C* , NLW2_NLFeeder2_C* nlf_c);

/// Solve.
///
/// @param solver: solver executable, such as "gurobi".
/// @param solver_opts: string of solver options,
///   such as "outlev=1 writeprob=model.lp".
///
/// @return Non-0 if all ok.
int NLW2_NLSOL_C_Solve(NLW2_NLSOL_C* ,
                       const char* solver, const char* solver_opts);

/// Read solution.
///
/// @param solh_c: solution handler.
///
/// @return Non-0 if all ok.
int NLW2_NLSOL_C_ReadSolution(NLW2_NLSOL_C* ,
                              NLW2_SOLHandler2_C* solh_c);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLSOLC_H

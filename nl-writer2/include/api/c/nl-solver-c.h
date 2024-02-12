/**
 * C API: extern "C" wrapper for the NLSOL class
 *
 */

#ifndef NLSOLVER_C_H
#define NLSOLVER_C_H

#include "api/c/nl-model-c.h"
#include "api/c/nl-feeder2-c.h"
#include "api/c/sol-handler2-c.h"
#include "api/c/nl-writer2-misc-c.h"


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// extern "C" wrapper of mp::NLSolver.
///
/// Manager for solving optimization models via NL files.
/// It performs zero-overhead model/solution transmission.
/// In particular, it does not store any intermediate
/// model/solution representation.
///
/// This wrapper offers both "easy" but limited,
/// as well as full NL functionality.
///
/// To create / destroy,
/// use NLW2_MakeNLSolver_C() / NLW2_DestroyNLSolver_C().
///
/// To manipulate, use NLW2_LoadModel_C(),
/// etc, see the below API.
///
/// @see C API tests/examples.
typedef struct NLW2_NLSolver_C {
  /// Internal data
  void *p_nlsol_;
  void *p_utl_;
  void *p_sol_;      // for NLModel
} NLW2_NLSolver_C;


//////////// NLSolver C API //////////////

/// Construct.
/// @param p_utils: utils object, can be NULL.
NLW2_NLSolver_C NLW2_MakeNLSolver_C(NLW2_NLUtils_C* p_utils);

/// Destroy NLSolver_C.
void NLW2_DestroyNLSolver_C(NLW2_NLSolver_C* p_nlsol);

/// Set file stub [OPTIONAL].
///
/// Used for filename base of .nl, .col, row, etc. input files,
/// as well as .sol output files.
///
/// @note If not provided, a temporary filename is used;
/// then, .nl is deleted upon object desruction.
void NLW2_SetFileStub_C(NLW2_NLSolver_C* , const char* stub);

/// Retrieve file stub.
const char* NLW2_GetFileStub_C(NLW2_NLSolver_C* );

/// Set NL options for loading NLW2_NLModel_C [OPTIONAL].
///
/// If not provided, default is used.
///
/// @note Not used for loading with NLW2_NLFeeder2_C.
void NLW2_SetNLOptions_C(NLW2_NLSolver_C* ,
                         NLW2_NLOptionsBasic_C );

/// Get NLOptions used for loading NLModel.
NLW2_NLOptionsBasic_C NLW2_GetNLOptions_C(NLW2_NLSolver_C* );

/// Get error message.
/// Nonempty iff error occurred.
const char* NLW2_GetErrorMessage_C(NLW2_NLSolver_C* );

/// Load and solve NLW2_NLModel_C and return result.
///
/// @return Solution object with computed obj_value_.
///   Valid as long as the NLW2_NLSolver_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @see NLW2_LoadNLModel_C(), NLW2_RunSolver_C(),
///   NLW2_ReadSolution_C() for details.
NLW2_Solution_C NLW2_SolveNLModel_C(NLW2_NLSolver_C* ,
                                    NLW2_NLModel_C* ,
                                    const char* solver,
                                    const char* solver_opts);

/// Load NLW2_NLFeeder2_C, solve,
/// and read into SOLHandler2_C.
///
/// @return Non-zero iff all ok.
///
/// @see NLW2_LoadNLFeed2_C(), NLW2_RunSolver_C(),
///   NLW2_Read2SOLHanlder2_C() for details.
int NLW2_SolveFeederHandler_C(NLW2_NLSolver_C* ,
                              NLW2_NLFeeder2_C* ,
                              NLW2_SOLHandler2_C* ,
                              const char* solver,
                              const char* solver_opts);

/// Load NLW2_NLModel_C.
/// NL file name base and some options
/// can be provided, if non-defaults desired,
/// via SetFileStub() and SetNLOptions().
///
/// @return true if all ok, otherwise see
///   GetErrorMessage().
int NLW2_LoadNLModel_C(NLW2_NLSolver_C* ,
                       NLW2_NLModel_C* );

/// Write NL and any accompanying files.
///
/// @param nlf_c: NL feeder.
///
/// @return true if all ok, otherwise see
///   GetErrorMessage().
/// @return Non-0 if all ok.
int NLW2_LoadNLFeed2_C(NLW2_NLSolver_C* , NLW2_NLFeeder2_C* nlf_c);

/// Run a solver after loading model.
///
/// @param solver: solver executable, such as "gurobi".
/// @param solver_opts: string of solver options,
///   such as "outlev=1 writeprob=model.lp".
///
/// @return Non-0 if all ok.
int NLW2_RunSolver_C(NLW2_NLSolver_C* ,
                     const char* solver, const char* solver_opts);

/// Read solution after solving a model
/// loaded with NLW2_LoadNLModel_C().
///
/// @return Solution object.
///   Valid as long as the NLW2_NLSolver_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @note To compute objective value,
///   execute NLW2_ComputeObjValue_C()
///   if x_ available.
NLW2_Solution_C NLW2_ReadSolution_C(NLW2_NLSolver_C* );

/// Read solution to SOLHandler2.
///
/// @param solh_c: solution handler.
///
/// @return Non-0 if all ok.
int NLW2_Read2SOLHandler2_C(NLW2_NLSolver_C* ,
                            NLW2_SOLHandler2_C* solh_c);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLSOLVER_C_H

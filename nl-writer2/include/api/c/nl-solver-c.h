/*
 C API: extern "C" wrapper for mp::NLSolver.

 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */

#ifndef NLSOLVER_C_H
#define NLSOLVER_C_H

#include "api/c/nl-model-c.h"
#include "api/c/nl-feeder-c.h"
#include "api/c/sol-handler-c.h"
#include "api/c/nl-writer2-misc-c.h"


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Manager for solving optimization models via NL files.
/// It performs zero-overhead model/solution transmission.
/// In particular, it does not store any intermediate
/// model/solution representation.
///
/// This wrapper offers both full NL functionality,
/// as well as "easy" but limited interface via NLW2_NLModel_C.
///
/// @note To manipulate,
///   use NLW2_MakeNLSolver_C() / NLW2_DestroyNLSolver_C(),
///   NLW2_LoadNLModel_C(),
///   etc, see the below API.
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
/// @note Not used for loading with NLW2_NLFeeder_C.
void NLW2_SetNLOptions_C(NLW2_NLSolver_C* ,
                         NLW2_NLOptionsBasic_C );

/// Get NLOptions used for loading NLModel.
NLW2_NLOptionsBasic_C NLW2_GetNLOptions_C(NLW2_NLSolver_C* );

/// Get error message.
/// Nonempty iff error occurred.
const char* NLW2_GetErrorMessage_C(NLW2_NLSolver_C* );

/// Load and solve NLW2_NLModel_C and return result.
///
/// @return NLW2_Solution_C object with computed obj_value_.
///   Valid as long as the NLW2_NLSolver_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @see NLW2_LoadNLModel_C(), NLW2_RunSolver_C(),
///   NLW2_ReadSolution_C() for details.
NLW2_NLSolution_C NLW2_SolveNLModel_C(NLW2_NLSolver_C* ,
                                    NLW2_NLModel_C* ,
                                    const char* solver,
                                    const char* solver_opts);

/// Load NLW2_NLFeeder_C, solve,
/// and read into SOLHandler_C.
///
/// @return Non-zero iff all ok.
///
/// @see NLW2_LoadNLFeed2_C(), NLW2_RunSolver_C(),
///   NLW2_Read2SOLHanlder2_C() for details.
int NLW2_SolveFeederHandler_C(NLW2_NLSolver_C* ,
                              NLW2_NLFeeder_C* ,
                              NLW2_SOLHandler_C* ,
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
int NLW2_LoadNLFeed2_C(NLW2_NLSolver_C* , NLW2_NLFeeder_C* nlf_c);

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
/// @return NLW2_Solution_C object.
///   Valid as long as the NLW2_NLSolver_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @note To compute objective value,
///   execute NLW2_ComputeObjValue_C()
///   if x_ available.
NLW2_NLSolution_C NLW2_ReadSolution_C(NLW2_NLSolver_C* );

/// Read solution to SOLHandler.
///
/// @param solh_c: solution handler.
///
/// @return Non-0 if all ok.
int NLW2_Read2SOLHandler_C(NLW2_NLSolver_C* ,
                            NLW2_SOLHandler_C* solh_c);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLSOLVER_C_H

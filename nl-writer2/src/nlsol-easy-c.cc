/**
 NL Solver "Easy", for special model classes

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

#include "api/c/nlsol-easy-c.h"
#include "api/c/nl-writer2-misc-c-impl.h"
#include "mp/nlsol-easy.h"

/// Cast & check
template <class Type>
Type* CastNZ(void* p) {
  auto result = (Type*)p;
  assert(result);
  return result;
}

#ifdef __cplusplus  // Implementing C API from C++
extern "C" {
#endif

/// Construct NLW2_NLModel_Easy_C
///
/// @param probname: can be NULL.
NLW2_NLModel_Easy_C NLW2_Make_NLModel_Easy_C(const char* probname) {
  NLW2_NLModel_Easy_C nlme;
  nlme.p_data_ = new mp::NLModel_Easy(probname);
  return nlme;
}

/// Destroy NLW2_NLModel_Easy_C
void NLW2_Destroy_NLModel_Easy_C(NLW2_NLModel_Easy_C* nlme) {
  delete CastNZ<mp::NLModel_Easy>(nlme->p_data_);
  nlme->p_data_ = nullptr;
}

/// Add variables (all at once.)
void NLME_SetCols_C(NLW2_NLModel_Easy_C* nlme, NLW2_ColData_C vd)
{ CastNZ<mp::NLModel_Easy>(nlme->p_data_)->SetCols(vd); }

/// Add variable names
void NLME_SetColNames_C(NLW2_NLModel_Easy_C* nlme, const char *const *nm)
{ CastNZ<mp::NLModel_Easy>(nlme->p_data_)->SetColNames(nm); }

/// Add linear constraints (all at once).
/// Only rowwise matrix supported.
void NLME_SetRows_C(NLW2_NLModel_Easy_C* nlme,
    int nr, const double* rlb, const double* rub,
    NLW2_SparseMatrix_C A)
{ CastNZ<mp::NLModel_Easy>(nlme->p_data_)->SetRows(nr, rlb, rub, A); }

/// Add constraint names
void NLME_SetRowNames_C(NLW2_NLModel_Easy_C* nlme, const char *const *nm)
{ CastNZ<mp::NLModel_Easy>(nlme->p_data_)->SetRowNames(nm); }

/// Add linear objective (only single objective supported.)
/// Sense: NLW2_ObjSenseM....
/// Coefficients: dense vector.
void NLME_SetLinearObjective_C(NLW2_NLModel_Easy_C* nlme,
                               int sense, double c0, const double* c) {
  CastNZ<mp::NLModel_Easy>(nlme->p_data_)
      ->SetLinearObjective(sense, c0, c);
}

/// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
/// Format: NLW2_HessianFormat...
void NLME_SetHessian_C(NLW2_NLModel_Easy_C* nlme,
                       int format, NLW2_SparseMatrix_C Q)
{ CastNZ<mp::NLModel_Easy>(nlme->p_data_)->SetHessian(format, Q); }

/// Set obj name
void NLME_SetObjName_C(NLW2_NLModel_Easy_C* nlme, const char* nm)
{ CastNZ<mp::NLModel_Easy>(nlme->p_data_)->SetObjName(nm); }

/// Compute objective value
double NLME_ComputeObjValue_C(NLW2_NLModel_Easy_C* nlme,
                            const double* x) {
  return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ComputeObjValue(x);
}


/// Get problem name
const char* NLME_ProbName_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ProbName(); }
/// Get variables
NLW2_ColData_C NLME_ColData_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ColData(); }
/// Get var names
const char *const *NLME_ColNames_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ColNames(); }
/// Get var name [i]
const char *NLME_ColName_C(NLW2_NLModel_Easy_C* nlme, int i)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ColName(i); }
/// Lin con matrix
NLW2_SparseMatrix_C NLME_GetA_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->GetA(); }
/// N cols
int NLME_NumCols_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->NumCols(); }
/// N rows
int NLME_NumRows_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->NumRows(); }
/// Row lb
const double *NLME_RowLowerBounds_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->RowLowerBounds(); }
/// Row ub
const double *NLME_RowUpperBounds_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->RowUpperBounds(); }
/// Row names
const char *const *NLME_RowNames_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->RowNames(); }
/// Row name [i]
const char *NLME_RowName_C(NLW2_NLModel_Easy_C* nlme, int i)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->RowName(i); }
/// Obj sense
int NLME_ObjSense_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ObjSense(); }
/// Obj offset
double NLME_ObjOffset_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ObjOffset(); }
/// Obj coefs
const double *NLME_ObjCoefficients_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ObjCoefficients(); }
/// Hessian format NLW2_HessianFormat...
int NLME_HessianFormat_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->HessianFormat(); }
/// Hessian matrix
NLW2_SparseMatrix_C NLME_Hessian_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->Hessian(); }
/// Obj name
const char* NLME_ObjName_C(NLW2_NLModel_Easy_C* nlme)
{ return CastNZ<mp::NLModel_Easy>(nlme->p_data_)->ObjName(); }


/// Construct.
NLW2_NLSOL_Easy_C NLW2_Make_NLSOL_Easy_C(NLW2_NLUtils_C* utl) {
  NLW2_NLSOL_Easy_C result {};         // p_nlutl_ = p_sol_ = 0
  result.p_nlse_ = new mp::NLSOL_Easy;
  if (utl)
    result.p_nlutl_ = new mp::NLUtils_C_Impl(utl);
  return result;
}
/// Destruct
void NLW2_Destroy_NLSOL_Easy_C(NLW2_NLSOL_Easy_C* nlse) {
  delete CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_);
  nlse->p_nlse_ = nullptr;
  if (nlse->p_nlutl_) {
    delete CastNZ<mp::NLUtils_C_Impl>(nlse->p_nlutl_);
    nlse->p_nlutl_ = nullptr;
  }
}


/// Set file stub [OPTIONAL].
///
/// Used for filename base of .nl, .col, row, etc. input files,
/// as well as .sol output files.
///
/// If not provided, a temporary filename is used;
/// then, .nl is deleted upon object destruction.
void NLSE_SetFileStub_C(NLW2_NLSOL_Easy_C* nlse, const char* stub)
{ CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_)->SetFileStub(stub); }

/// Retrieve file stub.
const char* NLSE_GetFileStub_C(NLW2_NLSOL_Easy_C* nlse)
{ return CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_)->GetFileStub().c_str(); }

/// Set NL options [OPTIONAL].
///
/// If not provided, default is used.
void NLSE_SetNLOptions_C(NLW2_NLSOL_Easy_C* nlse,
                         NLW2_NLOptionsBasic_C nlo)
{ CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_)->SetNLOptions(nlo); }

/// Get NLOptions
NLW2_NLOptionsBasic_C NLSE_GetNLOptions_C(NLW2_NLSOL_Easy_C* nlse)
{ return CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_)->GetNLOptions(); }

/// Get error message.
/// Nonempty iff error occurred.
const char* NLSE_GetErrorMessage_C(NLW2_NLSOL_Easy_C* nlse)
{ return CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_)->GetErrorMessage(); }

/// Load and solve model and return result.
///
/// @return Solution object.
///   Valid as long as the NLW2_NLSOL_Easy_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @see NLSE_LoadModel_C(), NLSE_RunSolver_C(),
///   NLSE_ReadSolution_C() for details.
NLSE_Solution_C NLSE_Solve_C(NLW2_NLSOL_Easy_C* nlse,
                             NLW2_NLModel_Easy_C* nlme,
                             const char* solver,
                             const char* solver_opts) {
  NLSE_Solution_C result;

  auto sol = CastNZ<mp::NLSOL_Easy>(nlse->p_nlse_)->Solve(
        *CastNZ<mp::NLModel_Easy>(nlme->p_data_), solver, solver_opts);

  return result;
}

/// Write NL and any accompanying files.
/// NL file name base and some options
/// can be provided, if non-defaults desired,
/// via SetFileStub() and SetNLOptions().
///
/// @return true if all ok, otherwise see
///   GetErrorMessage().
int NLSE_LoadModel_C(NLW2_NLSOL_Easy_C* ,
                     NLW2_NLModel_Easy_C* );

/// RunSolver: run the given solver after loading the model.
///
/// @param solver: solver executable, such as "gurobi".
/// @param solver_opts: string of solver options,
///   such as "outlev=1 writeprob=model.lp".
///
/// @return true if all ok.
int NLSE_RunSolver_C(NLW2_NLSOL_Easy_C* ,
                     const char* solver,
                     const char* solver_opts);

/// Read solution.
///
/// @return Solution object.
///   Valid as long as the NLW2_NLSOL_Easy_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @note To compute objective value,
///   execute NLME_ComputeObjValue_C()
///   if x_ available.
NLSE_Solution_C NLSE_ReadSolution_C(NLW2_NLSOL_Easy_C* );

#ifdef __cplusplus
}  // extern "C"
#endif

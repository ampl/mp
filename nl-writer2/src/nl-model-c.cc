/*
 NL Model C API, for special model classes

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

#include "api/c/nl-model-c.h"
#include "mp/nl-model.h"

namespace {

/// Cast & check
template <class Type>
Type* CastNZ(void* p) {
  auto result = (Type*)p;
  assert(result);
  return result;
}

}

#ifdef __cplusplus  // Implementing C API from C++
extern "C" {
#endif

/// Construct NLW2_NLModel_C
///
/// @param probname: can be NULL.
NLW2_NLModel_C NLW2_MakeNLModel_C(const char* probname) {
  NLW2_NLModel_C nlme;
  nlme.p_data_ = new mp::NLModel(probname);
  return nlme;
}

/// Destroy NLW2_NLModel_C
void NLW2_DestroyNLModel_C(NLW2_NLModel_C* nlme) {
  delete CastNZ<mp::NLModel>(nlme->p_data_);
  nlme->p_data_ = nullptr;
}

/// Add variables (all at once.)
void NLW2_SetCols_C(NLW2_NLModel_C* nlme,
                    int num_col,
                    const double *lower,
                    const double *upper,
                    const int *type)
{ CastNZ<mp::NLModel>(nlme->p_data_)
      ->SetCols({num_col, lower, upper, type}); }

/// Add variable names
void NLW2_SetColNames_C(NLW2_NLModel_C* nlme, const char *const *nm)
{ CastNZ<mp::NLModel>(nlme->p_data_)->SetColNames(nm); }

/// Add linear constraints (all at once).
/// Only rowwise matrix supported.
void NLW2_SetRows_C(NLW2_NLModel_C* nlme,
                    int nr, const double* rlb, const double* rub,
                    NLW2_MatrixFormat format,
                    size_t num_nz,
                    const size_t *start,
                    const int *index,
                    const double *value)
{ CastNZ<mp::NLModel>(nlme->p_data_)
      ->SetRows(nr, rlb, rub,
                {nr, format, num_nz, start, index, value}); }

/// Add constraint names
void NLW2_SetRowNames_C(NLW2_NLModel_C* nlme, const char *const *nm)
{ CastNZ<mp::NLModel>(nlme->p_data_)->SetRowNames(nm); }

/// Add linear objective (only single objective supported.)
/// Sense: NLW2_ObjSenseM....
/// Coefficients: dense vector.
void NLW2_SetLinearObjective_C(NLW2_NLModel_C* nlme,
                               NLW2_ObjSense sense,
                               double c0, const double* c) {
  CastNZ<mp::NLModel>(nlme->p_data_)
      ->SetLinearObjective(sense, c0, c);
}

/// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
/// Format: NLW2_HessianFormat...
void NLW2_SetHessian_C(NLW2_NLModel_C* nlme,
                       NLW2_HessianFormat format,
                       int dim,
                       size_t num_nz_,
                       const size_t *start_,
                       const int *index_,
                       const double *value_)
{ CastNZ<mp::NLModel>(nlme->p_data_)
      ->SetHessian(format, {
                     dim, NLW2_MatrixFormatIrrelevant,
                     num_nz_, start_, index_, value_
                   }); }

/// Set obj name
void NLW2_SetObjName_C(NLW2_NLModel_C* nlme, const char* nm)
{ CastNZ<mp::NLModel>(nlme->p_data_)->SetObjName(nm); }


void NLW2_SetWarmstart_C(NLW2_NLModel_C* nlme,
                         NLW2_SparseVector_C ini_x)
{ CastNZ<mp::NLModel>(nlme->p_data_)->SetWarmstart(ini_x); }

void NLW2_SetDualWarmstart_C(NLW2_NLModel_C* nlme,
                             NLW2_SparseVector_C ini_y)
{ CastNZ<mp::NLModel>(nlme->p_data_)->SetWarmstart(ini_y); }

int NLW2_AddSuffix_C(NLW2_NLModel_C* nlme,
                     NLW2_NLSuffix_C suf_c) {
  mp::NLSuffix suf{
    suf_c.name_, suf_c.table_, suf_c.kind_,
    {suf_c.values_, suf_c.values_+suf_c.numval_}
  };
  return CastNZ<mp::NLModel>(nlme->p_data_)->AddSuffix(suf);
}


/// Compute objective value
double NLW2_ComputeObjValue_C(NLW2_NLModel_C* nlme,
                            const double* x) {
  return CastNZ<mp::NLModel>(nlme->p_data_)->ComputeObjValue(x);
}


/// Get problem name
const char* NLW2_ProbName_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ProbName(); }
/// Get variables
NLW2_ColData_C NLW2_Columns_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ColData(); }
/// Get var names
const char *const *NLW2_ColNames_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ColNames(); }
/// Get var name [i]
const char *NLW2_ColName_C(NLW2_NLModel_C* nlme, int i)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ColName(i); }
/// Lin con matrix
NLW2_SparseMatrix_C NLW2_GetA_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->GetA(); }
/// N cols
int NLW2_NumCols_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->NumCols(); }
/// N rows
int NLW2_NumRows_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->NumRows(); }
/// Row lb
const double *NLW2_RowLowerBounds_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->RowLowerBounds(); }
/// Row ub
const double *NLW2_RowUpperBounds_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->RowUpperBounds(); }
/// Row names
const char *const *NLW2_RowNames_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->RowNames(); }
/// Row name [i]
const char *NLW2_RowName_C(NLW2_NLModel_C* nlme, int i)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->RowName(i); }
/// Obj sense
int NLW2_ObjSense_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ObjSense(); }
/// Obj offset
double NLW2_ObjOffset_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ObjOffset(); }
/// Obj coefs
const double *NLW2_ObjCoefficients_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ObjCoefficients(); }
/// Hessian format NLW2_HessianFormat...
int NLW2_HessianFormat_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->HessianFormat(); }
/// Hessian matrix
NLW2_SparseMatrix_C NLW2_Hessian_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->Hessian(); }
/// Obj name
const char* NLW2_ObjName_C(NLW2_NLModel_C* nlme)
{ return CastNZ<mp::NLModel>(nlme->p_data_)->ObjName(); }


#ifdef __cplusplus
}  // extern "C"
#endif

/**
 * C API: extern "C" wrapper for mp::NLModel class
 *
 */

#ifndef NLMODEL_C_H
#define NLMODEL_C_H

#include "mp/nl-solver-basics-c.h"

#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif


/// extern "C" wrapper of mp::NLModel.
///
/// Intermediate representation for special model types:
/// (MI)LP, (MI)QP.
/// For full modeling capabilities
/// use NLW2_NLFeeder_C, NLW2_SOLHandler_C.
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
///
/// @see C API tests/examples.
typedef struct NLW2_NLModel_C {
  void *p_data_;
} NLW2_NLModel_C;

/// Construct NLW2_NLModel_C
///
/// @param probname: can be NULL.
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
NLW2_NLModel_C NLW2_MakeNLModel_C(const char* probname);

/// Destroy NLW2_NLModel_C
void NLW2_DestroyNLModel_C(NLW2_NLModel_C* );

/// Add variables (all at once.)
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetCols_C(NLW2_NLModel_C* ,
                    int num_col,
                    const double *lower,
                    const double *upper,
                    const int *type);

/// Add variable names
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetColNames_C(NLW2_NLModel_C* , const char *const *nm);

/// Add linear constraints (all at once).
/// Only rowwise matrix supported.
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetRows_C(NLW2_NLModel_C* ,
                    int nr, const double* rlb, const double* rub,
                    NLW2_MatrixFormat format_,
                    size_t num_nz_,
                    const size_t *start_,
                    const int *index_,
                    const double *value_);

/// Add constraint names
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetRowNames_C(NLW2_NLModel_C* , const char *const *nm);

/// Add linear objective (only single objective supported.)
/// Coefficients: dense vector.
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetLinearObjective_C(NLW2_NLModel_C* ,
                               NLW2_ObjSense sense,
                               double c0, const double* c);

/// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetHessian_C(NLW2_NLModel_C* ,
                       NLW2_HessianFormat format,
                       int dim,
                       size_t num_nz_,
                       const size_t *start_,
                       const int *index_,
                       const double *value_);

/// Set obj name
///
/// @note All model data pointers should stay valid until
/// loading the model into NLW2_NLSolver_C.
void NLW2_SetObjName_C(NLW2_NLModel_C* , const char* nm);

/// Compute objective value
double NLW2_ComputeObjValue_C(NLW2_NLModel_C* ,
                            const double* x);


/// Get problem name
const char* NLW2_ProbName_C(NLW2_NLModel_C* );
/// Get variables
NLW2_ColData_C NLW2_Columns_C(NLW2_NLModel_C* );
/// Get var names
const char *const *NLW2_ColNames_C(NLW2_NLModel_C* );
/// Get var name [i]
const char *NLW2_ColName_C(NLW2_NLModel_C* , int i);
/// Lin con matrix
NLW2_SparseMatrix_C NLW2_GetA_C(NLW2_NLModel_C* );
/// N cols
int NLW2_NumCols_C(NLW2_NLModel_C* );
/// N rows
int NLW2_NumRows_C(NLW2_NLModel_C* );
/// Row lb
const double *NLW2_RowLowerBounds_C(NLW2_NLModel_C* );
/// Row ub
const double *NLW2_RowUpperBounds_C(NLW2_NLModel_C* );
/// Row names
const char *const *NLW2_RowNames_C(NLW2_NLModel_C* );
/// Row name [i]
const char *NLW2_RowName_C(NLW2_NLModel_C* , int i);
/// Obj sense
int NLW2_ObjSense_C(NLW2_NLModel_C* );
/// Obj offset
double NLW2_ObjOffset_C(NLW2_NLModel_C* );
/// Obj coefs
const double *NLW2_ObjCoefficients_C(NLW2_NLModel_C* );
/// Hessian format NLW2_HessianFormat...
int NLW2_HessianFormat_C(NLW2_NLModel_C* );
/// Hessian matrix
NLW2_SparseMatrix_C NLW2_Hessian_C(NLW2_NLModel_C* );
/// Obj name
const char* NLW2_ObjName_C(NLW2_NLModel_C* );


/// NL suffix type
typedef struct NLW2_Suffix_C {
  /// Name
  const char* name_;
  /// Suffix table
  const char* table_;
  /// Kind
  int kind_;
  /// Values. Always double precision.
  const double* values_;
} NLW2_Suffix_C;


/// NL solution
typedef struct NLW2_Solution_C {
  /// Solve result.
  /// If >-2, solver interaction successful.
  /// Then:
  /// -1      'unknown' unexpected termination
  /// 0- 99   'solved' optimal solution found
  /// 100-199 'solved?' optimal solution indicated, but error likely
  /// 200-299 'infeasible' constraints cannot be satisfied
  /// 300-399 'unbounded' objective can be improved without limit
  /// 400-499 'limit' stopped by a limit that you set (such as on iterations)
  /// 500-999 'failure' stopped by an error condition in the solver
  ///
  /// Individual solvers may have more specific values,
  /// see https://ampl.com/products/solvers/solvers-we-sell/.
  int solve_result_;
  /// Number of solve_message's initial characters
  /// already printed on the screen
  int nbs_;
  /// Solve message
  const char* solve_message_;
  /// Objective value.
  /// Only returned by NLW2_Solve_C().
  /// Otherwise, after NLW2_ReadSolution...,
  /// should be manually computed, e.g.,
  /// by NLW2_ComputeObjValue_C().
  double obj_val_;
  /// N primal values
  int n_primal_values_;
  /// Primals
  const double* x_;
  /// N dual values
  int n_dual_values_;
  /// Duals
  const double* y_;
  /// Num suffixes
  int nsuf_;
  /// Suffixes
  const NLW2_Suffix_C* suffixes_;
} NLW2_Solution_C;


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLMODEL_C_H

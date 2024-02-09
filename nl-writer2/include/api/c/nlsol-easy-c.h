/**
 * C API: extern "C" wrapper for the NLSOL_Easy class
 *
 */

#ifndef NLSOLEASYC_H
#define NLSOLEASYC_H

#include "mp/basic-defs-c.h"
#include "api/c/nl-writer2-misc-c.h"

#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif


/// extern "C" wrapper of mp::NLModel_Easy.
///
/// Intermediate representation for special model types:
/// (MI)LP, (MI)QP.
///
/// All model data pointers should stay valid until
/// loading the model into NLW2_NLSOL_Easy_C.
///
/// @see C API tests/examples.
typedef struct NLW2_NLModel_Easy_C {
  void *p_data_;
} NLW2_NLModel_Easy_C;

/// Construct NLW2_NLModel_Easy_C
///
/// @param probname: can be NULL.
NLW2_NLModel_Easy_C NLW2_MakeNLModel_Easy_C(const char* probname);

/// Destroy NLW2_NLModel_Easy_C
void NLW2_DestroyNLModel_Easy_C(NLW2_NLModel_Easy_C* );

/// Add variables (all at once.)
void NLME_SetCols_C(NLW2_NLModel_Easy_C* ,
                    int num_col,
                    const double *lower,
                    const double *upper,
                    const int *type);

/// Add variable names
void NLME_SetColNames_C(NLW2_NLModel_Easy_C* , const char *const *nm);

/// Add linear constraints (all at once).
/// Only rowwise matrix supported.
void NLME_SetRows_C(NLW2_NLModel_Easy_C* ,
                    int nr, const double* rlb, const double* rub,
                    int format_,
                    size_t num_nz_,
                    const size_t *start_,
                    const int *index_,
                    const double *value_);

/// Add constraint names
void NLME_SetRowNames_C(NLW2_NLModel_Easy_C* , const char *const *nm);

/// Add linear objective (only single objective supported.)
/// Sense: NLW2_ObjSenseM....
/// Coefficients: dense vector.
void NLME_SetLinearObjective_C(NLW2_NLModel_Easy_C* ,
                               int sense, double c0, const double* c);

/// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
/// Format: NLW2_HessianFormat...
void NLME_SetHessian_C(NLW2_NLModel_Easy_C* ,
                       int format,
                       int dim,
                       size_t num_nz_,
                       const size_t *start_,
                       const int *index_,
                       const double *value_);

/// Set obj name
void NLME_SetObjName_C(NLW2_NLModel_Easy_C* , const char* nm);

/// Compute objective value
double NLME_ComputeObjValue_C(NLW2_NLModel_Easy_C* ,
                            const double* x);


/// Get problem name
const char* NLME_ProbName_C(NLW2_NLModel_Easy_C* );
/// Get variables
NLW2_ColData_C NLME_ColData_C(NLW2_NLModel_Easy_C* );
/// Get var names
const char *const *NLME_ColNames_C(NLW2_NLModel_Easy_C* );
/// Get var name [i]
const char *NLME_ColName_C(NLW2_NLModel_Easy_C* , int i);
/// Lin con matrix
NLW2_SparseMatrix_C NLME_GetA_C(NLW2_NLModel_Easy_C* );
/// N cols
int NLME_NumCols_C(NLW2_NLModel_Easy_C* );
/// N rows
int NLME_NumRows_C(NLW2_NLModel_Easy_C* );
/// Row lb
const double *NLME_RowLowerBounds_C(NLW2_NLModel_Easy_C* );
/// Row ub
const double *NLME_RowUpperBounds_C(NLW2_NLModel_Easy_C* );
/// Row names
const char *const *NLME_RowNames_C(NLW2_NLModel_Easy_C* );
/// Row name [i]
const char *NLME_RowName_C(NLW2_NLModel_Easy_C* , int i);
/// Obj sense
int NLME_ObjSense_C(NLW2_NLModel_Easy_C* );
/// Obj offset
double NLME_ObjOffset_C(NLW2_NLModel_Easy_C* );
/// Obj coefs
const double *NLME_ObjCoefficients_C(NLW2_NLModel_Easy_C* );
/// Hessian format NLW2_HessianFormat...
int NLME_HessianFormat_C(NLW2_NLModel_Easy_C* );
/// Hessian matrix
NLW2_SparseMatrix_C NLME_Hessian_C(NLW2_NLModel_Easy_C* );
/// Obj name
const char* NLME_ObjName_C(NLW2_NLModel_Easy_C* );


/// extern "C" wrapper of mp::NLSOL_Easy.
///
/// Provides "easy" API to use AMPL solvers
/// for (MI)QP models.
typedef struct NLW2_NLSOL_Easy_C {
  void *p_nlse_;
  void *p_nlutl_;
  void *p_sol_;
} NLW2_NLSOL_Easy_C;

/// Construct.
/// @param nlu can be NULL.
NLW2_NLSOL_Easy_C NLW2_MakeNLSOL_Easy_C(NLW2_NLUtils_C* nlu);
/// Destruct
void NLW2_DestroyNLSOL_Easy_C(NLW2_NLSOL_Easy_C* );


/// Set file stub [OPTIONAL].
///
/// Used for filename base of .nl, .col, row, etc. input files,
/// as well as .sol output files.
///
/// If not provided, a temporary filename is used;
/// then, .nl is deleted upon object destruction.
void NLSE_SetFileStub_C(NLW2_NLSOL_Easy_C* , const char* );

/// Retrieve file stub.
const char* NLSE_GetFileStub_C(NLW2_NLSOL_Easy_C* );

/// Set NL options [OPTIONAL].
///
/// If not provided, default is used.
void NLSE_SetNLOptions_C(NLW2_NLSOL_Easy_C* ,
                         NLW2_NLOptionsBasic_C );

/// Get NLOptions
NLW2_NLOptionsBasic_C NLSE_GetNLOptions_C(NLW2_NLSOL_Easy_C* );

/// Get error message.
/// Nonempty iff error occurred.
const char* NLSE_GetErrorMessage_C(NLW2_NLSOL_Easy_C* );

/// Suffix type
typedef struct NLSE_Suffix_C {
  /// Name
  const char* name_;
  /// Suffix table
  const char* table_;
  /// Kind
  int kind_;
  /// Values. Always double precision.
  const double* values_;
} NLSE_Suffix_C;


/// Solution
typedef struct NLSE_Solution_C {
  /// Solve result.
  /// >=-1 if contains data.
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
  const NLSE_Suffix_C* suffixes_;
} NLSE_Solution_C;

/// Load and solve model and return result.
///
/// @return Solution object.
///   Valid as long as the NLW2_NLSOL_Easy_C object lives,
///   and until the next Solve() or ReadSolution().
///
/// @see NLSE_LoadModel_C(), NLSE_RunSolver_C(),
///   NLSE_ReadSolution_C() for details.
NLSE_Solution_C NLSE_Solve_C(NLW2_NLSOL_Easy_C* ,
                             NLW2_NLModel_Easy_C* ,
                             const char* solver,
                             const char* solver_opts);

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

#endif // NLSOLEASYC_H

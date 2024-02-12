#ifndef NLW2_BASIC_DEFS_C_H
#define NLW2_BASIC_DEFS_C_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

const int NLW2_ObjSenseMinimize = 0;
const int NLW2_ObjSenseMaximize = 1;

const int NLW2_VarTypeContinuous = 0;
const int NLW2_VarTypeInteger = 1;
//const int NLW2_VarTypeSemiContinuous = 2;
//const int NLW2_VarTypeSemiInteger = 3;
//const int NLW2_VarTypeImplicitInteger = 4;

//const int NLW2_MatrixFormatColwise = 1;
const int NLW2_MatrixFormatRowwise = 2;

const int NLW2_HessianFormatTriangular = 1;
const int NLW2_HessianFormatSquare = 2;

/// Variables' data by pointers
typedef struct NLW2_ColData_C {
  /// Num vars
  int num_col_;
  /// lower bounds
  const double *lower_;
  /// upper bounds
  const double *upper_;
  /// type: NLW2_VarType...
  /// Set to NULL if all continuous.
  const int *type_;
} NLW2_ColData_C;

/// Sparse matrix.
typedef struct NLW2_SparseMatrix_C {
  /// Size of the start_ array:
  /// N cols (for colwise) / N rows (for rowwise),
  /// depending on format_.
  int num_colrow_;
  /// Format (NLW2_MatrixFormat...).
  /// Only rowwise supported.
  int format_;
  /// Nonzeros
  size_t num_nz_;
  /// Row / col starts
  const size_t *start_;
  /// Entry index
  const int *index_;
  /// Entry value
  const double *value_;
} NLW2_SparseMatrix_C;


/// Basic NL options for NLModel_Easy.
/// Prefer to create by NLW2_MakeNLOptionsBasic_Default().
typedef struct NLW2_NLOptionsBasic_C {
  /// NL text mode?
  int n_text_mode_;
  /// NL comments in text mode?
  int want_nl_comments_;
  /// Flags (1== want output suffixes)
  int flags_;
} NLW2_NLOptionsBasic_C;

/// Use this to create default NL options for NLModel_Easy.
NLW2_NLOptionsBasic_C NLW2_MakeNLOptionsBasic_C_Default();

/// NLW2_WriteNLResultCode enum.
enum NLW2_WriteNLResultCode {
  NLW2_WriteNL_Unset = 0,
  NLW2_WriteNL_OK = 1,
  NLW2_WriteNL_CantOpen,
  NLW2_WriteNL_Failed
};

/// SOL read result code
enum NLW2_SOLReadResultCode {
  NLW2_SOLRead_Result_Not_Set = -1,
  NLW2_SOLRead_OK = 0,
  NLW2_SOLRead_Fail_Open,
  NLW2_SOLRead_Early_EOF,
  NLW2_SOLRead_Bad_Format,
  NLW2_SOLRead_Bad_Line,
  NLW2_SOLRead_Bad_Options,
  NLW2_SOLRead_Vector_Not_Finished,
  NLW2_SOLRead_Bad_Suffix
};


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLW2_BASIC_DEFS_C_H

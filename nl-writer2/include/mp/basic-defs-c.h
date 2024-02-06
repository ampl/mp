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

const int NLW2_MatrixFormatColwise = 1;
// const int NLW2_MatrixFormatRowwise = 2;

const int NLW2_HessianFormatTriangular = 1;
const int NLW2_HessianFormatSquare = 2;

/// Variables' data by pointers
struct NLW2_ColData {
  /// Num vars
  int num_col_;
  /// lower bounds
  const double *lower_;
  /// upper bounds
  const double *upper_;
  /// type: NLW2_VarType...
  /// Set to NULL if all continuous.
  const int *type_;
};

/// Sparse matrix.
struct NLW2_SparseMatrix {
  /// N rows (for rowwise) / N cols (colwise),
  /// depending on format_.
  int num_row_;
  /// Format (NLW2_MatrixFormat...).
  /// Only rowwise supported.
  int format_;
  /// Nonzeros
  size_t num_nz_;
  /// Row / col starts
  const int *start_;
  /// Entry index
  const int *index_;
  /// Entry value
  const double *value_;
};


/// Basic NL options for NLModel_Easy.
/// Prefer to create by NLW2_MakeNLOptionsBasic_Default().
struct NLW2_NLOptionsBasic {
  /// NL text mode?
  int n_text_mode_;
  /// NL comments in text mode?
  int want_nl_comments_;
  /// Flags (1== want output suffixes)
  int flags_;
};

/// Use this to create default NL options for NLModel_Easy.
NLW2_NLOptionsBasic NLW2_MakeNLOptionsBasic_Default();

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLW2_BASIC_DEFS_C_H

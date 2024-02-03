/**
 NL Solver "Easy", for special model classes

 Copyright (C) 2023 AMPL Optimization Inc.

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
#ifndef NLSOLEASY_H
#define NLSOLEASY_H

#include "mp/basic-defs-c.h"

namespace mp {

/// Class NLModelEasy.
///
/// Intermediate representation for special model types:
/// LP, QP.
///
/// All pointers should stay valid until
/// loading the model into NLSOL_Easy.
class NLModelEasy {
public:
  /**
   *  Build an (MI)QP model in a single function call.
   *  For (MI)LP, pass NULL for Q and/or integrality.
   *
   * @param num_col     The number of columns.
   * @param num_row     The number of rows.
   * @param num_nz      The number of elements in the constraint matrix.
   * @param q_num_nz    The number of elements in the Hessian matrix.
   * @param a_format    The format of the constraint matrix to use in the form
   *                    of a `NLW2_MatrixFormat` constant.
   * @param q_format    The format of the Hessian matrix to use in the form of a
   *                    `NLW2_HessianFormat` constant.
   * @param sense       The optimization sense in the form of a `NLW2_ObjSense`
   *                    constant.
   * @param offset      The constant term in the objective function.
   * @param col_cost    An array of length [num_col] with the objective
   *                    coefficients.
   * @param col_lower   An array of length [num_col] with the lower column
   *                    bounds.
   * @param col_upper   An array of length [num_col] with the upper column
   *                    bounds.
   * @param row_lower   An array of length [num_row] with the upper row bounds.
   * @param row_upper   An array of length [num_row] with the upper row bounds.
   * @param a_start     The constraint matrix is provided in compressed
   *                    sparse column form (if `a_format` is
   *                    `NLW2_MatrixFormatColwise`, otherwise compressed sparse
   *                    row form). The sparse matrix consists of three arrays,
   *                    `a_start`, `a_index`, and `a_value`. `a_start` is an
   * array of length [num_col] containing the starting index of each column in
   * `a_index`. If `a_format` is `NLW2_MatrixFormatRowwise` the array is of
   * length [num_row] corresponding to each row.
   * @param a_index     An array of length [num_nz] with indices of matrix
   *                    entries.
   * @param a_value     An array of length [num_nz] with values of matrix
   *                    entries.
   * @param q_start     The Hessian matrix is provided in the same format as the
   *                    constraint matrix, using `q_start`, `q_index`, and
   *                    `q_value` in the place of `a_start`, `a_index`, and
   *                    `a_value`. If the model is linear, pass NULL.
   * @param q_index     An array of length [q_num_nz] with indices of matrix
   *                    entries. If the model is linear, pass NULL.
   * @param q_value     An array of length [q_num_nz] with values of matrix
   *                    entries. If the model is linear, pass NULL.
   * @param integrality An array of length [num_col] containing a
   *                    `NLW2_VarType` constant for each column.
   *
   * @return true iff ok.
   */
  bool BuildMIQP(void *highs, const int num_col, const int num_row,
                 const int num_nz, const int q_num_nz, const int a_format,
                 const int q_format, const int sense, const double offset,
                 const double *col_cost, const double *col_lower,
                 const double *col_upper, const double *row_lower,
                 const double *row_upper, const int *a_start,
                 const int *a_index, const double *a_value, const int *q_start,
                 const int *q_index, const double *q_value,
                 const int *integrality);
};

}  // namespace mp

#endif // NLSOLEASY_H

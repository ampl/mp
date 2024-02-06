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

#include <string>

#include "mp/basic-defs-c.h"
#include "mp/nl-utils2.h"

namespace mp {

/// Class NLModel_Easy.
///
/// Intermediate representation for special model types:
/// LP, QP.
///
/// All pointers should stay valid until
/// loading the model into NLSOL_Easy.
class NLModel_Easy {
public:
  /// Construct
  NLModel_Easy(const char* probname = nullptr)
    : prob_name_(probname) { }

  /// Add variables (all at once.)
  void SetCols(NLW2_ColData vd) { vars_ = vd; }

  /// Add variable names
  void SetColNames(const char *const *nm) { var_names_=nm; }

  /// Add linear constraints (all at once).
  /// Only rowwise matrix supported.
  void SetRows(
      int nr, const double* rlb, const double* rub,
      NLW2_SparseMatrix A)
  { num_row_=nr; row_lb_=rlb; row_ub_=rub; A_=A; }

  /// Add constraint names
  void SetRowNames(const char *const *nm) { row_names_=nm; }

  /// Add linear objective (only single objective supported.)
  /// Sense: NLW2_ObjSenseM....
  /// Coefficients: dense vector.
  void SetLinearObjective(int sense, double c0, const double* c);

  /// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
  /// Format: NLW2_HessianFormat...
  void SetHessian(int format, NLW2_SparseMatrix Q)
  { Q_format_ = format; Q_ = Q; }

  /// Write to NL file.
  /// Recommended usage via class NLSOL_Easy.
  /// @return empty string iff ok.
  std::string WriteNL(const std::string& fln_base,
                      NLW2_NLOptionsBasic opts,
                      NLUtils &ut);


  /// Get problem name
  const char* ProbName() const { return prob_name_; }
  /// Get variables
  NLW2_ColData ColData() const { return vars_; }
  /// Get var names
  const char *const *ColNames() const { return var_names_; }
  /// Lin con matrix
  NLW2_SparseMatrix GetA() const { return A_; }
  /// N cols
  int NumCols() const { return vars_.num_col_; }
  /// N rows
  int NumRows() const { return num_row_; }
  /// Row lb
  const double *RowLowerBounds() const { return row_lb_; }
  /// Row ub
  const double *RowUpperBounds() const { return row_ub_; }
  /// Row names
  const char *const *RowNames() const { return row_names_; }
  /// Obj sense
  int ObjSense() const { return obj_sense_; }
  /// Obj offset
  double ObjOffset() const { return obj_c0_; }
  /// Obj coefs
  const double *ObjCoefficients() const { return obj_c_; }
  /// Hessian format NLW2_HessianFormat...
  int HessianFormat() const { return Q_format_; }
  /// Hessian matrix
  NLW2_SparseMatrix Hessian() const { return Q_; }
  /// Obj name
  const char* ObjName() const { return obj_name_; }

private:
  const char* prob_name_ {"NLSOL_Easy_model"};
  NLW2_ColData vars_ {};
  const char *const *var_names_ {};
  NLW2_SparseMatrix A_ {};
  int num_row_ {};
  const double *row_lb_ {};
  const double *row_ub_ {};
  const char *const *row_names_ {};
  int obj_sense_ {};
  double obj_c0_ {};
  const double *obj_c_ {};
  int Q_format_ {};
  NLW2_SparseMatrix Q_ {};
  const char* obj_name_ {"obj[1]"};
};


/// Class NLSOL_Easy.
///
/// A wrapper for mp::NLSOL to use AMPL solvers
/// for (MI)QP models.
class NLSOL_Easy {
public:
};

}  // namespace mp

#endif // NLSOLEASY_H

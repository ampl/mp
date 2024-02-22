/*
 NL model, for special model classes

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

#ifndef NLMODEL_H
#define NLMODEL_H

#include <string>
#include <vector>
#include <set>
#include <memory>
#include <cassert>

#include "mp/nl-solver-basics-c.h"

namespace mp {


/// NL suffix type
struct NLSuffix {
  /// Construct
  NLSuffix(std::string name, int kind, std::vector<double> v={})
    : name_(std::move(name)), table_({}), kind_(kind),
      values_(std::move(v)) { }
  /// Construct
  NLSuffix(std::string name, std::string table, int kind,
           std::vector<double> v={})
    : name_(std::move(name)), table_(std::move(table)),
      kind_(kind), values_(std::move(v)) { }

  /// Name
  std::string name_{};
  /// Suffix table
  std::string table_{};
  /// Kind.
  ///
  ///   VAR     =    0,  /**< Applies to variables. */
  ///   CON     =    1,  /**< Applies to constraints. */
  ///   OBJ     =    2,  /**< Applies to objectives. */
  ///   PROBLEM =    3   /**< Applies to problems. */
  ///
  /// If the suffix should be delivered as real-valued,
  /// the kind_ should be bitwise-OR'ed with 0x4.
  int kind_{-1};
  /// Values. Always double precision. Dense vector.
  std::vector<double> values_ {};

  /// operator<
  bool operator<(const NLSuffix& s) const {
    return std::make_pair(name_, kind_&3)
        < std::make_pair(s.name_, s.kind_&3);
  }
};

/// NL suffix set.
class NLSuffixSet : private std::set<NLSuffix> {
protected:
  using Base = std::set<NLSuffix>;
public:
  /// Add suffix.
  /// @return true iff new suffix (not existed before.)
  bool Add(NLSuffix suf)
  { return this->insert(suf).second; }

  /// Find suffix.
  /// @param k: kind (only the first 2 bits are used.)
  /// @return NLSuffix*, nullptr iff not found.
  const NLSuffix* Find(const std::string& nm, int k) const {
    NLSuffix tmp {nm, {}, k&3};
    auto it = this->find(tmp);
    return (this->end()!=it) ? &*it : nullptr;
  }

  /// Expose size, empty
  using Base::size;
  using Base::empty;
  using Base::clear;

  /// Expose begin, end
  using Base::begin;
  using Base::end;
  using Base::rbegin;
  using Base::rend;
  using Base::cbegin;
  using Base::cend;
  using Base::crbegin;
  using Base::crend;
};

/// Declare
class NLUtils;


/// Class NLModel.
///
/// Intermediate representation for special model types:
/// (MI)LP, (MI)QP.
/// For fully nonlinear models with expression trees,
/// use NLSolver with `NLFeeder`/`SOLHandler`.
///
/// @note All pointers should stay valid until
/// loading the model into NLSolver.
class NLModel {
public:
  /// Construct.
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  NLModel(const char* probname = nullptr)
    : prob_name_(probname ? probname : "NLModelInstance") { }

  /// Add variables (all at once.)
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetCols(NLW2_ColData_C vd) { vars_ = vd; }

  /// Add variable names
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetColNames(const char *const *nm) { var_names_=nm; }

  /// Add linear constraints (all at once).
  /// Only rowwise matrix supported.
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetRows(
      int nr, const double* rlb, const double* rub,
      NLW2_SparseMatrix_C A)
  { num_row_=nr; row_lb_=rlb; row_ub_=rub; A_=A; }

  /// Add constraint names
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetRowNames(const char *const *nm) { row_names_=nm; }

  /// Add linear objective (only single objective supported.)
  /// Coefficients: dense vector.
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetLinearObjective(NLW2_ObjSense sense, double c0,
                          const double* c = nullptr)
  { obj_sense_=sense; obj_c0_=c0; obj_c_=c; }

  /// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetHessian(NLW2_HessianFormat format, NLW2_SparseMatrix_C Q)
  { Q_format_ = format; Q_ = Q; }

  /// Set obj name
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetObjName(const char* nm) { obj_name_=(nm ? nm : ""); }

  /// Set initial solution.
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetWarmstart(NLW2_SparseVector_C ini_x)
  { ini_x_ = ini_x; }

  /// Set dual initial solution.
  ///
  /// @note All pointers should stay valid until
  /// loading the model into NLSolver.
  void SetDualWarmstart(NLW2_SparseVector_C ini_y)
  { ini_y_ = ini_y; }

  /// Add suffix, e.g., basis statuses.
  /// @return true iff new suffix added (vs replaced.)
  /// @note SOS constraints can be modeled as suffixes
  ///   for some AMPL solvers.
  bool AddSuffix(NLSuffix suf) {
    assert(
          (size_t)SufSize(suf.kind_) == suf.values_.size());
    return suffixes_.Add(std::move(suf));
  }

  /// Information exported by WriteNL()
  struct PreprocessData {
    /// var permutation
    std::vector<int> vperm_;
    /// var inverse permutation
    std::vector<int> vperm_inv_;
  };

  /// Write to NL file.
  /// Recommended usage via class NLSolver.
  /// @return empty string iff ok.
  std::string WriteNL(const std::string& file_stub,
                      NLW2_NLOptionsBasic_C opts,
                      NLUtils &ut, PreprocessData &pd);

  /// Compute objective value
  double ComputeObjValue(const double* x) const;


  /// Get problem name
  const char* ProbName() const { return prob_name_; }
  /// Get variables
  NLW2_ColData_C ColData() const { return vars_; }
  /// Get var names
  const char *const *ColNames() const { return var_names_; }
  /// Get var name [i]
  const char *ColName(int i) const {
    assert(0<=i && i<NumCols());
    return ColNames() ? ColNames()[i] : "";
  }
  /// Lin con matrix
  NLW2_SparseMatrix_C GetA() const { return A_; }
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
  /// Row name [i]
  const char *RowName(int i) const {
    assert(0<=i && i<NumRows());
    return RowNames() ? RowNames()[i] : "";
  }
  /// Obj sense
  int ObjSense() const { return obj_sense_; }
  /// Obj offset
  double ObjOffset() const { return obj_c0_; }
  /// Obj coefs
  const double *ObjCoefficients() const { return obj_c_; }
  /// Hessian format NLW2_HessianFormat...
  int HessianFormat() const { return Q_format_; }
  /// Hessian matrix
  NLW2_SparseMatrix_C Hessian() const { return Q_; }
  /// Obj name
  const char* ObjName() const { return obj_name_; }

  /// Warm start
  NLW2_SparseVector_C Warmstart() const
  { return ini_x_; }
  /// Dual warm start
  NLW2_SparseVector_C DualWarmstart() const
  { return ini_y_; }
  /// Suffixes
  const NLSuffixSet Suffixes() const { return suffixes_; }

protected:
  int SufSize(int kind) const {
    switch (kind&3) {
    case 0: return NumCols();
    case 1: return NumRows();
    default: return 1;         // objectives / problem
    }
  }

private:
  const char* prob_name_ {"mp::NLModel"};
  NLW2_ColData_C vars_ {};
  const char *const *var_names_ {};
  NLW2_SparseMatrix_C A_ {};
  int num_row_ {};
  const double *row_lb_ {};
  const double *row_ub_ {};
  const char *const *row_names_ {};
  int obj_sense_ {};
  double obj_c0_ {};
  const double *obj_c_ {};
  int Q_format_ {};
  NLW2_SparseMatrix_C Q_ {};
  const char* obj_name_ {"obj[1]"};

  NLW2_SparseVector_C ini_y_ {}, ini_x_{};
  NLSuffixSet suffixes_;
};


/// Returned solution.
struct NLSolution {
  /// Any result obtained from the solver?
  operator bool() const { return solve_result_ > -2; }
  /**
   Solve result.
   If >-2, solver interaction successful. Then:

   - -1      *unknown* - unexpected termination
   - 0- 99   *solved* - optimal solution found
   - 100-199 *solved?* - optimal solution indicated, but error likely
   - 200-299 *infeasible* - constraints cannot be satisfied
   - 300-399 *unbounded* - objective can be improved without limit
   - 400-499 *limit* - stopped by a limit that you set (such as on iterations)
   - 500-999 *failure* - stopped by an error condition in the solver

   @note NLSolution is feasible (not proven optimal)
     if *unbounded* or *limit* and \a x_ populated.

   @note Individual solvers may have more specific values,
     see https://ampl.com/products/solvers/solvers-we-sell/.
  */
  int solve_result_ {-2};   // -2: "not set"
  /// Number of solve_message's initial characters
  /// already printed on the screen
  int nbs_{};
  /// Solve message
  std::string solve_message_;
  /// Objective value.
  /// Only returned by Solve(NLModel).
  /// Otherwise, after ReadSolution(),
  /// should be manually computed, e.g.,
  /// by NLModel::ComputeObjValue().
  double obj_val_ {};
  /// Primals
  std::vector<double> x_;
  /// Duals
  std::vector<double> y_;
  /// Suffixes
  NLSuffixSet suffixes_;
};

}  // namespace mp

#endif // NLMODEL_H

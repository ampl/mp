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

#ifndef NLSOLEASY_H
#define NLSOLEASY_H

#include <string>
#include <vector>
#include <set>
#include <memory>
#include <cassert>

#include "mp/basic-defs-c.h"
#include "mp/nl-utils2.h"

namespace mp {

/// Class NLModel_Easy.
///
/// Intermediate representation for special model types:
/// (MI)LP, (MI)QP.
///
/// All pointers should stay valid until
/// loading the model into NLSOL_Easy.
class NLModel_Easy {
public:
  /// Construct
  NLModel_Easy(const char* probname = nullptr)
    : prob_name_(probname ? probname : "NLModelInstance") { }

  /// Add variables (all at once.)
  void SetCols(NLW2_ColData_C vd) { vars_ = vd; }

  /// Add variable names
  void SetColNames(const char *const *nm) { var_names_=nm; }

  /// Add linear constraints (all at once).
  /// Only rowwise matrix supported.
  void SetRows(
      int nr, const double* rlb, const double* rub,
      NLW2_SparseMatrix_C A)
  { num_row_=nr; row_lb_=rlb; row_ub_=rub; A_=A; }

  /// Add constraint names
  void SetRowNames(const char *const *nm) { row_names_=nm; }

  /// Add linear objective (only single objective supported.)
  /// Sense: NLW2_ObjSenseM....
  /// Coefficients: dense vector.
  void SetLinearObjective(int sense, double c0, const double* c)
  { obj_sense_=sense; obj_c0_=c0; obj_c_=c; }

  /// Add Q for the objective quadratic part 0.5 @ x.T @ Q @ x.
  /// Format: NLW2_HessianFormat...
  void SetHessian(int format, NLW2_SparseMatrix_C Q)
  { Q_format_ = format; Q_ = Q; }

  /// Set obj name
  void SetObjName(const char* nm) { obj_name_=(nm ? nm : ""); }

  /// Information exported by WriteNL()
  struct PreprocessData {
    /// var permutation
    std::vector<int> vperm_;
    /// var inverse permutation
    std::vector<int> vperm_inv_;
  };

  /// Write to NL file.
  /// Recommended usage via class NLSOL_Easy.
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

private:
  const char* prob_name_ {"NLSOL_Easy_model"};
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
};


/// Declare NLSOL
class NLSOL;

/// Declare NLHeader
class NLHeader;


/// Class NLSOL_Easy.
///
/// A wrapper for mp::NLSOL to use AMPL solvers
/// for (MI)QP models.
class NLSOL_Easy {
public:
  /// Construct.
  NLSOL_Easy();
  /// Destruct
  ~NLSOL_Easy();

  /// Set NLUtils [OPTIONAL].
  ///
  /// If not provided, default is used.
  void SetNLUtils(mp::NLUtils* put);

  /// Get NLUtils
  NLUtils* GetNLUtils() const;

  /// Set file stub [OPTIONAL].
  ///
  /// Used for filename base of .nl, .col, row, etc. input files,
  /// as well as .sol output files.
  ///
  /// If not provided, a temporary filename is used;
  /// then, .nl is deleted upon object destruction.
  void SetFileStub(std::string stub);

  /// Retrieve file stub.
  const std::string& GetFileStub() const;

  /// Set NL options [OPTIONAL].
  ///
  /// If not provided, default is used.
  void SetNLOptions(NLW2_NLOptionsBasic_C nlo) { nl_opts_=nlo; }

  /// Get NLOptions
  NLW2_NLOptionsBasic_C GetNLOptions() const { return nl_opts_; }

  /// Get error message.
  /// Nonempty iff error occurred.
  const char* GetErrorMessage() const;

  /// Suffix type
  struct Suffix {
    /// Name
    std::string name_;
    /// Suffix table
    std::string table_;
    /// Kind
    int kind_;
    /// Values. Always double precision.
    std::vector<double> values_;

    /// operator<
    bool operator<(const Suffix& s) const {
      return std::make_pair(name_, kind_)
          < std::make_pair(s.name_, s.kind_);
    }
  };

  /// Suffix set.
  using SuffixSet = std::set<Suffix>;

  /// Solution
  struct Solution {
    /// Any result obtained from the solver?
    operator bool() const { return solve_result_ > -2; }
    /// Solve result
    int solve_result_ {-2};   // "unset"
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
    SuffixSet suffixes_;
  };

  /// Solve model and return result.
  ///
  /// @return Solution object
  ///    (has operator bool() for checking
  ///     if any result was obtained.)
  ///
  /// See LoadModel(), Solve(), ReadSolution()
  /// for details.
  Solution Solve(const NLModel_Easy& mdl,
                 const std::string& solver,
                 const std::string& solver_opts) {
    Solution sol;
    if (LoadModel(mdl)
        && Solve(solver, solver_opts)) {
      sol = ReadSolution();
      if (sol.x_.size())
        sol.obj_val_ = mdl.ComputeObjValue(sol.x_.data());
    }
    return sol;
  }

  /// Write NL and any accompanying files.
  /// NL file name base and some options
  /// can be provided, if non-defaults desired,
  /// via SetFileStub() and SetNLOptions().
  ///
  /// @return true if all ok, otherwise see
  ///   GetErrorMessage().
  bool LoadModel(const NLModel_Easy& mdl);

  /// Solve.
  ///
  /// @param solver: solver executable, such as "gurobi".
  /// @param solver_opts: string of solver options,
  ///   such as "outlev=1 writeprob=model.lp".
  ///
  /// @return true if all ok.
  bool Solve(const std::string& solver,
             const std::string& solver_opts);

  /// Read solution.
  ///
  /// @return Solution object
  ///    (has operator bool() for checking
  ///     if any result was obtained.)
  ///
  /// @note To compute objective value,
  ///   execute NLModel_Easy::ComputeObjValue()
  ///   if x_ available.
  Solution ReadSolution();

private:
  std::unique_ptr<NLSOL> p_nlsol_;
  std::unique_ptr<NLHeader> p_nlheader_;
  NLModel_Easy::PreprocessData pd_;
  NLW2_NLOptionsBasic_C nl_opts_;
};

}  // namespace mp

#endif // NLSOLEASY_H

/*
 Optimization problem

 Copyright (C) 2014 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#ifndef MP_PROBLEM_H_
#define MP_PROBLEM_H_

#include <cstddef>  // for std::size_t
#include <limits>
#include <vector>

#include "mp/expr.h"
#include "mp/suffix.h"

namespace mp {

// An optimization problem.
class Problem : public ExprFactory, public SuffixManager {
 private:
  // A variable.
  struct Var {
    double lb;
    double ub;
    Var(double lb, double ub) : lb(lb), ub(ub) {}
  };
  std::vector<Var> vars_;

  // Packed variable types.
  std::vector<bool> var_types_;

  class LinearExpr {
   private:
    struct Term {
      int var_index;
      double coef;
      Term(int var_index, double coef) : var_index(var_index), coef(coef) {}
    };
    std::vector<Term> terms_;

   public:
    void Reserve(int num_terms) {
      terms_.reserve(num_terms);
    }

    void AddTerm(int var_index, double coef) {
      terms_.push_back(Term(var_index, coef));
    }

    // TODO: accessors
  };

  class LinearExprBuilder {
   private:
    LinearExpr *expr_;

   public:
    explicit LinearExprBuilder(LinearExpr *expr) : expr_(expr) {}

    void AddTerm(int var_index, double coef) {
      expr_->AddTerm(var_index, coef);
    }
  };

  // Packed objective types.
  std::vector<bool> obj_types_;

  // Linear parts of objective expessions.
  std::vector<LinearExpr> linear_objs_;

  // Nonlinear parts of objective expressions.
  // The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_objs_;

  // An algebraic constraint.
  struct AlgebraicCon {
    // Linear part of an algebraic constraint expression.
    // Nonlinear parts are stored in nonlinear_cons_ to avoid overhead
    // for linear problems.
    LinearExpr linear_expr;
    double lb;
    double ub;
    AlgebraicCon(double lb, double ub) : lb(lb), ub(ub) {}
  };
  std::vector<AlgebraicCon> algebraic_cons_;

  // Information about complementarity conditions.
  // compl_vars_[i] > 0 means constraint i complements variable
  // compl_vars_[i] - 1. The array can be empty if there are no
  // complementarity conditions.
  std::vector<unsigned> compl_vars_;

  // Nonlinear parts of algebraic constraint expressions.
  // The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_cons_;

  // Logical constraint expressions.
  std::vector<LogicalExpr> logical_cons_;

  // Linear parts of common expressions.
  std::vector<LinearExpr> linear_exprs_;

  // Nonlinear parts of common expressions.
  std::vector<NumericExpr> nonlinear_exprs_;

  // Initial values for variables.
  std::vector<double> initial_values_;

  // Initial values for dual variables.
  std::vector<double> initial_dual_values_;

  std::vector<Function> funcs_;

  static std::size_t max_index() {
    return static_cast<std::size_t>(std::numeric_limits<int>::max());
  }

 public:
  // Returns the number of variables.
  int num_vars() const { return static_cast<int>(vars_.size()); }

  // Returns the number of objectives.
  int num_objs() const { return static_cast<int>(linear_objs_.size()); }

  // Returns the number of algebraic constraints.
  int num_algebraic_cons() const {
    return static_cast<int>(algebraic_cons_.size());
  }

  // Returns the number of logical constraints.
  int num_logical_cons() const {
    return static_cast<int>(logical_cons_.size());
  }

  // Adds a variable.
  void AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    MP_ASSERT(vars_.size() < max_index(), "too many variables");
    vars_.push_back(Var(lb, ub));
    var_types_.push_back(type != var::CONTINUOUS);
  }

  typedef LinearExprBuilder LinearObjBuilder;

  // Adds an objective.
  // Returns a handler for receiving linear terms in the objective.
  LinearObjBuilder AddObj(obj::Type type, NumericExpr expr,
                          int num_linear_terms = 0);

  typedef LinearExprBuilder LinearConBuilder;

  // Adds an algebraic constraint.
  // Returns a handler for receiving linear terms in the constraint.
  LinearConBuilder AddCon(NumericExpr expr, double lb, double ub,
                          int num_linear_terms = 0);

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr) {
    MP_ASSERT(logical_cons_.size() < max_index(),
              "too many logical constraints");
    logical_cons_.push_back(expr);
  }

  // Begins building a common expression (defined variable).
  // Returns a handler for receiving linear terms in the common expression.
  LinearExprBuilder BeginCommonExpr(int num_linear_terms) {
    linear_exprs_.push_back(LinearExpr());
    LinearExpr &linear = linear_exprs_.back();
    linear.Reserve(num_linear_terms);
    return LinearExprBuilder(&linear);
  }

  // Ends building a common expression.
  void EndCommonExpr(LinearExprBuilder, NumericExpr expr, int) {
    nonlinear_exprs_.push_back(expr);
  }

  // Sets a complementarity relation.
  void SetComplement(int con_index, int var_index, int flags);

  // Sets the initial value for a variable.
  void SetInitialValue(int var_index, double value) {
    MP_ASSERT(0 <= var_index && var_index <= num_vars(), "invalid index");
    if (initial_values_.size() <= var_index) {
      initial_values_.reserve(vars_.capacity());
      initial_values_.resize(num_vars());
    }
    initial_values_[var_index] = value;
  }

  // Sets the initial value for a dual variable.
  void SetInitialDualValue(int con_index, double value) {
    MP_ASSERT(0 <= con_index && con_index <= num_algebraic_cons(),
              "invalid index");
    if (initial_dual_values_.size() <= con_index) {
      initial_dual_values_.reserve(algebraic_cons_.capacity());
      initial_dual_values_.resize(num_algebraic_cons());
    }
    initial_dual_values_[con_index] = value;
  }

  struct ColumnSizeHandler {
    void Add(int) {
      // Ignore column sizes as the constraints are stored row-wise.
    }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    return ColumnSizeHandler();
  }

  typedef Suffix *SuffixPtr;
  typedef mp::SuffixSet SuffixSet;

  class SuffixHandler {
   private:
    Suffix *suffix_;

   public:
    explicit SuffixHandler(Suffix *s) : suffix_(s) {}

    // Sets an integer suffix value.
    void SetValue(int index, int value) {
      suffix_->set_value(index, value);
    }

    // Sets a double suffix value.
    void SetValue(int index, double value) {
      suffix_->set_value(index, value);
    }
  };

  // Adds a suffix.
  // name: Suffix name that may not be null-terminated.
  SuffixHandler AddSuffix(fmt::StringRef name, int kind, int num_values);

  // Sets problem information and reserves memory for problem components.
  void SetInfo(const ProblemInfo &info);
};
}  // namespace mp

#endif  // MP_PROBLEM_H_

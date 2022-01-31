/*
 A minimal implementation of the ProblemBuilder concept.

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

#ifndef MP_PROBLEM_BUILDER_H_
#define MP_PROBLEM_BUILDER_H_

#include <cassert>
#include <cstring>
#include <memory>
#include <set>
#include <vector>

#include "mp/common.h"
#include "mp/error.h"
#include "mp/nl-reader.h"
#include "mp/problem.h"

namespace mp {

/// A minimal implementation of the ProblemBuilder concept.
template <typename Impl, typename ExprType>
class ProblemBuilder : public SuffixManager {
 private:
  struct ExprBuilder {
    void AddArg(ExprType arg) { internal::Unused(&arg); }
  };

  template <typename T>
  struct SuffixHandler {
    void SetValue(int index, T value) { internal::Unused(index, value); }
  };

 public:
  typedef ExprType Expr;
  typedef Expr NumericExpr;
  typedef Expr LogicalExpr;
  typedef Expr CountExpr;
  typedef Expr Reference;

  static void ReportUnhandledConstruct(fmt::CStringRef name) {
    throw MakeUnsupportedError(name);
  }

  void SetInfo(const ProblemInfo &) {}

  // Adds a variable.
  void AddVar(double lb, double ub, var::Type type) {
    internal::Unused(lb, ub, type);
    MP_DISPATCH(ReportUnhandledConstruct("variable"));
  }

  struct LinearExprBuilder {
    void AddTerm(int var_index, double coef) {
      internal::Unused(var_index, coef);
    }
  };

  typedef LinearExprBuilder LinearObjBuilder;

  // Adds an objective.
  // Returns a builder for the linear part of the objective expression.
  LinearObjBuilder AddObj(
      obj::Type type, NumericExpr expr, int num_linear_terms) {
    internal::Unused(type, &expr, num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("objective"));
    return LinearObjBuilder();
  }

  typedef LinearExprBuilder LinearConBuilder;

  // Adds an algebraic constraint.
  // Returns a builder for the linear part of the constraint expression.
  LinearConBuilder AddCon(double lb, double ub, NumericExpr expr,
                          int num_linear_terms) {
    internal::Unused(lb, ub, &expr, num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("algebraic constraint"));
    return LinearConBuilder();
  }

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr) {
    internal::Unused(&expr);
    MP_DISPATCH(ReportUnhandledConstruct("logical constraint"));
  }

  struct CommonExpr {
    LinearExprBuilder set_linear_expr(int num_linear_terms) const {
      internal::Unused(num_linear_terms);
      return LinearExprBuilder();
    }
    void set_nonlinear_expr(NumericExpr expr) const {
      internal::Unused(&expr);
    }
    void set_position(int position) const {
      internal::Unused(position);
    }
  };

  CommonExpr common_expr(int expr_index) {
    internal::Unused(expr_index);
    MP_DISPATCH(ReportUnhandledConstruct("common expression"));
    return CommonExpr();
  }

  // Adds a common expression (defined variable).
  CommonExpr AddCommonExpr(NumericExpr expr) {
    internal::Unused(&expr);
    MP_DISPATCH(ReportUnhandledConstruct("common expression"));
    return CommonExpr();
  }

  // Sets a complementarity relation.
  void SetComplementarity(int con_index, int var_index, ComplInfo info) {
    internal::Unused(con_index, var_index, &info);
    MP_DISPATCH(ReportUnhandledConstruct("complementarity constraint"));
  }

  struct Variable {
    void set_value(double value) {
      internal::Unused(value);
      // Initial values are ignored by default.
    }
  };

  Variable var(int index) {
    internal::Unused(index);
    return Variable();
  }

  struct AlgebraicCon {
    void set_dual(double value) {
      internal::Unused(value);
      // Initial dual values are ignored by default.
    }
  };

  AlgebraicCon algebraic_con(int index) {
    internal::Unused(index);
    return AlgebraicCon();
  }

  class Function {
   private:
    // Safe bool type.
    typedef void (Function::*SafeBool)() const;

   public:
    // Returns a value convertible to bool that can be used in conditions but
    // not in comparisons and evaluates to "true" if this function is not null
    // and "false" otherwise.
    // Example:
    //   if (f) {
    //     // Do something if f is not null.
    //   }
    operator SafeBool() const { return 0; }
  };

  void AddFunctions(int num_funcs) {
    function(num_funcs);
  }

  Function DefineFunction(int index, fmt::StringRef name,
                          int num_args, func::Type type) {
    internal::Unused(&name, num_args, type);
    return function(index);
  }

  // Adds a function.
  Function AddFunction(fmt::StringRef name, int num_args, func::Type type) {
    internal::Unused(&name, num_args, type);
    return function(0);
  }

  Function function(int index) {
    internal::Unused(index);
    MP_DISPATCH(ReportUnhandledConstruct("function"));
    return Function();
  }

  typedef SuffixHandler<int> IntSuffixHandler;

  // Adds a suffix.
  IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int num_values) {
    internal::Unused(&name, kind, num_values);
    return IntSuffixHandler();
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  // Adds a suffix.
  DblSuffixHandler AddDblSuffix(fmt::StringRef name, int kind, int num_values) {
    internal::Unused(&name, kind, num_values);
    return DblSuffixHandler();
  }

  typedef ExprBuilder NumericExprBuilder;
  typedef ExprBuilder IteratedExprBuilder;
  typedef ExprBuilder CallExprBuilder;
  typedef ExprBuilder NumberOfExprBuilder;
  typedef ExprBuilder CountExprBuilder;
  typedef ExprBuilder IteratedLogicalExprBuilder;

  NumericExpr MakeNumericConstant(double value) {
    internal::Unused(value);
    MP_DISPATCH(ReportUnhandledConstruct(
                  "numeric constant in nonlinear expression"));
    return NumericExpr();
  }

  Reference MakeVariable(int var_index) {
    internal::Unused(var_index);
    MP_DISPATCH(ReportUnhandledConstruct("variable in nonlinear expression"));
    return Reference();
  }

  NumericExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    internal::Unused(kind, &arg);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return NumericExpr();
  }

  NumericExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    internal::Unused(kind, &lhs, &rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return NumericExpr();
  }

  NumericExpr MakeIf(LogicalExpr condition,
      NumericExpr then_expr, NumericExpr else_expr) {
    internal::Unused(&condition, &then_expr, &else_expr);
    MP_DISPATCH(ReportUnhandledConstruct("if expression"));
    return NumericExpr();
  }

  struct PLTermBuilder {
    void AddSlope(double slope) { internal::Unused(slope); }
    void AddBreakpoint(double breakpoint) { internal::Unused(breakpoint); }
  };

  PLTermBuilder BeginPLTerm(int num_breakpoints) {
    internal::Unused(num_breakpoints);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return PLTermBuilder();
  }
  NumericExpr EndPLTerm(PLTermBuilder builder, Reference arg) {
    internal::Unused(&builder, &arg);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return NumericExpr();
  }

  CallExprBuilder BeginCall(Function func, int num_args) {
    internal::Unused(&func, num_args);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return CallExprBuilder();
  }
  NumericExpr EndCall(CallExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return NumericExpr();
  }

  IteratedExprBuilder BeginIterated(expr::Kind kind, int num_args) {
    internal::Unused(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return IteratedExprBuilder();
  }
  NumericExpr EndIterated(IteratedExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("vararg expression"));
    return NumericExpr();
  }

  NumericExprBuilder BeginSum(int num_args) {
    internal::Unused(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericExprBuilder();
  }
  NumericExpr EndSum(NumericExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericExpr();
  }

  NumberOfExprBuilder BeginNumberOf(int num_args, NumericExpr arg0) {
    internal::Unused(num_args, &arg0);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumberOfExprBuilder();
  }
  NumericExpr EndNumberOf(NumberOfExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumericExpr();
  }

  typedef ExprBuilder SymbolicNumberOfExprBuilder;

  SymbolicNumberOfExprBuilder BeginSymbolicNumberOf(int num_args, Expr arg0) {
    internal::Unused(num_args, &arg0);
    MP_DISPATCH(ReportUnhandledConstruct("symbolic numberof expression"));
    return SymbolicNumberOfExprBuilder();
  }
  NumericExpr EndSymbolicNumberOf(SymbolicNumberOfExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("symbolic numberof expression"));
    return NumericExpr();
  }

  CountExprBuilder BeginCount(int num_args) {
    internal::Unused(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
    return CountExprBuilder();
  }
  NumericExpr EndCount(CountExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
    return NumericExpr();
  }

  LogicalExpr MakeLogicalConstant(bool value) {
    internal::Unused(value);
    MP_DISPATCH(ReportUnhandledConstruct("logical constant"));
    return LogicalExpr();
  }

  LogicalExpr MakeNot(LogicalExpr arg) {
    internal::Unused(&arg);
    MP_DISPATCH(ReportUnhandledConstruct("logical not"));
    return LogicalExpr();
  }

  LogicalExpr MakeBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    internal::Unused(kind, &lhs, &rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    internal::Unused(kind, &lhs, &rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    internal::Unused(kind, &lhs, &rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeImplication(
      LogicalExpr condition, LogicalExpr then_expr, LogicalExpr else_expr) {
    internal::Unused(&condition, &then_expr, &else_expr);
    MP_DISPATCH(ReportUnhandledConstruct("implication expression"));
    return LogicalExpr();
  }

  IteratedLogicalExprBuilder BeginIteratedLogical(
      expr::Kind kind, int num_args) {
    internal::Unused(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return IteratedLogicalExprBuilder();
  }
  LogicalExpr EndIteratedLogical(IteratedLogicalExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("iterated logical expression"));
    return LogicalExpr();
  }

  typedef ExprBuilder PairwiseExprBuilder;

  PairwiseExprBuilder BeginPairwise(expr::Kind kind, int num_args) {
    internal::Unused(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return PairwiseExprBuilder();
  }
  LogicalExpr EndPairwise(PairwiseExprBuilder builder) {
    internal::Unused(&builder);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return LogicalExpr();
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  Expr MakeStringLiteral(fmt::StringRef value) {
    internal::Unused(&value);
    MP_DISPATCH(ReportUnhandledConstruct("string literal"));
    return Expr();
  }

  Expr MakeSymbolicIf(LogicalExpr condition, Expr then_expr, Expr else_expr) {
    internal::Unused(&condition, &then_expr, &else_expr);
    MP_DISPATCH(ReportUnhandledConstruct("symbolic if expression"));
    return Expr();
  }
};

/// An optimization problem with a column-wise constraint matrix.
class ColProblem : public Problem {
 private:
  // Column-wise constraint matrix.
  std::vector<int> col_starts_;
  std::vector<int> row_indices_;
  std::vector<double> coefs_;

  friend class ColProblemBuilder;

 public:
  ColProblem() {}

  template <class Solver>
  explicit ColProblem(const Solver &) {}

  int col_start(int col_index) const { return col_starts_[col_index]; }
  int row_index(int elt_index) const { return row_indices_[elt_index]; }
  double value(int elt_index) const { return coefs_[elt_index]; }

  const int *col_starts() const { return col_starts_.data(); }
  const int *row_indices() const { return row_indices_.data(); }
  const double *values() const { return coefs_.data(); }

  // Returns the built problem. This is used for compatibility with the problem
  // builder API.
  ColProblem &problem() { return *this; }
};

// An NL handler that builds a problem with a column-wise constraint matrix.
class ColProblemBuilder : public internal::NLProblemBuilder<ColProblem> {
 private:
  typedef internal::NLProblemBuilder<ColProblem> Base;

 public:
  explicit ColProblemBuilder(ColProblem &p)
    : internal::NLProblemBuilder<ColProblem>(p) {}

  void OnHeader(const NLHeader &h) {
    Base::OnHeader(h);
    ColProblem &problem = builder();
    problem.col_starts_.reserve(h.num_vars + 1);
    problem.col_starts_.resize(2);
    problem.row_indices_.resize(h.num_con_nonzeros);
    problem.coefs_.resize(h.num_con_nonzeros);
  }

  class ColumnSizeHandler {
   private:
    ColProblem *problem_;

    explicit ColumnSizeHandler(ColProblem &p) : problem_(&p) {}

    friend class ColProblemBuilder;

   public:
    void Add(int size) {
      // Convert column size to column offset.
      problem_->col_starts_.push_back(problem_->col_starts_.back() + size);
    }
  };

  ColumnSizeHandler OnColumnSizes() {
    return ColumnSizeHandler(builder());
  }

  class LinearConHandler {
   private:
    ColProblem *problem_;
    int con_index_;

    friend class ColProblemBuilder;

    LinearConHandler(ColProblem &p, int con_index)
      : problem_(&p), con_index_(con_index) {}

   public:
    void AddTerm(int var_index, double coef) {
      int index = problem_->col_starts_[var_index + 1]++;
      problem_->row_indices_[index] = con_index_;
      problem_->coefs_[index] = coef;
    }
  };

  LinearConHandler OnLinearConExpr(int con_index, int) {
    // Pass zero as the number of linear terms as we store constraints
    // column-wise rather than row-wise.
    Base::OnLinearConExpr(con_index, 0);
    return LinearConHandler(builder(), con_index);
  }
};
}  // namespace mp

#endif  // MP_PROBLEM_BUILDER_H_

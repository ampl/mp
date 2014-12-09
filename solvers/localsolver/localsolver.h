/*
 AMPL solver interface to LocalSolver.

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

#ifndef MP_SOLVERS_LOCALSOLVER_H_
#define MP_SOLVERS_LOCALSOLVER_H_

#include <cassert>
#include <climits>
#include <vector>

#include <localsolver.h>

#include "mp/problem-builder.h"
#include "mp/solver.h"

namespace mp {

namespace ls = localsolver;

class LocalSolver;

// Converts int to lsint.
inline ls::lsint AsLSInt(int value) { return value; }

// This class provides methods for building a problem in LocalSolver format.
class LSProblemBuilder :
    public ProblemBuilder<LSProblemBuilder, ls::LSExpression> {
 private:
  // LocalSolver only supports one model per solver.
  ls::LocalSolver solver_;
  ls::LSModel model_;
  int num_objs_;
  int num_cons_;

  std::vector<ls::LSExpression> vars_;
  std::vector<ls::LSExpression> common_exprs_;
  std::vector<double> initial_values_;

  typedef ProblemBuilder<LSProblemBuilder, ls::LSExpression> Base;

  ls::LSExpression Negate(ls::LSExpression arg) {
    return model_.createExpression(ls::O_Sub, AsLSInt(0), arg);
  }

  // Returns true if e is a zero constant.
  static bool IsConst(ls::LSExpression e, double value) {
    return e.isConstant() && e.getDoubleValue() == value;
  }

  static void RequireZero(ls::LSExpression e, const char *func_name) {
    if (!IsConst(e, 0)) {
      Base::ReportUnhandledConstruct(
            fmt::format("{} with nonzero second parameter", func_name));
    }
  }

  struct HyperbolicTerms {
    ls::LSExpression exp_x;
    ls::LSExpression exp_minus_x;
  };

  HyperbolicTerms MakeHyperbolicTerms(ls::LSExpression arg);

  // Makes a binary expression.
  template <typename LHS, typename RHS>
  ls::LSExpression MakeBinary(ls::LSOperator op, LHS lhs, RHS rhs) {
    return model_.createExpression(op, lhs, rhs);
  }

  // Makes an expression representing arg / 2.
  ls::LSExpression Half(ls::LSExpression arg) {
    return MakeBinary(ls::O_Div, arg, AsLSInt(2));
  }

  // Makes an expression representing arg + 1.
  ls::LSExpression Plus1(ls::LSExpression arg) {
    return MakeBinary(ls::O_Sum, arg, AsLSInt(1));
  }

  // Makes an expression representing lhs div rhs.
  template <typename RHS>
  ls::LSExpression IntDiv(ls::LSExpression lhs, RHS rhs) {
    return MakeBinary(ls::O_Div, MakeBinary(
                        ls::O_Sub, lhs, MakeBinary(ls::O_Mod, lhs, rhs)), rhs);
  }

  void CheckBounds(int index, std::size_t ub) {
    assert(0 <= index && static_cast<std::size_t>(index) <= ub);
    MP_UNUSED2(index, ub);
  }

 public:
  explicit LSProblemBuilder(LocalSolver &);

  ls::LocalSolver &solver() { return solver_; }

  int num_vars() const { return vars_.size(); }

  // Returns the number of objectives.
  // The return value may be different from the one returned by
  // model_.getNbObjectives() because a dummy LocalSolver objective is
  // added if the problem doesn't containt objectives.
  int num_objs() const { return num_objs_; }

  // Returns the number of constraints.
  // The return value may be different from the one returned by
  // model_.getNbConstraints() because two LocalSolver constaints are
  // added for a single range constraint.
  int num_algebraic_cons() const { return num_cons_; }

  ls::LSExpression *vars() { return &vars_[0]; }

  const double *initial_values() const {
    return !initial_values_.empty() ? &initial_values_[0] : 0;
  }

  void SetInfo(const ProblemInfo &info);

  class LinearExprBuilder {
   private:
    LSProblemBuilder &builder_;
    ls::LSExpression expr_;

   public:
    LinearExprBuilder(LSProblemBuilder &builder,
                      ls::LSExpression &expr, ls::LSExpression sum)
      : builder_(builder), expr_(sum) {
      if (expr != ls::LSExpression()) {
        // Add nonlinear expression.
        sum.addOperand(expr);
        expr = sum;
      } else {
        expr = sum;
      }
    }

    ls::LSExpression expr() { return expr_; }

    void AddTerm(int var_index, double coef) {
      if (coef == 0)
        return;
      ls::LSExpression var = builder_.vars_[var_index];
      expr_.addOperand(coef == 1 ? var : builder_.model_.createExpression(
                                     ls::O_Prod, coef, var));
    }
  };

  void AddVar(double lb, double ub, var::Type type);

  typedef LinearExprBuilder LinearObjBuilder;

  LinearObjBuilder AddObj(obj::Type type, ls::LSExpression expr, int);

  typedef LinearExprBuilder LinearConBuilder;

  // Adds an algebraic constraint.
  LinearConBuilder AddCon(double lb, double ub, ls::LSExpression expr, int);

  // Adds a logical constraint.
  void AddCon(ls::LSExpression expr) {
    model_.addConstraint(expr);
  }

  LinearExprBuilder BeginCommonExpr(int num_terms) {
    ls::LSExpression expr;
    return LinearExprBuilder(
          *this, expr,
          num_terms != 0 ? model_.createExpression(ls::O_Sum) : expr);
  }

  void EndCommonExpr(LinearExprBuilder builder, ls::LSExpression expr, int) {
    ls::LSExpression result = builder.expr();
    if (result == ls::LSExpression())
      result = expr;
    if (expr != ls::LSExpression())
      result.addOperand(expr);
    common_exprs_.push_back(result);
  }

  void SetInitialValue(int var_index, double value) {
    // Store initial values because it is not possible to assign them
    // until the model is closed.
    if (initial_values_.empty())
      initial_values_.resize(vars_.size());
    initial_values_[var_index] = value;
  }

  // Ignore Jacobian column sizes.
  ColumnSizeHandler GetColumnSizeHandler() { return ColumnSizeHandler(); }

  ls::LSExpression MakeNumericConstant(double value) {
    ls::lsint int_value = value;
    return int_value == value ?
          model_.createConstant(int_value) : model_.createConstant(value);
  }

  ls::LSExpression MakeVariable(int var_index) {
    CheckBounds(var_index, vars_.size());
    return vars_[var_index];
  }

  ls::LSExpression MakeCommonExpr(int index) {
    return common_exprs_[index];
  }

  ls::LSExpression MakeUnary(expr::Kind kind, ls::LSExpression arg);
  ls::LSExpression MakeBinary(
      expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs);

  ls::LSExpression MakeIf(ls::LSExpression condition,
      ls::LSExpression true_expr, ls::LSExpression false_expr) {
    if (IsConst(true_expr, 1) && IsConst(false_expr, 0))
      return condition;
    return model_.createExpression(ls::O_If, condition, true_expr, false_expr);
  }

  // LocalSolver doesn't support piecewise-liner terms and arbitrary functions.

  class ExprBuilder {
   private:
    ls::LSExpression expr_;

   public:
    explicit ExprBuilder(ls::LSExpression expr) : expr_(expr) {}

    ls::LSExpression expr() const { return expr_; }

    void AddArg(ls::LSExpression arg) { expr_.addOperand(arg); }
  };

  typedef ExprBuilder NumericExprBuilder;
  typedef ExprBuilder VarArgExprBuilder;

  VarArgExprBuilder BeginVarArg(expr::Kind kind, int num_args);
  ls::LSExpression EndVarArg(VarArgExprBuilder builder) {
    return builder.expr();
  }

  ExprBuilder BeginSum(int) {
    return ExprBuilder(model_.createExpression(ls::O_Sum));
  }
  ls::LSExpression EndSum(ExprBuilder builder) { return builder.expr(); }

  class NumberOfExprBuilder {
   private:
    ls::LSModel model_;
    ls::LSExpression numberof_;
    ls::LSExpression value_;

   public:
    NumberOfExprBuilder(ls::LSModel model, ls::LSExpression value)
      : model_(model), numberof_(model.createExpression(ls::O_Sum)),
        value_(value) {}

    ls::LSExpression numberof() const { return numberof_; }

    void AddArg(ls::LSExpression arg) {
      numberof_.addOperand(model_.createExpression(ls::O_Eq, arg, value_));
    }
  };

  NumberOfExprBuilder BeginNumberOf(int, ls::LSExpression value) {
    return NumberOfExprBuilder(model_, value);
  }
  ls::LSExpression EndNumberOf(NumberOfExprBuilder builder) {
    return builder.numberof();
  }

  typedef ExprBuilder CountExprBuilder;

  ExprBuilder BeginCount(int num_args) { return BeginSum(num_args); }
  NumericExpr EndCount(ExprBuilder builder) { return EndSum(builder); }

  ls::LSExpression MakeLogicalConstant(bool value) {
    return model_.createConstant(AsLSInt(value));
  }

  ls::LSExpression MakeNot(ls::LSExpression arg) {
    return model_.createExpression(ls::O_Not, arg);
  }

  ls::LSExpression MakeBinaryLogical(
      expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs);

  ls::LSExpression MakeRelational(
      expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs);

  ls::LSExpression MakeLogicalCount(
      expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs);

  ls::LSExpression MakeImplication(
      ls::LSExpression condition, ls::LSExpression true_expr,
      ls::LSExpression false_expr) {
    return MakeIf(condition, true_expr, false_expr);
  }

  typedef ExprBuilder IteratedLogicalExprBuilder;

  ExprBuilder BeginIteratedLogical(expr::Kind kind, int num_args);
  LogicalExpr EndIteratedLogical(ExprBuilder builder) { return builder.expr(); }

  struct PairwiseExprBuilder {
    expr::Kind kind;
    std::vector<ls::LSExpression> args;

    PairwiseExprBuilder(expr::Kind k, int num_args) : kind(k) {
      args.reserve(num_args);
    }
    void AddArg(ls::LSExpression arg) { args.push_back(arg); }
  };

  PairwiseExprBuilder BeginPairwise(expr::Kind kind, int num_args) {
    return PairwiseExprBuilder(kind, num_args);
  }

  ls::LSExpression EndPairwise(PairwiseExprBuilder builder);
};

class LocalSolver : public SolverImpl<LSProblemBuilder> {
 private:
  // Integer options.
  enum Option {
    SEED,
    THREADS,
    ANNEALING_LEVEL,
    VERBOSITY,
    TIME_BETWEEN_DISPLAYS,
    TIMELIMIT,
    NUM_OPTIONS
  };

  int options_[NUM_OPTIONS];
  fmt::LongLong iterlimit_;
  std::string logfile_;

  int DoGetIntOption(const SolverOption &, Option id) const {
    return options_[id];
  }

  template <int LB, int UB>
  void DoSetIntOption(const SolverOption &opt, int value, Option id) {
    if (value < LB || value > UB)
      throw InvalidOptionValue(opt, value);
    options_[id] = value;
  }

  void SetNonnegativeIntOption(const SolverOption &opt, int value, Option id) {
    return DoSetIntOption<0, INT_MAX>(opt, value, id);
  }

  std::string GetVerbosity(const SolverOption &opt) const;
  void SetVerbosity(const SolverOption &opt, fmt::StringRef value);

  std::string GetLogFile(const SolverOption &) const { return logfile_; }
  void SetLogFile(const SolverOption &, fmt::StringRef value) {
    logfile_.assign(value.c_str(), value.size());
  }

  fmt::LongLong GetIterLimit(const SolverOption &) const { return iterlimit_; }
  void SetIterLimit(const SolverOption &opt, fmt::LongLong value) {
    if (value < 1)
      throw InvalidOptionValue(opt, value);
    iterlimit_ = value;
  }

 protected:
  virtual void DoSolve(ls::LocalSolver &s) {
    s.solve();
  }

 public:
  LocalSolver();

  void Solve(ProblemBuilder &builder, SolutionHandler &sh);

  LocalSolver &GetProblemBuilder(fmt::StringRef) { return *this; }
};
}  // namespace mp

#endif  // MP_SOLVERS_LOCALSOLVER_H_

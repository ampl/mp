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
 public:
  class Bound {
   private:
    int index_;
    bool isint_;
    union {
      int int_value_;
      double dbl_value_;
    };

   public:
    Bound(int index, int value)
      : index_(index), isint_(true), int_value_(value) {}
    Bound(int index, double value)
      : index_(index), isint_(false), dbl_value_(value) {}

    int index() const { return index_; }

    ls::lsint int_value() const {
      return isint_ ? AsLSInt(int_value_) : static_cast<ls::lsint>(dbl_value_);
    }
    double dbl_value() const {
      return isint_ ? static_cast<double>(int_value_) : dbl_value_;
    }
  };

 private:
  // LocalSolver only supports one model per solver.
  ls::LocalSolver solver_;
  ls::LSModel model_;
  int num_cons_;
  double pl_bigm_;

  std::vector<ls::LSExpression> vars_;
  std::vector<ls::LSExpression> common_exprs_;
  std::vector<double> initial_values_;

  struct ObjInfo {
    obj::Type type;
    ls::LSExpression expr;
    ObjInfo(obj::Type t, ls::LSExpression e) : type(t), expr(e) {}
  };

  std::vector<ObjInfo> objs_;
  std::vector<Bound> obj_bounds_;

  std::vector<LogicalExpr> logical_cons_;

  static const double LS_INF;

  static ls::lsint ConvertToInt(double value) {
    localsolver::lsint int_value = static_cast<localsolver::lsint>(value);
    if (int_value != value)
      throw mp::Error("value {} can't be represented as int", value);
    return int_value;
  }

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
  ls::LSExpression TruncDiv(ls::LSExpression lhs, RHS rhs) {
    return MakeBinary(ls::O_Div, MakeBinary(
                        ls::O_Sub, lhs, MakeBinary(ls::O_Mod, lhs, rhs)), rhs);
  }

  void CheckBounds(int index, std::size_t ub) {
    assert(0 <= index && static_cast<std::size_t>(index) <= ub);
    internal::Unused(index, ub);
  }

  void SetInitialValue(int var_index, double value) {
    // Store initial values because it is not possible to assign them
    // until the model is closed.
    if (initial_values_.empty())
      initial_values_.resize(vars_.size());
    initial_values_[var_index] = value;
  }

  template <typename T>
  class SuffixHandler {
   private:
    std::vector<Bound> *bounds_;

   public:
    explicit SuffixHandler(std::vector<Bound> *bounds = 0) : bounds_(bounds) {}

    void SetValue(int index, T value) {
      if (bounds_)
        bounds_->push_back(Bound(index, value));
    }
  };

  template <typename T>
  SuffixHandler<T> AddSuffix(fmt::StringRef name, suf::Kind kind,
                             int num_values) {
    if (name != "bound" || kind != suf::OBJ)
      return SuffixHandler<T>();
    obj_bounds_.reserve(num_values);
    return SuffixHandler<T>(&obj_bounds_);
  }

 public:
  explicit LSProblemBuilder(LocalSolver &, fmt::StringRef = "");

  ls::LocalSolver &solver() { return solver_; }

  int num_vars() const { return static_cast<int>(vars_.size()); }

  // Returns the number of objectives.
  // The return value may be different from the one returned by
  // model_.getNbObjectives() because a dummy LocalSolver objective is
  // added if the problem doesn't containt objectives.
  int num_objs() const { return static_cast<int>(objs_.size()); }

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
                      ls::LSExpression term, ls::LSExpression sum)
      : builder_(builder), expr_(sum) {
      if (term.getIndex() != -1)
        sum.addOperand(term);
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

  class Variable {
   private:
    LSProblemBuilder *builder_;
    int index_;

    friend class LSProblemBuilder;

    Variable(LSProblemBuilder *b, int index) : builder_(b), index_(index) {}

   public:
    void set_lb(double lb) const {
      ls::LSExpression &var = builder_->vars_[index_];
      double inf = std::numeric_limits<double>::infinity();
      if (var.getOperator() == ls::O_Float) {
        var.setOperand(0, lb == -inf ? -LS_INF : lb);
      } else {
        ls::lsint int_lb = lb == -inf ?
              std::numeric_limits<int>::min() : ConvertToInt(lb);
        var.setOperand(0, int_lb);
      }
    }

    void set_ub(double ub) const {
      ls::LSExpression &var = builder_->vars_[index_];
      double inf = std::numeric_limits<double>::infinity();
      if (var.getOperator() == ls::O_Float) {
        var.setOperand(1, ub == inf ? LS_INF : ub);
      } else {
        ls::lsint int_ub = ub == inf ?
              std::numeric_limits<int>::max() : ConvertToInt(ub);
        var.setOperand(1, int_ub);
      }
    }

    void set_value(double value) const {
      builder_->SetInitialValue(index_, value);
    }
  };

  typedef Variable MutVariable;

  Variable var(int index) {
    CheckBounds(index, vars_.size());
    return Variable(this, index);
  }

  void AddVar(double lb, double ub, var::Type type);

  typedef LinearExprBuilder LinearObjBuilder;

  void AddObj(obj::Type type, ls::LSExpression expr = ls::LSExpression()) {
    objs_.push_back(ObjInfo(type, expr));
  }

  class Objective {
   private:
    LSProblemBuilder *builder_;
    int index_;

    friend class LSProblemBuilder;

    Objective(LSProblemBuilder *b, int index) : builder_(b), index_(index) {}

   public:
    void set_type(obj::Type type) const { builder_->objs_[index_].type = type; }

    void set_nonlinear_expr(ls::LSExpression expr) const {
      builder_->objs_[index_].expr = expr;
    }

    LinearObjBuilder set_linear_expr(int) const {
      ls::LSExpression sum = builder_->model_.createExpression(ls::O_Sum);
      ls::LSExpression &expr = builder_->objs_[index_].expr;
      LinearObjBuilder builder(*builder_, expr, sum);
      expr = sum;
      return builder;
    }
  };

  Objective obj(int index) {
    CheckBounds(index, objs_.size());
    return Objective(this, index);
  }

  class LogicalCon {
   private:
    LSProblemBuilder *builder_;
    int index_;

    friend class LSProblemBuilder;

    LogicalCon(LSProblemBuilder *b, int index) : builder_(b), index_(index) {}

   public:
    void set_expr(LogicalExpr expr) {
      builder_->logical_cons_[index_] = expr;
    }
  };

  LogicalCon logical_con(int index) {
    CheckBounds(index, logical_cons_.size());
    return LogicalCon(this, index);
  }

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
    if (result.getIndex() == -1)
      result = expr;
    else if (expr.getIndex() != -1)
      result.addOperand(expr);
    common_exprs_.push_back(result);
  }

  ls::LSExpression MakeNumericConstant(double value) {
    ls::lsint int_value = static_cast<ls::lsint>(value);
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

  struct PLTermBuilder {
    std::vector<double> slopes;
    std::vector<double> breakpoints;

    explicit PLTermBuilder(int num_breakpoints) {
      slopes.reserve(num_breakpoints + 1);
      breakpoints.reserve(num_breakpoints);
    }

    void AddSlope(double slope) { slopes.push_back(slope); }
    void AddBreakpoint(double breakpoint) { breakpoints.push_back(breakpoint); }
  };

  PLTermBuilder BeginPLTerm(int num_breakpoints) {
    if (localsolver::LSVersion::getMajorVersionNumber() < 5)
      Base::BeginPLTerm(num_breakpoints);
    return PLTermBuilder(num_breakpoints);
  }
  ls::LSExpression EndPLTerm(
        const PLTermBuilder &builder, ls::LSExpression arg);

  class ExprBuilder {
   private:
    ls::LSExpression expr_;

   public:
    explicit ExprBuilder(ls::LSExpression expr) : expr_(expr) {}

    ls::LSExpression expr() const { return expr_; }

    void AddArg(ls::LSExpression arg) { expr_.addOperand(arg); }
  };

  typedef ExprBuilder NumericExprBuilder;
  typedef ExprBuilder IteratedExprBuilder;

  ExprBuilder BeginIterated(expr::Kind kind, int num_args);
  ls::LSExpression EndIterated(ExprBuilder builder) {
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

  const std::vector<Bound> &obj_bounds() const { return obj_bounds_; }

  typedef SuffixHandler<int> IntSuffixHandler;
  typedef SuffixHandler<double> DblSuffixHandler;

  IntSuffixHandler AddIntSuffix(fmt::StringRef name, suf::Kind kind,
                                int num_values) {
    return AddSuffix<int>(name, kind, num_values);
  }

  DblSuffixHandler AddDblSuffix(fmt::StringRef name, suf::Kind kind,
                                int num_values) {
    return AddSuffix<double>(name, kind, num_values);
  }

  // Returns the built problem.
  LSProblemBuilder &problem() {
    for (std::vector<ObjInfo>::const_iterator
         i = objs_.begin(), end = objs_.end(); i != end; ++i) {
      ls::LSObjectiveDirection dir =
          i->type == obj::MIN ? ls::OD_Minimize : ls::OD_Maximize;
      model_.addObjective(i->expr, dir);
    }
    for (std::vector<ls::LSExpression>::const_iterator
         i = logical_cons_.begin(), end = logical_cons_.end(); i != end; ++i) {
      model_.addConstraint(*i);
    }
    return *this;
  }
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
  double pl_bigm_;
  std::string logfile_;
  std::string envfile_;
  mp::OptionValueInfo verbosities_[4];

  struct OptionInfo {
    Option opt;
    int lb;
    OptionInfo(Option opt, int lb) : opt(opt), lb(lb) {}
  };

  int DoGetIntOption(const SolverOption &, OptionInfo info) const {
    return options_[info.opt];
  }

  template <int UB>
  void DoSetIntOption(const SolverOption &opt, int value, OptionInfo info) {
    if (value < info.lb || value > UB)
      throw InvalidOptionValue(opt, value);
    options_[info.opt] = value;
  }

  std::string GetVerbosity(const SolverOption &opt) const;
  void SetVerbosity(const SolverOption &opt, fmt::StringRef value);

  std::string GetLogFile(const SolverOption &) const { return logfile_; }
  void SetLogFile(const SolverOption &, fmt::StringRef value) {
    logfile_.assign(value.data(), value.size());
  }

  std::string GetEnvFile(const SolverOption &) const { return envfile_; }
  void SetEnvFile(const SolverOption &, fmt::StringRef value) {
    envfile_.assign(value.data(), value.size());
  }

  fmt::LongLong GetIterLimit(const SolverOption &) const { return iterlimit_; }
  void SetIterLimit(const SolverOption &opt, fmt::LongLong value) {
    int lb = localsolver::LSVersion::getMajorVersionNumber() < 5 ? 1 : 0;
    if (value < lb)
      throw InvalidOptionValue(opt, value);
    iterlimit_ = value;
  }

  double GetPLBigM(const SolverOption &) const { return pl_bigm_; }
  void SetPLBigM(const SolverOption &, double value) { pl_bigm_ = value; }

 protected:
  virtual void DoSolve(ls::LocalSolver &s) {
    s.solve();
  }

 public:
  LocalSolver();

  void Solve(ProblemBuilder &builder, SolutionHandler &sh);

  double pl_bigm() const { return pl_bigm_; }
};
}  // namespace mp

#endif  // MP_SOLVERS_LOCALSOLVER_H_

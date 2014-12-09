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
#include "mp/suffix.h"

namespace mp {

// A minimal implementation of the ProblemBuilder concept.
template <typename Impl, typename ExprType>
class ProblemBuilder : public SuffixManager {
 private:
  struct ExprBuilder {
    void AddArg(ExprType arg) { MP_UNUSED(arg); }
  };

  template <typename T>
  struct SuffixHandler {
    void SetValue(int index, T value) { MP_UNUSED2(index, value); }
  };

 public:
  typedef ExprType Expr;
  typedef Expr NumericExpr;
  typedef Expr LogicalExpr;
  typedef Expr CountExpr;
  typedef Expr Variable;

  static void ReportUnhandledConstruct(fmt::StringRef name) {
    throw MakeUnsupportedError(name);
  }

  void SetInfo(const ProblemInfo &) {}
  void EndBuild() {}

  // Adds a variable.
  void AddVar(double lb, double ub, var::Type type) {
    MP_UNUSED3(lb, ub, type);
    MP_DISPATCH(ReportUnhandledConstruct("variable"));
  }

  struct LinearExprBuilder {
    void AddTerm(int var_index, double coef) { MP_UNUSED2(var_index, coef); }
  };

  typedef LinearExprBuilder LinearObjBuilder;

  // Adds an objective.
  // Returns a builder for the linear part of the objective expression.
  LinearObjBuilder AddObj(
      obj::Type type, NumericExpr expr, int num_linear_terms) {
    MP_UNUSED3(type, expr, num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("objective"));
    return LinearObjBuilder();
  }

  typedef LinearExprBuilder LinearConBuilder;

  // Adds an algebraic constraint.
  // Returns a builder for the linear part of the constraint expression.
  LinearConBuilder AddCon(double lb, double ub, NumericExpr expr,
                          int num_linear_terms) {
    MP_UNUSED3(expr, lb, ub); MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("algebraic constraint"));
    return LinearConBuilder();
  }

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr) {
    MP_UNUSED(expr);
    MP_DISPATCH(ReportUnhandledConstruct("logical constraint"));
  }

  // Adds a common expression (defined variable).
  // Returns a builder for the linear part of the common expression.
  LinearExprBuilder BeginCommonExpr(int num_linear_terms) {
    MP_UNUSED(num_linear_terms);
    MP_DISPATCH(ReportUnhandledConstruct("common expression"));
    return LinearExprBuilder();
  }

  void EndCommonExpr(LinearExprBuilder builder,
                     NumericExpr expr, int position) {
    MP_UNUSED3(builder, expr, position);
  }

  // Sets a complementarity relation.
  void SetComplement(int con_index, int var_index, int flags) {
    MP_UNUSED3(con_index, var_index, flags);
    MP_DISPATCH(ReportUnhandledConstruct("complementarity constraint"));
  }

  void SetInitialValue(int var_index, double value) {
    MP_UNUSED2(var_index, value);
    // Initial values are ignored by default.
  }
  void SetInitialDualValue(int con_index, double value) {
    MP_UNUSED2(con_index, value);
    // Initial dual values are ignored by default.
  }

  struct ColumnSizeHandler {
    void Add(int size) { MP_UNUSED(size); }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    MP_DISPATCH(ReportUnhandledConstruct("Jacobian column size"));
    return ColumnSizeHandler();
  }

  struct Function {};

  // Adds a function.
  Function AddFunction(fmt::StringRef name, int num_args, func::Type type) {
    MP_UNUSED3(name, num_args, type);
    MP_DISPATCH(ReportUnhandledConstruct("function"));
    return Function();
  }

  typedef SuffixHandler<int> IntSuffixHandler;

  // Adds a suffix.
  IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int num_values) {
    MP_UNUSED3(kind, num_values, name);
    MP_DISPATCH(ReportUnhandledConstruct("integer suffix"));
    return IntSuffixHandler();
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  // Adds a suffix.
  DblSuffixHandler AddDblSuffix(fmt::StringRef name, int kind, int num_values) {
    MP_UNUSED3(kind, num_values, name);
    MP_DISPATCH(ReportUnhandledConstruct("double suffix"));
    return DblSuffixHandler();
  }

  typedef ExprBuilder NumericExprBuilder;
  typedef ExprBuilder VarArgExprBuilder;
  typedef ExprBuilder CallExprBuilder;
  typedef ExprBuilder NumberOfExprBuilder;
  typedef ExprBuilder CountExprBuilder;
  typedef ExprBuilder IteratedLogicalExprBuilder;

  NumericExpr MakeNumericConstant(double value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct(
                  "numeric constant in nonlinear expression"));
    return NumericExpr();
  }

  Variable MakeVariable(int var_index) {
    MP_UNUSED(var_index);
    MP_DISPATCH(ReportUnhandledConstruct("variable in nonlinear expression"));
    return Variable();
  }

  NumericExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    MP_UNUSED2(kind, arg);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return NumericExpr();
  }

  NumericExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return NumericExpr();
  }

  NumericExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("if expression"));
    return NumericExpr();
  }

  struct PLTermBuilder {
    void AddSlope(double slope) { MP_UNUSED(slope); }
    void AddBreakpoint(double breakpoint) { MP_UNUSED(breakpoint); }
  };

  PLTermBuilder BeginPLTerm(int num_breakpoints) {
    MP_UNUSED(num_breakpoints);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return PLTermBuilder();
  }
  NumericExpr EndPLTerm(PLTermBuilder builder, Variable var) {
    MP_UNUSED2(builder, var);
    MP_DISPATCH(ReportUnhandledConstruct("piecewise-linear term"));
    return NumericExpr();
  }

  CallExprBuilder BeginCall(Function func, int num_args) {
    MP_UNUSED2(func, num_args);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return CallExprBuilder();
  }
  NumericExpr EndCall(CallExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("function call"));
    return NumericExpr();
  }

  VarArgExprBuilder BeginVarArg(expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return VarArgExprBuilder();
  }
  NumericExpr EndVarArg(VarArgExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("vararg expression"));
    return NumericExpr();
  }

  NumericExprBuilder BeginSum(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericExprBuilder();
  }
  NumericExpr EndSum(NumericExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("sum"));
    return NumericExpr();
  }

  CountExprBuilder BeginCount(int num_args) {
    MP_UNUSED(num_args);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
    return CountExprBuilder();
  }
  NumericExpr EndCount(CountExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("count expression"));
    return NumericExpr();
  }

  NumberOfExprBuilder BeginNumberOf(int num_args, NumericExpr value) {
    MP_UNUSED2(num_args, value);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumberOfExprBuilder();
  }
  NumericExpr EndNumberOf(NumberOfExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("numberof expression"));
    return NumericExpr();
  }

  LogicalExpr MakeLogicalConstant(bool value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("logical constant"));
    return LogicalExpr();
  }

  LogicalExpr MakeNot(LogicalExpr arg) {
    MP_UNUSED(arg);
    MP_DISPATCH(ReportUnhandledConstruct("logical not"));
    return LogicalExpr();
  }

  LogicalExpr MakeBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return LogicalExpr();
  }

  LogicalExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    MP_DISPATCH(ReportUnhandledConstruct("implication expression"));
    return LogicalExpr();
  }

  IteratedLogicalExprBuilder BeginIteratedLogical(
      expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct(str(kind)));
    return IteratedLogicalExprBuilder();
  }
  LogicalExpr EndIteratedLogical(IteratedLogicalExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("iterated logical expression"));
    return LogicalExpr();
  }

  typedef ExprBuilder PairwiseExprBuilder;

  PairwiseExprBuilder BeginPairwise(expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return PairwiseExprBuilder();
  }
  LogicalExpr EndPairwise(PairwiseExprBuilder builder) {
    MP_UNUSED(builder);
    MP_DISPATCH(ReportUnhandledConstruct("alldiff expression"));
    return LogicalExpr();
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  Expr MakeStringLiteral(fmt::StringRef value) {
    MP_UNUSED(value);
    MP_DISPATCH(ReportUnhandledConstruct("string literal"));
    return Expr();
  }
};
}  // namespace mp

#endif  // MP_PROBLEM_BUILDER_H_
